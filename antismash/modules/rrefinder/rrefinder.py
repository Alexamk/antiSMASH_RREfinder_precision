import logging
from collections import defaultdict

from typing import Any, Dict, List, Optional, Set

from antismash.common.secmet import Record
from antismash.common.module_results import ModuleResults
from antismash.config import ConfigType, get_config

from antismash.common.hmmer import run_hmmer_copy, HmmerResults
from antismash.common import path

from antismash.common.secmet import Region
from antismash.common.secmet.features import Domain, Feature, FeatureLocation, RRE

from antismash.common.secmet.locations import location_from_string

class RREFinderResults(ModuleResults):
    """ Results class for the RREFinder analysis"""
    schema_version = 1  # when the data format in the results changes, this needs to be incremented

    def __init__(self, record_id: str, bitscore_cutoff: float, min_length: int, hits_per_protocluster: Dict[int, Dict[str, List[Dict[str, Any]]]]) -> None:
        super().__init__(record_id)
        # The cutoff used for hmmscan
        self.bitscore_cutoff = bitscore_cutoff
        self.min_length = min_length
        # All the hits per locus_tag, per protocluster ID
        self.hits_per_protocluster = hits_per_protocluster
        self.features = [] # type: List[Feature] # features created for RREs
        self.tool = 'rrefinder'
        self.database = 'RREFam.hmm'
        self.detection = 'hmmscan'
        self.convert_hits_to_features()

    def convert_hits_to_features(self) -> None:
        locus_tags_added = set() # type: Set[str]
        domain_counts = defaultdict(int) # type: Dict[str, int]

        for protocluster_number, hits_per_locus_tag in self.hits_per_protocluster.items():
            for locus_tag, hits in hits_per_locus_tag.items():
                if locus_tag in locus_tags_added:
                    continue
                for hit in hits:
                    location = location_from_string(hit['location'])
                    protein_location = FeatureLocation(hit['protein_start'], hit['protein_end'])
                    rre_feature = RRE(location, hit['description'], protein_location, tool=self.tool,
                                      locus_tag=locus_tag, domain=hit['domain'])
                    # Set additional properties
                    for attr in ['score', 'evalue', 'label', 'translation']:
                        setattr(rre_feature, attr, hit[attr])

                    rre_feature.database = self.database
                    rre_feature.detection = self.detection

                    domain_counts[hit['domain']] += 1  # 1-indexed, so increment before use
                    rre_feature.domain_id = "{}_{}_{:04d}".format(self.tool, rre_feature.locus_tag, domain_counts[hit['domain']])

                    self.features.append(rre_feature)
                locus_tags_added.add(locus_tag)

    def add_to_record(self, record: Record) -> None:
        """ Adds the analysis results to the record """
        if record.id != self.record_id:
            raise ValueError("Record to store in and record analysed don't match")

        for feature in self.features:
            record.add_feature(feature)

    def to_json(self) -> Dict[str, Any]:
        """ Constructs a JSON representation of this instance """

        return {
            "schema_version": self.schema_version,
            "bitscore_cutoff": self.bitscore_cutoff,
            "hits_per_protocluster": self.hits_per_protocluster,
            "min_length": self.min_length,
            "record_id": self.record_id,
        }

    @staticmethod
    def from_json(json: Dict[str, Any], record: Record, min_length: int, bitscore_cutoff: float) -> Optional["RREFinderResults"]: # type: ignore  # pylint: disable=arguments-differ
        """ Regenerates the RREFinderResults from json.
            If less strict settings are given (e.g. a lower cutoff
            or a shorter minimum length), the results are discarded.
        """

        # check that the previous data version is the same as current, if not, discard the results
        if json["schema_version"] != RREFinderResults.schema_version:
            return None

        if record.id != json.get("record_id"):
            logging.warning("RREFinder results are for different record, discarding previous results")
            return None

        prev_min_length = json.get("min_length")
        prev_bitscore_cutoff = json.get("bitscore_cutoff")

        if prev_min_length == None or prev_bitscore_cutoff == None:
            raise ValueError('Invalid RREfinderResults json dictionary')
        assert isinstance(prev_min_length, int) and isinstance(prev_bitscore_cutoff, float)

        # Check if the cutoff options set in this run are the same or stricter than in the previous run
        if bitscore_cutoff < prev_bitscore_cutoff:
            logging.debug("RREFinderResults bitscore cutoff has changed, discarding previous results")
            return None
        elif min_length < prev_min_length:
            logging.debug("RREFinderResults minimum length has changed, discarding previous results")
            return None

        # Refilter the hits (in case the cutoff is now more stringent)
        filtered_hits = filter_hits(json['hits_per_protocluster'], min_length, bitscore_cutoff)
        RRE_results = RREFinderResults(record.id, bitscore_cutoff, min_length, filtered_hits)
        return RRE_results

def is_ripp(product: str) -> bool:
    ripp_products = ['bacteriocin','cyanobactin','lanthipeptide',
                     'lassopeptide','linaridin','thiopeptide','sactipeptide',
                      'proteusin','glycocin','bottromycin','microcin']
    return product in ripp_products
    
def run_hmmscan_rrefinder(record: Record, bitscore_cutoff: float) -> Dict[int, Dict[str, List[Dict[str, Any]]]]:
    max_evalue = 1 # Mandatory argument for hmmscan, but hits are only filtered on score here
    hmm_database = path.get_full_path(__file__, 'data', 'RREFam.hmm')
    
    # Scan each CDS in each protocluster
    # This is somewhat redundant in the case of overlapping RiPP protoclusters,
    # as some locus_tags may be scanned twice

    hmm_hits = defaultdict(lambda: defaultdict(list)) # type: Dict[int, Dict[str, List[Dict[str, Any]]]]
    for region in record.get_regions():
        for protocluster in region.get_unique_protoclusters():
            if is_ripp(protocluster.product):
                protocluster_number = protocluster.get_protocluster_number()
                hmm_result = run_hmmer_copy(record, protocluster.cds_children, max_evalue, bitscore_cutoff, hmm_database, 'rrefinder')
                for hit in hmm_result.hits:
                    hmm_hits[protocluster_number][hit['locus_tag']].append(hit)
    return hmm_hits

def filter_hits(hits: Dict[int, Dict[str, List[Dict[str, Any]]]], min_length: int, bitscore_cutoff: float) -> Dict[int, Dict[str, List[Dict[str, Any]]]]:
    filtered_hits = defaultdict(lambda: defaultdict(list)) # type: Dict[int, Dict[str, List[Dict[str, Any]]]]
    for protocluster, locus_tags in hits.items():
        for locus_tag, hits_locus_tag in locus_tags.items():
            for hit in hits_locus_tag:
                if check_hmm_hit(hit, min_length, bitscore_cutoff):
                    filtered_hits[protocluster][locus_tag].append(hit)
    return filtered_hits

def check_hmm_hit(hit: Dict[str, Any], min_length: int, bitscore_cutoff: float) -> Dict[str, Any]:
    return (hit['protein_end'] - hit['protein_start']) >= min_length and (hit['score'] >= bitscore_cutoff)
    # Locations come from BioPython's HSPs, so they are pythonic (zero-based, half-open)

def run_rrefinder(record: Record, bitscore_cutoff: float, min_length: int) -> RREFinderResults:
    # Run hmmscan per protocluster and gather the hits
    hmm_hits = run_hmmscan_rrefinder(record, bitscore_cutoff)
    # Filter the hits
    filtered_hits = filter_hits(hmm_hits, min_length, bitscore_cutoff)
    # Convert to RREFinderResults object
    RRE_results = RREFinderResults(record.id, bitscore_cutoff, min_length, filtered_hits)
    return RRE_results
