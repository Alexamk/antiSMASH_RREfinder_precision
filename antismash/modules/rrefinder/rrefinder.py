import logging
from collections import defaultdict

from typing import Any, Dict, List, Optional

from antismash.common.secmet import Record
from antismash.common.module_results import ModuleResults
from antismash.config import ConfigType, get_config

from antismash.common.hmmer import run_hmmer_copy, HmmerResults
from antismash.common import path

from antismash.common.secmet import Region
from antismash.common.secmet.features import Domain, Feature, FeatureLocation

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

    def add_to_record(self, record: Record) -> None:
        """ Adds the analysis results to the record """
        if record.id != self.record_id:
            raise ValueError("Record to store in and record analysed don't match")
        raise NotImplementedError()  # remove this when completed

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

        # Check if the cutoff options set in this run are the same or stricter than in the previous run
        if bitscore_cutoff < prev_bitscore_cutoff:
            logging.debug("RREFinderResults bitscore cutoff has changed, discarding previous results")
            return None
        elif min_length < prev_min_length:
            logging.debug("RREFinderResults minimum length has changed, discarding previous results")
            return None

        # Refilter the hits (in case the cutoff is now more stringent)
        filtered_hits = filter_hits(json['hits_per_protocluster'], min_length, bitscore_cutoff)

        return RREFinderResults(record.id, bitscore_cutoff, min_length, filtered_hits)

def is_ripp(product: str) -> bool:
    ripp_products = ['bacteriocin','cyanobactin','lanthipeptide',
                     'lassopeptide','linaridin','thiopeptide','sactipeptide',
                      'proteusin','glycocin','bottromycin','microcin']
    return product in ripp_products
    
def run_hmmscan_rrefinder(record: Record, bitscore_cutoff: float) -> Dict[int, List[Dict[str, Any]]]:
    max_evalue = 1 # Mandatory argument for hmmscan, but hits are only filtered on score here
    hmm_database = path.get_full_path(__file__, 'data', 'RREFinder.hmm')
    
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
    return (hit['protein_end'] - hit['protein_start'] + 1) >= min_length and (hit['score'] >= bitscore_cutoff)
    # +1 to correct for normal (non-pythonic) indeces, which are used in biopython feature locations
    # double check this to be sure

def run_rrefinder(record: Record, bitscore_cutoff: float, min_length: int) -> RREFinderResults:
    # Run hmmscan per protocluster and gather the hits
    hmm_hits = run_hmmscan_rrefinder(record, bitscore_cutoff)
    # Filter the hits
    filtered_hits = filter_hits(hmm_hits, min_length, bitscore_cutoff)
    # Convert to RREFinderResults object
    RRE_results = RREFinderResults(record.id, bitscore_cutoff, min_length, filtered_hits)

    return RRE_results
