# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import logging
from collections import defaultdict

from typing import Any, Dict, List, Optional, Set, Tuple

from antismash.config import ConfigType, get_config

from antismash.common import path
from antismash.common.hmmer import run_hmmer_copy, HmmerResults
from antismash.common.module_results import ModuleResults
from antismash.common.secmet import Record, Region
from antismash.common.secmet.features import CDSFeature, FeatureLocation, RRE
from antismash.common.secmet.locations import location_from_string

class RREFinderResults(ModuleResults):
    """ Results class for the RREFinder analysis"""
    schema_version = 1  # when the data format in the results changes, this needs to be incremented

    def __init__(self, record_id: str, bitscore_cutoff: float, min_length: int, hits_per_protocluster: Dict[int, List[str]],
                 hit_info: Dict[str, List[Dict[str, Any]]]) -> None:
        super().__init__(record_id)
        # The cutoffs used for hmmscan
        self.bitscore_cutoff = bitscore_cutoff
        self.min_length = min_length
        # All the locus tags that are an RRE hit per protocluster
        self.hits_per_protocluster = hits_per_protocluster
        # All the hit info per locus tag
        self.hit_info = hit_info
        self.features = [] # type: List[Feature] # features created for RREs
        self.tool = 'rrefinder'
        self.database = 'RREFam.hmm'
        self.detection = 'hmmscan'
        self.convert_hits_to_features()

    def convert_hits_to_features(self) -> None:
        '''Convert all the hits found to features'''
        domain_counts = defaultdict(int) # type: Dict[str, int]
        for locus_tag, hits in self.hit_info.items():
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
            "hit_info": self.hit_info,
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
        filtered_hit_info, filtered_hits_per_protocluster = filter_hits(json['hit_info'], json['hits_per_protocluster'], min_length, bitscore_cutoff)
        RRE_results = RREFinderResults(record.id, bitscore_cutoff, min_length, filtered_hits_per_protocluster, filtered_hit_info)
        return RRE_results

def is_ripp(product: str) -> bool:
    ripp_products = ['bacteriocin','cyanobactin','lanthipeptide',
                     'lassopeptide','linaridin','thiopeptide','sactipeptide',
                      'proteusin','glycocin','bottromycin','microcin']
    return product in ripp_products
    
def gather_rre_candidates(record: Record) -> Tuple[Dict[int, List[str]], Dict[str, CDSFeature]]:
    '''Gather all RRE candidates that need to be analyzed with hmmscan
       and all unique candidates (by CDS name) to prevent double analysis
       and features in the case of overlapping RiPP protoclusters.
    '''
    rre_candidates_per_protocluster = defaultdict(list) # type: Dict[int, List[str]]
    cds_info = defaultdict(CDSFeature) # type: Dict[str, CDSFeature]
    
    for region in record.get_regions():
        for protocluster in region.get_unique_protoclusters():
            if is_ripp(protocluster.product):
                protocluster_number = protocluster.get_protocluster_number()
                for cds in protocluster.cds_children:
                    cds_name = cds.get_name()
                    rre_candidates_per_protocluster[protocluster_number].append(cds_name)
                    cds_info[cds_name] = cds
    return rre_candidates_per_protocluster, cds_info

def extract_rre_hits(hmm_result: HmmerResults) -> Dict[str, List[Dict[str, Any]]]:
    '''Extract the hits per locus_tag from a HmmerResults object'''
    hit_info = defaultdict(list) # type: Dict[str, List[Dict[str, Any]]]
    for hit in hmm_result.hits:
        hit_info[hit['locus_tag']].append(hit)
    return hit_info

def filter_hits(hit_info: Dict[str, List[Dict[str, Any]]], candidates_per_protocluster: Dict[int, List[str]],
                min_length: int, bitscore_cutoff: float) -> Tuple[Dict[str, List[Dict[str, Any]]], Dict[int, List[str]]]:
    '''Filter the hits based on the bitscore and length criteria'''
    filtered_tags_per_protocluster = defaultdict(list) # type: Dict[List[str]]
    filtered_hit_info = defaultdict(list) # type: Dict[str, List[Dict[str, Any]]]
    for protocluster, locus_tags in candidates_per_protocluster.items():
        for locus_tag in locus_tags:
            hits = hit_info[locus_tag]
            for hit in hits:
                if check_hmm_hit(hit, min_length, bitscore_cutoff):
                    filtered_hit_info[locus_tag].append(hit)
            if filtered_hit_info.get(locus_tag):
                filtered_tags_per_protocluster[protocluster].append(locus_tag)
    return filtered_hit_info, filtered_tags_per_protocluster

def check_hmm_hit(hit: Dict[str, Any], min_length: int, bitscore_cutoff: float) -> Dict[str, Any]:
    return (hit['protein_end'] - hit['protein_start']) >= min_length and (hit['score'] >= bitscore_cutoff)
    # Locations come from BioPython's HSPs, so they are pythonic (zero-based, half-open)

def run_rrefinder(record: Record, bitscore_cutoff: float, min_length: int) -> RREFinderResults:
    # Gather all RRE candidates
    candidates_per_protocluster, cds_info = gather_rre_candidates(record)
    # Run hmmscan per protocluster and gather the hits
    hmm_database = path.get_full_path(__file__, 'data', 'RREFam.hmm')
    hmm_results = run_hmmer_copy(record, cds_info.values(), max_evalue=1, min_score=bitscore_cutoff, database=hmm_database, tool='rrefinder')
    # Extract the RRE hits
    hit_info = extract_rre_hits(hmm_results)
    # Filter the hits
    filtered_hit_info, filtered_hits_per_protocluster = filter_hits(hit_info, candidates_per_protocluster, min_length, bitscore_cutoff)
    # Convert to RREFinderResults object
    RRE_results = RREFinderResults(record.id, bitscore_cutoff, min_length, filtered_hits_per_protocluster, filtered_hit_info)
    return RRE_results
