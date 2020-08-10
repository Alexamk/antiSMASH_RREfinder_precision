import logging
from collections import defaultdict

from typing import Any, Dict, List, Optional

from antismash.common.secmet import Record
from antismash.common.module_results import ModuleResults
from antismash.config import ConfigType, get_config

from antismash.common.hmmer import run_hmmer_copy, HmmerResults
from antismash.common import path

from antismash.common.secmet import Region

class RREFinderResults(ModuleResults):
    """ Results class for the RREFinder analysis"""
    schema_version = 1  # when the data format in the results changes, this needs to be incremented

    # define whatever construction arguments are needed, record_id is required by the superclass
    # it's good to keep any command line option values here to know when they're changed for --reuse-results
    def __init__(self, record_id: str, cutoff: float, min_length: int, hits: Dict[int, Dict[str, List[Dict[str, Any]]]]) -> None:
        super().__init__(record_id)
        # The cutoff used for hmmscan
        self.cutoff = cutoff
        self.min_length = min_length
        # All the hits per locus_tag, per protocluster ID
        self.hits = hits
#        # Keep protoclusters with RRE hits
#        self.clusters = defaultdict(set)  # type: Dict[int, Set[str]]
#        # e.g. self.clusters[cluster_number] = {gene1_locus, gene2_locus}

    # implement a conversion to a JSON-compatible dictionary
    # all elements must one of: str, int, float, list, or a dict of those types (this can recurse)
    def to_json(self) -> Dict[str, Any]:
        """ Constructs a JSON representation of this instance """

        return {
            "schema_version": self.schema_version,
            "cutoff": self.cutoff,
            "hits": self.hits,
#            "protoclusters": {key: list(val) for key, val in self.clusters.items()}
        }

    # once _all_ analysis modules have completed, their results are added with this method
    # adding to the record during the analysis will cause issues
    def add_to_record(self, record: Record) -> None:
        """ Adds the analysis results to the record """
        if record.id != self.record_id:
            raise ValueError("Record to store in and record analysed don't match")
        # any results would be added here
        # for an example of new features, see antismash.modules.tta
        # for an example of qualifiers, see antismash.modules.t2pks
        # any new feature types or qualifiers would be implemented in antismash.common.secmet,
        #   and would need to be able to be converted to and from biopython's SeqFeature without loss
        # raise NotImplementedError()  # remove this when completed

    # implement a conversion from the JSON-compatible data returned by to_json()
    # this allows --results-reuse to avoid running the module again if not neccessary
    @staticmethod
    def from_json(json: Dict[str, Any], record: Record) -> Optional["RREFinderResults"]:
        """ Constructs a new results instance from a JSON format and the
            original record analysed.
        """
        # check that the previous data version is the same as current, if not, discard the results
        if json["schema_version"] != RREFinderResults.schema_version:
            return None

        # as an example, checking if the example cutoff option matches that of the previous run
        options = get_config()
        if options.template_cutoff != json["cutoff"]:
            # it's nice to log some decisions to the debug logger so that it's easier to follow
            logging.debug("RREFinderResults cutoff has changed, discarding previous results")
            return None

        # the exact reconstruction depends on what details are stored
        # to match the conversion to JSON that would be:
        results = RREFinderResults(json["record_id"], json["cutoff"])
        for other in json["other"]:
            results.some_other_information.append(other)

        return results


def is_ripp(product: str) -> bool:
    ripp_products = ['bacteriocin','cyanobactin','lanthipeptide',
                     'lassopeptide','linaridin','thiopeptide','sactipeptide',
                      'proteusin','glycocin','bottromycin','microcin']
    return product in ripp_products
    
def hmmscan_rrefinder(record: Record, bitscore_cutoff: float, min_length: int) -> Dict[int, Dict[str, List[Dict[str, Any]]]]:
    max_evalue = 1 # Mandatory argument for hmmscan, but hits are only filtered on score here
    hmm_database = path.get_full_path(__file__, 'data', 'RREFinder.hmm')
    
    # Analyse each CDS in each protocluster
    # This is redundant in the case of overlapping RiPP protoclusters
    hmm_hits = defaultdict(lambda: defaultdict(list))
    for region in record.get_regions():
        for protocluster in region.get_unique_protoclusters():
            if is_ripp(protocluster.product):
                protocluster_number = protocluster.get_protocluster_number()
                hmm_result = run_hmmer_copy(record, protocluster.cds_children, max_evalue, bitscore_cutoff, min_length, hmm_database, 'rrefinder')
                hmm_hits_filtered = filter_hmm_hits_min_length(hmm_result.hits, min_length)
                if hmm_hits_filtered:
                    for hit in hmm_hits_filtered:
                        hmm_hits[protocluster_number][hit['locus_tag']].append(hit)
    return hmm_hits
    
def filter_hmm_hits_min_length(hits: List[Dict[str, Any]], min_length: int) -> List[Dict[str, Any]]:
    filtered_hits = []
    for hit in hits:
        length = hit['protein_end'] - hit['protein_start']
        if length >= min_length:
            filtered_hits.append(hit)
    return filtered_hits
    
#def organize_hits_per_locustag(hits_per_protocluster: List[Dict[str, Any]]) -> Dict[str, Dict[str, Any]]):
#    hits_per_locus_tag = {}
#    for hit in hits_

#def collect_hmm_hits(hmm_results: List[HmmerResults]) -> List[Dict[str, Any]]:
#    hits = []
#    for hmm_result in hmm_results:
#        hits.extend(hmm_result.hits)
#    return hits

def run_rrefinder(record: Record, bitscore_cutoff: float) -> RREFinderResults:
    
#    nr_ripps = count_ripps(record)
#    logging.critical('Found %i RiPP protocluster(s).' %nr_ripps)
#    results = RREFinderResults(record.id, bitscore_cutoff, nr_ripps)
    min_length = 50 #TODO: Add this as a commandline argument for antiSMASH
    # Run hmmscan per protocluster
    hmm_hits = hmmscan_rrefinder(record,bitscore_cutoff,min_length)
    # Gather all the hits
    RRE_results = RREFinderResults(record.id,bitscore_cutoff,min_length,hmm_hits)    
    # Add the hits for each protocluster
    logging.critical('Found hits in %i protoclusters.' %len(RRE_results.hits))
    
    return RRE_results
    
    
    
