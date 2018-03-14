# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# initial test if approach even works... later, maybe check for Pfam ids of interest first and only work with these?

#from antismash.common.secmet.feature import Feature, FeatureLocation
import logging
from collections import defaultdict
from typing import Any, Dict, List

#gos = defaultdict(lambda: 0)
from antismash.common import path
from antismash.common.module_results import ModuleResults
from antismash.common.secmet.record import Record


class GeneOntology:
    def __init__(self, _id: str, description: str):
        # assert _id.startswith('GO:')
        if not _id.startswith('GO:'):
            raise ValueError('Invalid Gene Ontology ID')
        self.id = _id
        assert isinstance(description, str)
        self.description = description
    def __str__(self):
        return self.id

class GeneOntologies:
    def __init__(self, pfam: str, gos: List[GeneOntology]):
        self.pfam = str(pfam)
        assert isinstance(gos, list) and gos
        self.go_entries = gos

    def __str__(self) -> str:
        return str([str(go_entry) for go_entry in self.go_entries])


class Pfam2GoResults(ModuleResults):
    """Holds results for Pfam to Gene Ontology module."""
    def __init__(self, record_id: str, pfam_domains_with_gos):
        super().__init__(record_id)
        #  store mapping of PFAM domain ID and GO terms
        self.pfam_domains_with_gos = pfam_domains_with_gos

    def add_to_record(self, record: Record):
        print(self.pfam_domains_with_gos)

    def to_json(self) -> Dict[str, Any]:
        """ Construct a JSON representation of this instance """

    @staticmethod
    def from_json(json: Dict[str, Any], record) -> "Pfam2GoResults":
        """ Constructs a new Pfam2GoResults instance from a json format and the
            original record analysed.
        """


def parse_all_mappings(mapfile):
    go_ids_by_pfam = defaultdict(list)
    go_desc_by_id = {}
    with open(path.get_full_path(__file__, mapfile), 'r') as pfam_map:
        for line in pfam_map:
            if line.startswith('!'):
                continue
            pfam_infos = line.split(' > ')[0]
            go_infos = line.split(' > ')[1].strip()
            pfam_id = pfam_infos.split(' ')[0].replace('Pfam:','')
            go_id = go_infos.split(' ; ')[1]
            #gos[go_id] += 1
            go_readable = go_infos.split(' ; ')[0].replace('GO:','')
            go_ids_by_pfam[pfam_id].append(go_id)
            go_desc_by_id[go_id] = go_readable
    return go_ids_by_pfam,go_desc_by_id


def build_as_i_go(mapfile) -> Dict[str, GeneOntologies]:
    results = {}
    gene_ontology_per_pfam = defaultdict(list)
    with open(path.get_full_path(__file__, mapfile), 'r') as pfam_map:
        # replace with function?
        for line in pfam_map:
            if line.startswith('!'):
                continue
            pfam_info = line.split(' > ')[0]
            go_info = line.split(' > ')[1].strip()
            pfam_id = pfam_info.split(' ')[0].replace('Pfam:','')
            go_id = go_info.split(' ; ')[1]
            go_readable = go_info.split(' ; ')[0].replace('GO:','')
            gene_ontology_per_pfam[pfam_id].append(GeneOntology(go_id, go_readable))
            # if pfam_id not in results:
            #     current_ontologies = GeneOntologies(pfam_id,[GeneOntology(go_id,go_readable)])
            # else:
            #     current_ontologies.go_entries.append(GeneOntology(go_id,go_readable))
            # results[pfam_id] = current_ontologies
    for pfam, ontology_list in gene_ontology_per_pfam.items():
        results[pfam] = GeneOntologies(pfam,ontology_list)
    return results


def construct_gene_ontologies(pfam: str, go_ids: List, go_desc_by_id: Dict) -> GeneOntologies:
    gene_ontologies = []
    for go_id in go_ids:
        if go_id not in go_desc_by_id:
            raise ValueError('Gene Ontology ID has no associated description: %s.' % go_id)
        gene_ontology = GeneOntology(go_id, go_desc_by_id[go_id])
        gene_ontologies.append(gene_ontology)
    return GeneOntologies(pfam,gene_ontologies)


def build_at_the_end(go_ids_by_pfam: Dict, go_desc_by_id: Dict) -> Dict[str, GeneOntologies]:
    results = {}
    #go_ids_by_pfam, go_desc_by_id = parse_all_mappings(mapfile)
    for pfam in go_ids_by_pfam:
        construction = construct_gene_ontologies(pfam, go_ids_by_pfam[pfam], go_desc_by_id)
        results[pfam] = construction

    return results


def get_gos_for_pfams(record) -> Dict[str, List[GeneOntologies]]:
    pfam_domains_with_gos = defaultdict(list)
    pfams = record.get_pfam_domains()
    full_gomap_as_ontologies = build_as_i_go('pfam2go-march-2018.txt')
    if not pfams:
        logging.info('No Pfam domains found')
    for pfam in pfams:
        pfam_ids = pfam.db_xref
        if not pfam_ids:
            logging.info('No Pfam ids found')
        for pfam_id in pfam_ids:
            pfam_id = pfam_id.partition('.')[0] # strip out version number; supposedly faster than split
            gene_ontologies_for_pfam = full_gomap_as_ontologies.get(pfam_id, False)
            if gene_ontologies_for_pfam:
                pfam_domains_with_gos[pfam].append(gene_ontologies_for_pfam)
    return pfam_domains_with_gos


if __name__ == '__main__':
    go_ids_for_debug, go_ids_and_descs = parse_all_mappings('pfam2go-march-2018.txt')
    print('PF08494', build_at_the_end(go_ids_for_debug,go_ids_and_descs)['PF08494'])
    print('PF08494', build_as_i_go('pfam2go-march-2018.txt')['PF08494'])
#print('GO:0004559',go_desc_by_id['GO:0004559'])


#print(list(filter(lambda x: x[1] !=1, gos.items())))