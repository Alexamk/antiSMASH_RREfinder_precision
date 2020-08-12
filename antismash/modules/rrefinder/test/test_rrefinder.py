# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest
from unittest.mock import Mock, patch

from antismash.common.test.helpers import DummyHMMResult, DummyRecord

from antismash.common.hmmer import HmmerResults

from antismash.modules.rrefinder.html_output import will_handle
from antismash.modules.rrefinder.rrefinder import (
    is_ripp,
    check_hmm_hit,
    filter_hits,
    run_hmmscan_rrefinder,
    RREFinderResults,
    run_rrefinder,
    )


class TestRREResults(unittest.TestCase):

    def setUp(self):
        self.hit_a1 = {'locus_tag': 'locus_tag_a', 'location': '1000-2000',
                        'label': 'a_generic_label', 'domain': 'RRE_type_A', 'evalue': 0.1,
                        'protein_start': 1, 'protein_end': 50, 'score': 38.0,
                      }
        self.hit_b1 = {'locus_tag': 'locus_tag_b',
                        'label': 'a_generic_label', 'domain': 'RRE_type_B', 'evalue': 0.2,
                        'protein_start': 1, 'protein_end': 90, 'score': 25.0
                        }
        self.filtered_hits = {1: {'locus_tag_a': [self.hit_a1],
                                    },
                               2: {'locus_tag_b': [self.hit_b1],
                                    },
                              }
        self.min_length = 50
        self.bitscore_cutoff = 25.0

        self.record = Mock()
        self.record.id = self.record_id = 'test_record'
        self.res_object = RREFinderResults(self.record_id, self.bitscore_cutoff, self.min_length, self.filtered_hits)

        self.json_dict = {"schema_version": RREFinderResults.schema_version,
                          "bitscore_cutoff": self.bitscore_cutoff,
                          "hits_per_protocluster": self.filtered_hits,
                          "min_length": self.min_length,
                          "record_id": self.record_id, }

    def test_init(self):
        assert self.res_object.record_id == self.record_id
        assert self.res_object.bitscore_cutoff == self.bitscore_cutoff
        assert self.res_object.min_length == self.min_length
        assert self.res_object.hits_per_protocluster == self.filtered_hits

    def test_to_json(self):
        assert self.json_dict == self.res_object.to_json()

    def test_from_json(self):
        res_object = RREFinderResults.from_json(self.json_dict, self.record, self.min_length, self.bitscore_cutoff)
        assert res_object
        assert res_object.hits_per_protocluster == self.filtered_hits
        assert res_object.schema_version == RREFinderResults.schema_version
        assert res_object.bitscore_cutoff == self.bitscore_cutoff
        assert res_object.min_length == self.min_length
        assert res_object.record_id == self.record_id

        json_dict2 = dict(self.json_dict)
        json_dict2['record_id'] = 'another_record'
        assert not RREFinderResults.from_json(json_dict2, self.record, self.min_length, self.bitscore_cutoff)

        json_dict3 = dict(self.json_dict)
        json_dict3['schema_version'] = 'not_a_valid_schema'
        assert not RREFinderResults.from_json(json_dict3, self.record, self.min_length, self.bitscore_cutoff)

        # More lenient settings shouldn't return results
        assert not RREFinderResults.from_json(self.json_dict, self.record, 25, self.bitscore_cutoff)
        assert not RREFinderResults.from_json(self.json_dict, self.record, self.min_length, 15.0)

        # Stricter settings should filter
        stricter_score_cutoff_hits = {1: {'locus_tag_a': [self.hit_a1],}, }
        res_object2 = RREFinderResults.from_json(self.json_dict, self.record, self.min_length, 35.0)
        assert res_object2.hits_per_protocluster == stricter_score_cutoff_hits

        stricter_length_cutoff_hits = {2: {'locus_tag_b': [self.hit_b1],}, }
        res_object3 = RREFinderResults.from_json(self.json_dict, self.record, 80, self.bitscore_cutoff)
        assert res_object3.hits_per_protocluster == stricter_length_cutoff_hits

class TestRREFinder(unittest.TestCase):

    def setUp(self):
        self.ripps = ['bacteriocin','cyanobactin','lanthipeptide',
             'lassopeptide','linaridin','thiopeptide','sactipeptide',
              'proteusin','glycocin','bottromycin','microcin']

        self.hit_a1 = {'locus_tag': 'locus_tag_a', 'location': '1000-2000',
                        'label': 'a_generic_label', 'domain': 'RRE_type_A', 'evalue': 0.1,
                        'protein_start': 1, 'protein_end': 50, 'score': 38.0,

                        }
        self.hit_a2 = {'locus_tag': 'locus_tag_a',
                        'label': 'a_generic_label', 'domain': 'RRE_type_B', 'evalue': 0.1,
                        'protein_start': 1, 'protein_end': 40, 'score': 38.0
                        }
        self.hit_b1 = {'locus_tag': 'locus_tag_b',
                        'label': 'a_generic_label', 'domain': 'RRE_type_B', 'evalue': 0.2,
                        'protein_start': 1, 'protein_end': 90, 'score': 25.0
                        }
        self.hit_b2 = {'locus_tag': 'locus_tag_b',
                        'label': 'a_generic_label', 'domain': 'RRE_type_C', 'evalue': 0.2,
                        'protein_start': 15, 'protein_end': 50, 'score': 42.0
                        }

        self.hit_c =  {'locus_tag': 'locus_tag_c',
                        'label': 'a_generic_label', 'domain': 'RRE_type_D', 'evalue': 0.2,
                        'protein_start': 1, 'protein_end': 120, 'score': 6.0
                        }
        
        self.unfiltered_hits = {1: {'locus_tag_a': [self.hit_a1, self.hit_a2],
                                      'locus_tag_c': [self.hit_c],
                                     },
                                2: {'locus_tag_b': [self.hit_b1, self.hit_b2],
                                     }
                               }

        self.filtered_hits = {1: {'locus_tag_a': [self.hit_a1],
                                    },
                               2: {'locus_tag_b': [self.hit_b1],
                                    },
                              }

        self.bitscore_cutoff = 25.0
        self.max_evalue = 1
        self.min_length = 50
        self.record_id = 'test_record'
        self.mock_record = self.make_mock_record()
        self.database = 'rre_database'
        self.tool = 'rrefinder'

        self.hmm_res_all = {1: HmmerResults(self.record_id, self.max_evalue, self.bitscore_cutoff, self.database, self.tool, [self.hit_a1, self.hit_a2, self.hit_c]),
                            2: HmmerResults(self.record_id, self.max_evalue, self.bitscore_cutoff, self.database, self.tool, [self.hit_b1, self.hit_b2]) }
        self.res_object = RREFinderResults(self.record_id, self.bitscore_cutoff, self.min_length, self.filtered_hits)

    def make_mock_record(self):
        # Simple mock record containing a mock region, which contains three mock protoclusters
        # protocluster.cds_children is the same as the protocluster_number, so that results
        # can be retrieved in testing hmmscan
        region = Mock(name='my_region')
        record = Mock(name='my_record')

        protoclusters = []
        for protocluster_number in range(1,4):
            protocluster = Mock(name = 'my_protocluster_%i' %protocluster_number)
            protocluster.cds_children = protocluster_number
            protocluster.get_protocluster_number.return_value = protocluster_number
            protoclusters.append(protocluster)
            if protocluster_number == 3:
                protocluster.product = 'Not_a_RiPP'
            else:
                protocluster.product = self.ripps[0]

        region.get_unique_protoclusters.return_value = protoclusters
        record.get_regions.return_value = [region, ]
        return record

    def return_hmm_results(self, record, protocluster_cds_children, max_evalue, bitscore_cutoff, hmm_database, tool):
        # The protocluster.cds_children attribute has been mocked to return an integer so that the results can be retrieved in this test
        return self.hmm_res_all[protocluster_cds_children]

    def test_is_ripp(self):
        for ripp in self.ripps:
            assert is_ripp(ripp)
            assert not is_ripp(ripp[1:])

    def test_will_handle(self):
        assert will_handle(self.ripps)
        assert not will_handle([ripp[1:] for ripp in self.ripps])

    def test_check_hmm_hit(self):
        assert check_hmm_hit(self.hit_a1, self.min_length, self.bitscore_cutoff)
        assert not check_hmm_hit(self.hit_a2, self.min_length, self.bitscore_cutoff)
        assert check_hmm_hit(self.hit_b1, self.min_length, self.bitscore_cutoff)
        assert not check_hmm_hit(self.hit_b2, self.min_length, self.bitscore_cutoff)
        assert not check_hmm_hit(self.hit_c, self.min_length, self.bitscore_cutoff)

    def test_filter_hits(self):
        filtered_hits = filter_hits(self.unfiltered_hits, self.min_length, self.bitscore_cutoff)
        assert self.filtered_hits == filtered_hits
        assert not filter_hits(self.unfiltered_hits, 1000, 10)
        assert filter_hits(self.unfiltered_hits, 0, 0) == self.unfiltered_hits

    @patch('antismash.modules.rrefinder.rrefinder.run_hmmer_copy')
    def test_run_hmmscan_rrefinder(self, mocked_function):
        mocked_function.side_effect = self.return_hmm_results
        hmm_results = run_hmmscan_rrefinder(self.mock_record, self.bitscore_cutoff)
        assert hmm_results == self.unfiltered_hits

    def test_run_rrefinder(self):
        pass

