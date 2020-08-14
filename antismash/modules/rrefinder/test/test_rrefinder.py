import unittest

from collections import defaultdict
from unittest.mock import Mock, patch

from antismash.common.hmmer import HmmerResults
from antismash.common.secmet.features import RRE, FeatureLocation

from antismash.modules.rrefinder.html_output import will_handle
from antismash.modules.rrefinder.rrefinder import (
    is_ripp,
    check_hmm_hit,
    gather_rre_candidates,
    extract_rre_hits,
    filter_hits,
    RREFinderResults,
    run_rrefinder,
    )

class TestRREResults(unittest.TestCase):

    def setUp(self):
        self.hit_a1 = {'locus_tag': 'locus_tag_a', 'location': '[1000:2000]',
                        'label': 'a1_generic_label', 'domain': 'RRE_type_A', 'evalue': 0.1,
                        'protein_start': 0, 'protein_end': 50, 'score': 38.0,
                        'description': 'description_type_A', 'translation': 'FAKESEQ',
                      }
        self.hit_b1 = {'locus_tag': 'locus_tag_b', 'location': '[500:1500]',
                        'label': 'b1_generic_label', 'domain': 'RRE_type_B', 'evalue': 0.2,
                        'protein_start': 0, 'protein_end': 90, 'score': 25.0,
                        'description': 'description_type_B', 'translation': 'FAKESEQ',
                        }

        self.hit_info = dict((hit['locus_tag'], [hit]) for hit in [self.hit_a1, self.hit_b1])
        self.hits_per_protocluster = {1: ['locus_tag_a'],
                                      2: ['locus_tag_b']}
        self.tool = 'rrefinder'
        self.database = 'RREFam.hmm'
        self.detection = 'hmmscan'

        self.min_length = 50
        self.bitscore_cutoff = 25.0

        self.record = Mock()
        self.record.id = self.record_id = 'test_record'

        self.record2 = Mock()
        self.record2.record_id = 'another_record'

        self.json_dict = {"schema_version": RREFinderResults.schema_version,
                          "bitscore_cutoff": self.bitscore_cutoff,
                          "hits_per_protocluster": self.hits_per_protocluster,
                          "hit_info": self.hit_info,
                          "min_length": self.min_length,
                          "record_id": self.record_id, }

        self.res_object = RREFinderResults(self.record_id, self.bitscore_cutoff, self.min_length, 
                                           self.hits_per_protocluster, self.hit_info)

    def test_init(self):
        assert self.res_object.record_id == self.record_id
        assert self.res_object.bitscore_cutoff == self.bitscore_cutoff
        assert self.res_object.min_length == self.min_length
        assert self.res_object.hits_per_protocluster == self.hits_per_protocluster
        assert self.res_object.hit_info == self.hit_info
        assert self.res_object.tool == self.tool
        assert self.res_object.detection == self.detection
        assert self.res_object.database == self.database
        assert self.res_object.features

    def test_convert_hits_to_features(self):
        hit_info = self.hit_info.copy()
        _ = hit_info.pop('locus_tag_b')
        res_object = RREFinderResults(self.record_id, self.bitscore_cutoff, self.min_length, self.hits_per_protocluster, hit_info)
        assert len(res_object.features) == 1
        feature = res_object.features[0]
        assert isinstance(feature, RRE)
        assert feature.location == FeatureLocation(1000, 2000)
        assert feature.protein_location == FeatureLocation(0, 50)
        assert feature.description == 'description_type_A'
        assert feature.domain == 'RRE_type_A'
        assert feature.locus_tag == 'locus_tag_a'
        assert feature.domain_id == '{}_locus_tag_a_0001'.format(self.tool)
        assert feature.database == self.database
        assert feature.detection == self.detection
        assert feature.score == self.hit_a1['score']
        assert feature.evalue == self.hit_a1['evalue']
        assert feature.label == self.hit_a1['label']
        assert feature.translation == self.hit_a1['translation']

    def test_to_json(self):
        assert self.json_dict == self.res_object.to_json()

    def test_from_json(self):
        res_object = RREFinderResults.from_json(self.json_dict, self.record, self.min_length, self.bitscore_cutoff)
        assert self.res_object
        assert self.res_object.hits_per_protocluster == self.hits_per_protocluster
        assert self.res_object.hit_info == self.hit_info
        assert self.res_object.bitscore_cutoff == self.bitscore_cutoff
        assert self.res_object.min_length == self.min_length
        assert self.res_object.record_id == self.record_id

        json_dict2 = dict(self.json_dict)
        json_dict2['record_id'] = 'another_record'
        assert not RREFinderResults.from_json(json_dict2, self.record, self.min_length, self.bitscore_cutoff)

        json_dict3 = dict(self.json_dict)
        json_dict3['schema_version'] = 'not_a_valid_schema'
        assert not RREFinderResults.from_json(json_dict3, self.record, self.min_length, self.bitscore_cutoff)

        # Invalid arguments should return None
        json_dict4 = dict(self.json_dict)
        _ = json_dict4.pop('bitscore_cutoff')
        with self.assertRaises(ValueError):
            _ = RREFinderResults.from_json(json_dict4, self.record, self.min_length, self.bitscore_cutoff)

        json_dict5 = dict(self.json_dict)
        _ = json_dict5.pop('min_length')
        with self.assertRaises(ValueError):
            _ = RREFinderResults.from_json(json_dict5, self.record, self.min_length, self.bitscore_cutoff)

        # More lenient settings shouldn't return results
        assert not RREFinderResults.from_json(self.json_dict, self.record, 25, self.bitscore_cutoff)
        assert not RREFinderResults.from_json(self.json_dict, self.record, self.min_length, 15.0)

        # Stricter settings should filter
        res_object2 = RREFinderResults.from_json(self.json_dict, self.record, self.min_length, 35.0)
        assert res_object2.hit_info['locus_tag_a'] == [self.hit_a1]
        assert len(res_object2.hit_info) == 1
        assert res_object2.hits_per_protocluster[1] == ['locus_tag_a']
        assert len(res_object2.hits_per_protocluster) == 1

        res_object3 = RREFinderResults.from_json(self.json_dict, self.record, 80, self.bitscore_cutoff)
        assert res_object3.hits_per_protocluster[2] == ['locus_tag_b']
        assert len(res_object3.hits_per_protocluster) == 1
        assert res_object3.hit_info['locus_tag_b'] == [self.hit_b1]
        assert len(res_object3.hit_info) == 1
        
    def test_add_record(self):
        with self.assertRaisesRegex(ValueError, "Record to store in and record analysed don't match"):
            self.res_object.add_to_record(self.record2)
        self.res_object.add_to_record(self.record)
        assert self.record.add_feature.call_count == 2
        for arg in self.record.add_feature.call_args_list:
            # First argument of positional arguments
            assert isinstance(arg[0][0], RRE)

class TestRREFinder(unittest.TestCase):

    def setUp(self):
        self.ripps = ['bacteriocin','cyanobactin','lanthipeptide',
             'lassopeptide','linaridin','thiopeptide','sactipeptide',
              'proteusin','glycocin','bottromycin','microcin']

        self.hit_a1 = {'locus_tag': 'locus_tag_a', 'location': '[1000:2000]',
                        'label': 'a1_generic_label', 'domain': 'RRE_type_A', 'evalue': 0.1,
                        'protein_start': 0, 'protein_end': 50, 'score': 38.0,
                        'description': 'description_type_A', 'translation': 'FAKESEQ',
                      }
        self.hit_a2 = {'locus_tag': 'locus_tag_a', 'location': '[1100:2100]',
                        'label': 'a1_generic_label', 'domain': 'RRE_type_A', 'evalue': 0.1,
                        'protein_start': 0, 'protein_end': 40, 'score': 38.0,
                        'description': 'description_type_A', 'translation': 'FAKESEQ',
                      }
        self.hit_b1 = {'locus_tag': 'locus_tag_b', 'location': '[500:1500]',
                        'label': 'b1_generic_label', 'domain': 'RRE_type_B', 'evalue': 0.2,
                        'protein_start': 0, 'protein_end': 90, 'score': 25.0,
                        'description': 'description_type_B', 'translation': 'FAKESEQ',
                        }
        self.hit_b2 = {'locus_tag': 'locus_tag_b', 'location': '[600:1600]',
                        'label': 'b1_generic_label', 'domain': 'RRE_type_C', 'evalue': 0.2,
                        'protein_start': 14, 'protein_end': 50, 'score': 42.0,
                        'description': 'description_type_C', 'translation': 'FAKESEQ',
                      }
        self.hit_c =  {'locus_tag': 'locus_tag_c', 'location': '[200:400]',
                        'label': 'a_generic_label', 'domain': 'RRE_type_D', 'evalue': 0.2,
                        'protein_start': 0, 'protein_end': 120, 'score': 6.0,
                        'description': 'description_type_D', 'translation': 'FAKESEQ',
                      }

        self.candidates_per_protocluster = defaultdict(list)
        self.candidates_per_protocluster[1] = ['locus_tag_a', 'locus_tag_c']
        self.candidates_per_protocluster[2] = ['locus_tag_b']
                                           
        self.filtered_hits_per_protocluster = defaultdict(list)
        self.filtered_hits_per_protocluster[1] = ['locus_tag_a']
        self.filtered_hits_per_protocluster[2] = ['locus_tag_b']
        
        self.all_hits = [self.hit_a1, self.hit_a2, self.hit_b1, self.hit_b2, self.hit_c]
        self.fake_cds_info = dict((locus_tag, locus_tag) for locus_tag in ['locus_tag_a', 'locus_tag_b', 'locus_tag_c'])
        self.hit_info = defaultdict(list)
        for hit in self.all_hits:
            self.hit_info[hit['locus_tag']].append(hit)
        self.filtered_hit_info = dict((hit['locus_tag'], [hit]) for hit in [self.hit_a1, self.hit_b1])
        
        self.bitscore_cutoff = 25.0
        self.max_evalue = 1
        self.min_length = 50
        self.record_id = 'test_record'
        self.mock_record = self.make_mock_record()
        self.database = 'rre_database'
        self.tool = 'rrefinder'
        
        self.hmm_res = HmmerResults(self.record_id, self.max_evalue, self.bitscore_cutoff, self.database, self.tool, self.all_hits)
    
    def make_mock_record(self):
        # Simple mock record containing a mock region, which contains three mock protoclusters
        # protocluster.cds_children is the same as the protocluster_number, so that results
        # can be retrieved in testing hmmscan
        region = Mock(name='my_region')
        record = Mock(name='my_record')
        record.id = self.record_id

        protoclusters = []
        for protocluster_number in range(1,4):
            protocluster = Mock(name = 'my_protocluster_%i' %protocluster_number)
            protocluster.get_protocluster_number.return_value = protocluster_number
            if protocluster_number == 3:
                protocluster.product = 'Not_a_RiPP'
            else:
                protocluster.product = self.ripps[0]
                # mocked protocluster just return the locus tag as cds_children
                protocluster.cds_children = []
                for locus_tag in self.candidates_per_protocluster[protocluster_number]:
                    cds = Mock()
                    cds.get_name.return_value = locus_tag
                    protocluster.cds_children.append(cds)
            protoclusters.append(protocluster)

        region.get_unique_protoclusters.return_value = protoclusters
        record.get_regions.return_value = [region, ]
        return record

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
        filtered_hit_info, filtered_hits_per_protocluster = filter_hits(self.hit_info, self.candidates_per_protocluster, self.min_length, self.bitscore_cutoff)
        assert self.filtered_hit_info == filtered_hit_info
        assert self.filtered_hits_per_protocluster == filtered_hits_per_protocluster
        empty_res = filter_hits(self.hit_info, self.candidates_per_protocluster, 1000, 10)
        assert not(empty_res[0]) and not(empty_res[1]) 
        unfiltered_res = filter_hits(self.hit_info, self.candidates_per_protocluster, 0, 0)
        assert unfiltered_res[0] == self.hit_info
        assert unfiltered_res[1] == self.candidates_per_protocluster

    def test_extract_rre_hits(self):
        assert self.hit_info == extract_rre_hits(self.hmm_res)
        
    def test_gather_rre_candidates(self):
        candidates_per_protocluster, cds_info = gather_rre_candidates(self.mock_record)
        assert self.candidates_per_protocluster == candidates_per_protocluster
        for key in self.fake_cds_info:
            assert key in cds_info
            assert key == cds_info[key].get_name()
        
    @patch('antismash.modules.rrefinder.rrefinder.run_hmmer_copy')
    def test_run_rrefinder(self, mocked_function):
        mocked_function.return_value = self.hmm_res
        res_object = run_rrefinder(self.mock_record, self.bitscore_cutoff, self.min_length)
        assert res_object.record_id == self.record_id
        assert res_object.hits_per_protocluster == self.filtered_hits_per_protocluster
        assert res_object.hit_info == self.filtered_hit_info
        assert res_object.bitscore_cutoff == self.bitscore_cutoff
        assert res_object.min_length == self.min_length
        assert mocked_function.call_count == 1
        args = mocked_function.call_args
        cds_args = args[0][1]
        for cds in cds_args:
            assert cds.get_name.call_count == 1
