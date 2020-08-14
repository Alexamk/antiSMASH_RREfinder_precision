# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from tempfile import TemporaryDirectory

import antismash
from antismash.common.secmet.features import FeatureLocation
from antismash.common.test import helpers
from antismash.config import get_config, update_config, destroy_config, build_config
from antismash.modules import rrefinder

class RREFinderIntegration(unittest.TestCase):
    
    def setUp(self):
        self.options = build_config(["--minimal", "--rre-run"],
                       isolated=True, modules=antismash.get_all_modules())

    def tearDown(self):
        destroy_config()
        
    def test_options(self):
        assert not rrefinder.check_options(self.options)
        assert rrefinder.check_prereqs(self.options)
        assert rrefinder.is_enabled(self.options)
        
        options = build_config(["--minimal", "--rre-run", "--rre-cutoff", "-10"],
                  isolated=True, modules=antismash.get_all_modules())
        assert rrefinder.check_options(options)
        
    def test_nisin(self):
        record = Record.from_genbank(helpers.get_path_to_nisin_with_detection(), taxon="bacteria")[0]
        regions = record.get_regions()
        assert regions
        for region in regions:
            assert region.cds_children
        assert record.get_cds_features_within_regions()
        before_count = record.get_feature_count()
        
        prior_results = None
        results = rrefinder.run_on_record(record, prior_results, self.options)
        assert isinstance(results, rrefinder.RREFinderResults)
        assert results.bitscore_cutoff == self.options.rre_cutoff
        assert results.min_length == self.options.rre_min_length
        assert results.record_id == 'HM219853.1'
        assert len(results.features) == 1
        feature = results.features[0]
        assert feature.locus_tag == 'nisB'
        assert feature.score == 34.9
        assert feature.protein_location == FeatureLocation(141, 228)
        assert feature.domain == 'Lanthipeptide_RRE'
        assert record.get_feature_count() == before_count
        results.add_to_record(record)
        assert record.get_feature_count() == before_count + 1
        
    def test_regenerate_nisin(self):
        with TemporaryDirectory() as output_dir:
            args = ["--minimal", "--rre-run", "--output-dir", output_dir, helpers.get_path_to_nisin_genbank()]
            options = build_config(args, isolated=True, modules=antismash.get_all_modules())
            antismash.run_antismash(helpers.get_path_to_nisin_genbank(), options)

            # regen the results
            update_config({"reuse_results": os.path.join(output_dir, "nisin.json")})
            prior_results = read_data(None, options)
            record = prior_results.records[0]
            results = prior_results.results[0]
            rre_results = rrefinder.regenerate_previous_results(results.get("antismash.modules.rrefinder"), record, options)
            assert isinstance(rre_results, rrefinder.RREFinderResults)
            assert len(rre_results.features) == 1
        
