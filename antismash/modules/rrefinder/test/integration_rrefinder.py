# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import os
import unittest

from tempfile import TemporaryDirectory

import antismash
from antismash.common.secmet import Record
from antismash.common.secmet.features import FeatureLocation
from antismash.common.test import helpers
from antismash.config import get_config, update_config, destroy_config, build_config
from antismash.main import read_data
from antismash.modules import rrefinder

class RREFinderIntegration(unittest.TestCase):
    
    def setUp(self):
        self.options = build_config(["--minimal", "--rre-run"],
                       isolated=True, modules=antismash.get_all_modules())

    def tearDown(self):
        destroy_config()
        
    def test_options(self):
        assert rrefinder.check_options(self.options) == []
        assert rrefinder.check_prereqs(self.options) == []
        assert rrefinder.is_enabled(self.options)
        
        options = build_config(["--minimal", "--rre-run", "--rre-cutoff", "-10"],
                  isolated=True, modules=antismash.get_all_modules())
        issues = rrefinder.check_options(options)
        assert len(issues) == 1
        assert issues[0] == 'Supplied RREFinder cutoff is negative: %s' %options.rre_cutoff
        
    def test_nisin(self):
        genbank = helpers.get_path_to_nisin_with_detection()
        results = helpers.run_and_regenerate_results_for_module(genbank, rrefinder, self.options)
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
        assert feature.translation == 'FTRNNANYKFGDRVFQVYTINSSELEEVNIKYTNVYQIISEFCENDYQKYEDICETVTLCYGDEYRELSEQYLGSLIVNHYLISNLQK'