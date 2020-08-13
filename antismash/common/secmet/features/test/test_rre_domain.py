# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from antismash.common.secmet.features import FeatureLocation, RRE

class TestRRE(unittest.TestCase):

    def setUp(self):
        self.protein_location = FeatureLocation(1,5)
        self.location = FeatureLocation(6,10)
        self.tool = 'rrefinder_test'
        self.domain = 'RRE_type_a'
        self.description = 'This is a test RRE'
        self.locus_tag = 'locus_tag_a'
        self.rre = RRE(self.location, self.description, self.protein_location,
                       self.tool, self.locus_tag, self.domain)

    def test_init(self):
        assert self.rre.locus_tag == self.locus_tag
        assert self.rre.description == self.description
        assert self.rre.tool == self.tool
        assert self.rre.domain == self.domain
        assert self.rre.location == self.location
        assert self.rre.protein_location == self.protein_location

        # Test wrong description entries
        with self.assertRaises(TypeError):
            rre = RRE(self.location, 5, self.protein_location,
                      self.tool, self.locus_tag, self.domain)
        with self.assertRaisesRegex(ValueError, "RRE description cannot be empty"):
            rre = RRE(self.location, '', self.protein_location,
                      self.tool, self.locus_tag, self.domain)

    def test_conversion(self):
        bio = self.rre.to_biopython()

        assert len(bio) == 1
        assert bio[0].qualifiers['description'][0] == self.description

        rre = RRE.from_biopython(bio[0])
        assert rre.description == self.rre.description
        assert rre.tool == self.rre.tool
        assert rre.protein_location == self.rre.protein_location
        assert rre.locus_tag == self.rre.locus_tag
