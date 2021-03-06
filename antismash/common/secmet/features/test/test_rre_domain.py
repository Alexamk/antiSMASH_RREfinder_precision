# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from antismash.common.secmet.features import FeatureLocation, RRE

class TestRRE(unittest.TestCase):

    def setUp(self):
        self.protein_location = FeatureLocation(1, 5)
        self.location = FeatureLocation(6, 10)
        self.tool = 'rrefinder_test'
        self.domain = 'RRE_type_a'
        self.description = 'This is a test RRE'
        self.locus_tag = 'locus_tag_a'
        self.identifier = 'RREFam001'
        self.rre = RRE(self.location, self.description, self.protein_location,
                       self.identifier, self.tool, self.locus_tag, self.domain)

    def test_init(self):
        assert self.rre.locus_tag == self.locus_tag
        assert self.rre.description == self.description
        assert self.rre.tool == self.tool
        assert self.rre.domain == self.domain
        assert self.rre.location == self.location
        assert self.rre.protein_location == self.protein_location
        assert self.rre.identifier == self.identifier

        # Test wrong description entries
        with self.assertRaises(TypeError):
            _ = RRE(self.location, 5, self.protein_location, self.identifier,
                      self.tool, self.locus_tag, self.domain)
        with self.assertRaisesRegex(ValueError, "RRE description cannot be empty"):
            _ = RRE(self.location, '', self.protein_location, self.identifier,
                      self.tool, self.locus_tag, self.domain)
        with self.assertRaisesRegex(ValueError, "RREFam identifier cannot be empty"):
            _ = RRE(self.location, self.description, self.protein_location, '',
                      self.tool, self.locus_tag, self.domain)
        with self.assertRaises(ValueError):
            _ = RRE(self.location, self.description, self.protein_location, 'not_a_valid_identifier',
                      self.tool, self.locus_tag, self.domain)

    def test_conversion(self):
        bio = self.rre.to_biopython()

        assert len(bio) == 1
        assert bio[0].qualifiers['description'][0] == self.description
        assert bio[0].qualifiers['identifier'][0] == self.identifier

        rre = RRE.from_biopython(bio[0])
        assert rre.description == self.rre.description
        assert rre.tool == self.rre.tool
        assert rre.protein_location == self.rre.protein_location
        assert rre.locus_tag == self.rre.locus_tag
        assert rre.identifier == self.rre.identifier

        # Test with extra qualifiers
        bio = self.rre.to_biopython(qualifiers={'some_qualifier': ['some_value']})
        assert bio[0].qualifiers.get('some_qualifier')
        assert bio[0].qualifiers['some_qualifier'][0] == 'some_value'
