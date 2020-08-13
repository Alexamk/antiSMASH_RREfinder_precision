# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A feature to represent a RiPP Recognition Element (RRE) """

import logging

from collections import OrderedDict
from typing import Any, Dict, List, Optional, Type, TypeVar

from Bio.SeqFeature import SeqFeature

from .domain import Domain, generate_protein_location_from_qualifiers
from .feature import Feature, Location


T = TypeVar("T", bound="RRE")

class RRE(Domain):
    """ A feature representing an RRE within a CDS.
    """

    __slots__ = ['description']
    FEATURE_TYPE = 'RRE'

    def __init__(self, location: Location, description: str, protein_location: Location,
                 tool: str, locus_tag: str, domain: Optional[str] = None,
                 ) -> None:

        """ Arguments:
                location: the DNA location of the feature
                description: a string with a description of the RRE (e.g. 'RRE-containing protein in a pyrroloquinoline cluster')
                protein_location: the location within the parent CDS translation
                # identifier: the RREfam identifier (e.g. RREFam005)
                tool: the name of the tool used to find/create the feature (rrefinder)
                locus_tag: the name of the parent CDS feature
                domain: the name for the domain (e.g. Lanthipeptide_RRE)
        """
        super().__init__(location, self.FEATURE_TYPE, protein_location, locus_tag, domain,
                       tool=tool, created_by_antismash=True)
        if not isinstance(description, str):
            raise TypeError("RRE description must be a string, not %s" % type(description))
        if not description:
            raise ValueError("RRE description cannot be empty")
        self.description = description

    def to_biopython(self, qualifiers: Dict[str, List[str]] = None) -> List[SeqFeature]:
        mine = OrderedDict()  # type: Dict[str, List[str]]
        mine["description"] = [self.description]
        if qualifiers:
            mine.update(qualifiers)
        return super().to_biopython(mine)

    @classmethod
    def from_biopython(cls: Type[T], bio_feature: SeqFeature, feature: T = None,
                       leftovers: Dict[str, List[str]] = None, record: Any = None) -> T:

        if leftovers is None:
            leftovers = Feature.make_qualifiers_copy(bio_feature)

        tool = leftovers.pop("aSTool")[0]
        protein_location = generate_protein_location_from_qualifiers(leftovers, record)
        # Remove the protein_start and protein_end from the leftovers
        _ = leftovers.pop('protein_start')
        _ = leftovers.pop('protein_end')

        description = leftovers.pop('description')[0]
        locus_tag = leftovers.pop("locus_tag", ["(unknown)"])[0]
        feature = cls(bio_feature.location, description, protein_location, tool, locus_tag)

        # grab parent optional qualifiers
        super().from_biopython(bio_feature, feature=feature, leftovers=leftovers, record=record)

        return feature
