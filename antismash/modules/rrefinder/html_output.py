# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Manages HTML construction for the RREFinder module
"""

from typing import List

from antismash.common import path
from antismash.common.html_renderer import HTMLSections, FileTemplate
from antismash.common.layers import RegionLayer, RecordLayer, OptionsLayer

from .rrefinder import is_ripp, RREFinderResults


def will_handle(products: List[str]) -> bool:
    """ Returns True if the products provided are relevant to the module """
    return any([is_ripp(product) for product in products])


def generate_html(region_layer: RegionLayer, results: RREFinderResults,
                  _record_layer: RecordLayer, _options_layer: OptionsLayer
                  ) -> HTMLSections:
    """ Generates HTML output for the module """
    html = HTMLSections("rrefinder")

    side_tooltip = ("A beautiful description.")
    template = FileTemplate(path.get_full_path(__file__, "templates", "sidepanel.html"))
    
#    region_hits = results.get_subset_hits_for_region(region_layer.region_feature)
    
    html.add_sidepanel_section("RREFinder", template.render(results=results, region=region_layer.region_feature, tooltip=side_tooltip),"RREfinder")

    return html
