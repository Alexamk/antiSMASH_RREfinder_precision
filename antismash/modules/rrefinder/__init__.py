# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Implementation of RREFinder in antiSMASH.
    RiPP Recognition Elements (RREs) are protein domains of ~100
    amino acids in RiPP modifying enzymes that are vital for RiPP
    precursor peptide recognition. They are found in roughly half
    of all bacterial RiPP classes.
    RREFinder is a tool for the recognition of these domains.
    RREFinder has two modes for RRE detection: precision
    mode, which detects RREs of known RiPP classes with high
    confidence using a library of pHMMs, and exploratory mode,
    which uses HHPred to detect RREs. Exploratory mode detects
    more false positives, but is capable of detecting RREs for
    which no specific pHMM model has been designed.

    Only precision mode has been implemented here.

    For more details, see:
    https://www.biorxiv.org/content/10.1101/2020.03.14.992123v1
    and
    https://github.com/Alexamk/RREFinder
"""

import logging
import os

from typing import Any, Dict, List, Optional

from antismash.common.secmet import Record
from antismash.config import ConfigType, get_config
from antismash.config.args import ModuleArgs

from antismash.common.hmmer import ensure_database_pressed
from antismash.common import path

from .html_output import generate_html, will_handle
from .rrefinder import run_rrefinder, RREFinderResults

NAME = "rrefinder"
SHORT_DESCRIPTION = "Module for the pHMM-based detection of RiPP Recognition Elements (RREs)"

def get_arguments() -> ModuleArgs:
    """ Builds any commandline argument constructs that may be required

        Returns:
            an empty or populated ModuleArgs instance
    """

    args = ModuleArgs('Additional analysis', 'rre')
    args.add_analysis_toggle('run',     # the option as it appears on the command line
                             dest='run',  # the storage location in the antismash Config object
                             default=False,             # disabled by default
                             action='store_true',       # enabled if --template-analysis is given on the commandline
                             help="Run RREFinder precision mode on all RiPP gene clusters.")
    args.add_option('cutoff',     # the option as it appears on the command line
                    dest='cutoff',  # as it appears in the antismash Config object
                    type=float,              # the type of the option (int, str, float, ...)
                    default=25.0,            # the default value of the option
                    help="Bitscore cutoff for RRE pHMM detection.")
    args.add_option('min_length',
                    dest='min_length',
                    type=int,
                    default=50,
                    help='Minimum amino acid length of RRE domains.')
    return args


def check_options(options: ConfigType) -> List[str]:
    """ Checks that the provided options are compatible with each other

        Arguments:
            options: the current antismash config object

        Returns:
            a list of strings, each string being an issue with the given options
    """
    issues = []
    if options.rre_cutoff <= 0:
        issues.append("Supplied RREFinder cutoff is negative: %s" % options.rre_cutoff)
    return issues


def check_prereqs(options: ConfigType) -> List[str]:
    """ Check that all prerequisites are present

        Arguments:
            options: the current antismash config object

        Returns:
            a list of strings, each string being an issue with prerequisites
    """
    # Currently the only requirement is that the hmm_database is present and pressed
    failure_messages = []
    failure_messages.extend(prepare_data(logging_only=True))
    return failure_messages


def prepare_data(logging_only: bool = False) -> List[str]:
    """ Ensures packaged data is fully prepared

        Arguments:
            logging_only: whether to return error messages instead of raising exceptions

        Returns:
            a list of error messages (only if logging_only is True)
    """
    database = path.get_full_path(__file__, 'data', 'RREFam.hmm')
    return ensure_database_pressed(database, return_not_raise=logging_only)


def is_enabled(options: ConfigType) -> bool:
    """ Returns True if the module is enabled with the options provided
    """
    return options.rre_run


def regenerate_previous_results(previous: Dict[str, Any], record: Record,
                                options: ConfigType) -> Optional[RREFinderResults]:
    """ Regenerate the previous results from JSON format.

        Arguments:
            previous: the previous results as from JSON
            record: the Record these previous results were originally created from
            options: the current antismash config object

        Returns:
            an instance of the module's ModuleResults implementation,
            or None if the current options require the analysis to be rerun or cannot be regenerated
    """
    # if there isn't anything to work with, just return None
    if not previous:
        return None
    return RREFinderResults.from_json(previous, record, options.rre_min_length, options.rre_cutoff)


def run_on_record(record: Record, results: RREFinderResults, options: ConfigType) -> RREFinderResults:
    """ Run the analysis, unless the previous results apply to the given record

        Arguments:
            record: the Record being analysed
            results: an existing instance of the module's ModuleResults implementation (or None)
            options: the current antismash config object

        Returns:
            an instance of the module's ModuleResults implementation
    """
    # after a safety check that the results are the correct ones for the record, return them
    if isinstance(results, RREFinderResults) and results.record_id == record.id:
        return results
    # otherwise run the actual analysis and generate a results instance with your analysis results
    database = path.get_full_path(__file__, 'data', 'RREFam.hmm')
    results = run_rrefinder(record, options.rre_cutoff, options.rre_min_length, database)
    return results
