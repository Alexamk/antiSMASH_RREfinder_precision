# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Replace this text with a description of the module.
    It can also include references for the method implemented.
"""

# start with standard library imports
import logging
from typing import Any, Dict, List, Optional

# then any imports from external modules, e.g. biopython, if relevant

# then any imports from antismash
from antismash.common.secmet import Record
from antismash.config import ConfigType, get_config
from antismash.config.args import ModuleArgs

from antismash.common.hmmer import ensure_database_pressed
from antismash.common import path

from .html_output import generate_html, will_handle
from .rrefinder import run_rrefinder, RREFinderResults

# then any local file imports, e.g. from .somefile import..., if relevant

NAME = "rrefinder"
SHORT_DESCRIPTION = "Module for the pHMM-based detection of RiPP Recognition Elements (RREs)"


# define a results class, this is important as adding information to the record
# during analysis will cause issues
# for detailed examples, see any of the other analysis modules' implementations


def get_arguments() -> ModuleArgs:
    """ Builds any commandline argument constructs that may be required

        Returns:
            an empty or populated ModuleArgs instance
    """
    # construct the argument group, with section and prefix
    # the prefix will be enforced for all command line options for the module
    args = ModuleArgs('Additional analysis', 'rre')

    # an example toggle to turn on your analysis, if not set to always be enabled
    # can also be used to turn on/off extra features of your analysis
    args.add_analysis_toggle('run',     # the option as it appears on the command line
                             dest='run',  # the storage location in the antismash Config object
                             default=False,             # disabled by default
                             action='store_true',       # enabled if --template-analysis is given on the commandline
                             # and finally, text to show when the user runs with --help
                             help="Run RREFinder precision mode on all RiPP gene clusters.")

    # an example option setting a specific value
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
    # more complicated options are possible, for further information see antismash.config.args,
    # look at how other modules construct arguments, or ask for help
    return args


def check_options(options: ConfigType) -> List[str]:
    """ Checks that the provided options are compatible with each other

        Arguments:
            options: the current antismash config object

        Returns:
            a list of strings, each string being an issue with the given options
    """
    issues = []
    # test the options used in get_arguments here, if any
    # for example, enforcing that values are within a certain range
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
    # behaves similarly to check_options(), though checking for built databases,
    # that external programs are available, and so on
    # see antismash.detection.hmm_detection for an example of these
    
    hmm_database = path.get_full_path(__file__, 'data', 'RREFam.hmm')
    hmm_present = ensure_database_pressed(hmm_database, return_not_raise=True)
    # Currently the only error
    
    # if there are no external prerequisites, this can just return an empty list
    return hmm_present


def is_enabled(options: ConfigType) -> bool:
    """ Returns True if the module is enabled with the options provided
    """
    # the logic here depends on which command options you've created
    # using the example above, this is as simple as returning the toggle
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
    results = run_rrefinder(record, options.rre_cutoff, options.rre_min_length)
    
    
    # and return it
    return results
