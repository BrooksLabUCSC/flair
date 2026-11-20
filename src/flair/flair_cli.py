#! /usr/bin/env python3
import time
from datetime import timedelta
import argparse
import logging
from flair import VERSION
from flair.cli_help import FlairHelpFormatter
from flair.cli_retired import RetiredCheckParser
from flair.run_logging import configure_logging, log_run_start
from flair.run_signals import install_stop_handlers
from flair.output_compress import set_output_compression
from flair.pycbio.sys import cli, loggingOps
from flair import (flair_combine, flair_diffExp, flair_diffSplice, flair_fusion,
                   flair_isoalleles, flair_quantify, flair_spliceevents, flair_transcriptome)

# Modules in pipeline order.  Each one adds and then owns its own subparser and
# records its entry point with set_defaults(entry=), so no parser is filled in by
# one module and parsed by another.
#
# These are imported eagerly.  Importing all eight costs about 0.4 seconds, most of
# it numpy in flair_diffExp, which is not worth the pre-parsing needed to load
# them one at a time.
SUBCOMMAND_MODULES = (
    flair_transcriptome,
    flair_quantify,
    flair_combine,
    flair_isoalleles,
    flair_fusion,
    flair_diffExp,
    flair_diffSplice,
    flair_spliceevents,
)

class MisplacedGlobalOptionAction(argparse.Action):
    """Report how to correct a global option given after the subcommand.  The value
    of an option that takes one is repeated in the correction, since the correction
    without it would be read as taking the subcommand name for the value."""

    def __call__(self, parser, namespace, values, option_string=None):
        subcommand = parser.prog.split()[-1]
        given = f"{option_string} {values}" if values else option_string
        parser.error(f"{option_string} is a global option and must precede the subcommand:\n"
                     f"    flair {given} {subcommand} ...")

def _global_options(parser):
    "the global options to trap, as option string and the action the global parser has for it"
    return tuple((opt, action) for opt, action in parser._option_string_actions.items()
                 if opt not in ('-h', '--help'))

def _trap_misplaced_global_options(parser, subparsers):
    """Global options are only accepted before the subcommand.  Give each subparser
    a hidden copy of each one so that using it there explains the correct form
    instead of argparse's bare `unrecognized arguments'.  A flag takes no value
    here either, so the argument after it is still reported as the subcommand's."""
    for subparser in subparsers.choices.values():
        for opt, action in _global_options(parser):
            subparser.add_argument(opt, nargs=(0 if action.nargs == 0 else '?'),
                                   action=MisplacedGlobalOptionAction,
                                   default=argparse.SUPPRESS, help=argparse.SUPPRESS)

def _set_help_formatter(parser, subparsers):
    """Every parser formats help with FlairHelpFormatter, set here rather than in
    each subcommand's add_subparser"""
    parser.formatter_class = FlairHelpFormatter
    for subparser in subparsers.choices.values():
        subparser.formatter_class = FlairHelpFormatter

def build_parser():
    """Build the whole flair parser, every subcommand included.  This is the one
    entry point used both to run commands and to generate the command line
    documentation."""
    desc = """Run a FLAIR module.  This program is the main entry point
    for running FLAIR analysis."""
    parser = RetiredCheckParser(prog='flair', description=desc)
    parser.add_argument('--version', action='version', version='FLAIR ' + VERSION,
                        help="print FLAIR version")
    parser.add_argument('-z', '--compress', action='store_true',
                        help="gzip the text files this run writes, adding .gz to their "
                             "names, FASTA with bgzip so it can be indexed.  BAM, bigBed "
                             "and PDFs are not gzipped, nor are intermediate files.  "
                             "Compressed input is read with or without this option, "
                             "detected by file extension")
    loggingOps.addCmdOptions(parser, defaultLevel=logging.INFO)
    subparsers = parser.add_subparsers(title="subcommands", dest="subcommand",
                                       metavar="<subcommand>", required=True)
    for module in SUBCOMMAND_MODULES:
        module.add_subparser(subparsers)
    _trap_misplaced_global_options(parser, subparsers)
    _set_help_formatter(parser, subparsers)
    return parser

def run_subcommand(args):
    log_run_start()
    start_time = time.time()
    args.entry(args)
    elapsed = time.time() - start_time
    logging.info(f"flair {args.subcommand} took " + str(timedelta(seconds=round(elapsed))))

def main():
    args = cli.parseArgsWithLogging(build_parser())
    configure_logging()
    install_stop_handlers()
    set_output_compression(args.compress)
    with cli.ErrorHandler():
        run_subcommand(args)


if __name__ == '__main__':
    main()
