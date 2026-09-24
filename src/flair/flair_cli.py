#! /usr/bin/env python3
import time
from datetime import timedelta
import argparse
import logging
from flair import VERSION, set_unix_path
from flair.pycbio.sys import cli, loggingOps
from flair import (flair_align, flair_allelotyping, flair_combine, flair_diffExp,
                   flair_diffSplice, flair_fusion, flair_isoalleles, flair_quantify,
                   flair_transcriptome, flair_variantquant)

# Modules in pipeline order.  Each one adds and then owns its own subparser and
# records its entry point with set_defaults(entry=), so no parser is filled in by
# one module and parsed by another.
#
# These are imported eagerly.  Importing all ten costs about 0.4 seconds, most of
# it numpy in flair_diffExp, which is not worth the pre-parsing needed to load
# them one at a time.
SUBCOMMAND_MODULES = (
    flair_align,
    flair_transcriptome,
    flair_quantify,
    flair_combine,
    flair_variantquant,
    flair_fusion,
    flair_diffExp,
    flair_diffSplice,
    flair_allelotyping,
    flair_isoalleles,
)

class MisplacedGlobalOptionAction(argparse.Action):
    """Report how to correct a global option given after the subcommand."""

    def __call__(self, parser, namespace, values, option_string=None):
        subcommand = parser.prog.split()[-1]
        parser.error(f"{option_string} is a global option and must precede the subcommand:\n"
                     f"    flair {option_string} {subcommand} ...")

def _global_option_strings(parser):
    return tuple(opt for opt in parser._option_string_actions
                 if opt not in ('-h', '--help'))

def _trap_misplaced_global_options(parser, subparsers):
    """Global options are only accepted before the subcommand.  Give each subparser
    a hidden copy of each one so that using it there explains the correct form
    instead of argparse's bare `unrecognized arguments'."""
    for subparser in subparsers.choices.values():
        for opt in _global_option_strings(parser):
            subparser.add_argument(opt, nargs='?', action=MisplacedGlobalOptionAction,
                                   default=argparse.SUPPRESS, help=argparse.SUPPRESS)

def build_parser():
    """Build the whole flair parser, every subcommand included.  This is the one
    entry point used both to run commands and to generate the command line
    documentation."""
    desc = """Run a FLAIR module.  This program is the main entry point
    for running FLAIR analysis."""
    parser = argparse.ArgumentParser(prog='flair', description=desc)
    parser.add_argument('--version', action='version', version='FLAIR ' + VERSION,
                        help="print FLAIR version")
    loggingOps.addCmdOptions(parser, defaultLevel=logging.WARNING)
    subparsers = parser.add_subparsers(title="subcommands", dest="subcommand",
                                       metavar="<subcommand>", required=True)
    for module in SUBCOMMAND_MODULES:
        module.add_subparser(subparsers)
    _trap_misplaced_global_options(parser, subparsers)
    return parser

def run_subcommand(args):
    start_time = time.time()
    args.entry(args)
    elapsed = time.time() - start_time
    logging.info(f"Flair {args.subcommand} took " + str(timedelta(seconds=round(elapsed))))

def main():
    set_unix_path()
    args = cli.parseArgsWithLogging(build_parser())
    with cli.ErrorHandler():
        run_subcommand(args)


if __name__ == '__main__':
    main()
