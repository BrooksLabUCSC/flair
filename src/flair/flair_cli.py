#! /usr/bin/env python3
import sys
import time
from datetime import timedelta
import argparse
import logging
from flair import VERSION, set_unix_path
from flair.pycbio.sys import cli

VALID_MODULES = ('align', 'correct', 'transcriptome', 'collapse', 'quantify',
                 'combine', 'variants', 'fusion', 'diffexp', 'diffsplice')

def parse_args():
    """Argument parsing a module name for this and then returns the remaining arguments
    to pass on to
    Special handling is done for the logging options that are added to parser

    These are split out so they can be passed to modules when they passed to
    modules as they support it.
    """

    desc = '''Run a FLAIR module.  This program is the main entry point
    for running FLAIR analysis.'''
    parser = cli.ArgumentParserExtras(description=desc)
    parser.add_argument('--version', action='version', version='FLAIR ' + VERSION,
                        help="print FLAIR version")
    parser.add_argument("module", choices=VALID_MODULES, type=str.lower,
                        help="name of module to run")
    parser.add_argument('module_args', nargs=argparse.REMAINDER,
                        help="arguments to module")
    return parser.parse_opts_args()

def move_opts_to_argv(opts, module_argv):
    """move options back to module_argv for modules that have been converted
    to pycbio CLI."""
    nargv = []
    for n, v in opts.items():
        if n != 'version':
            nargv.extend([f'--{n}', v])
    return nargv + module_argv

def flair_module_run(opts, module, module_argv):  # noqa: C901
    start_time = time.time()
    sys.argv = [sys.argv[1]] + module_argv

    # delay import modules until needed to speed startup
    if module == 'align':
        from flair import flair_align
        flair_align.align()
    elif module == 'correct':
        from flair.flair_correct import correct
        correct()
    elif module == 'transcriptome':
        from flair.flair_transcriptome import collapsefrombam
        collapsefrombam()
    elif module == 'collapse':
        from flair.flair_collapse import collapse
        collapse()
    elif module == 'quantify':
        from flair.flair_quantify import quantify
        quantify()
    elif module == 'combine':
        from flair import flair_combine
        flair_combine.combine()
    elif module == 'variants':
        from flair.flair_variants import getvariants
        getvariants()
    elif module == 'fusion':
        from flair.flair_fusion import detectfusions
        detectfusions()
    elif module == 'diffexp':
        from flair import flair_diffExp
        flair_diffExp.diffExp()
    elif module == 'diffsplice':
        from flair import flair_diffSplice
        flair_diffSplice.diffSplice()

    elapsed = time.time() - start_time
    logging.info(f"Flair {module} took " + str(timedelta(seconds=round(elapsed))))

def main():
    set_unix_path()
    opts, args = parse_args()
    with cli.ErrorHandler():
        flair_module_run(opts, args.module, args.module_args)


if __name__ == '__main__':
    main()
