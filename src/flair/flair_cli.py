#! /usr/bin/env python3
import sys
import time
from datetime import timedelta
import argparse
import logging
from flair import VERSION, set_unix_path
from flair.pycbio.sys import cli

VALID_MODULES = ('align', 'transcriptome', 'quantify', 'combine',
                 'variantquant', 'fusion', 'diffexp', 'diffsplice', 'alleles', 'isoalleles')

def parse_args():
    """Argument parsing a module name for this and then returns the remaining arguments
    to pass on to
    Special handling is done for the logging options that are added to parser

    These are split out so they can be passed to modules when they passed to
    modules as they support it.
    """

    desc = '''Run a FLAIR module.  This program is the main entry point
    for running FLAIR analysis.'''
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('--version', action='version', version='FLAIR ' + VERSION,
                        help="print FLAIR version")
    parser.add_argument("module", choices=VALID_MODULES, type=str.lower,
                        help="name of module to run")
    parser.add_argument('module_args', nargs=argparse.REMAINDER,
                        help="arguments to module")
    # allow changing the default logging level here
    return cli.splitOptionsArgs(parser, cli.parseArgsWithLogging(parser, defaultLevel=logging.WARNING))

def flair_module_run(opts, module, module_argv):  # noqa: C901
    start_time = time.time()
    sys.argv = [sys.argv[1]] + module_argv

    # delay import modules until needed to speed startup
    if module == 'align':
        from flair import flair_align
        flair_align.align()
    elif module == 'transcriptome':
        from flair.flair_transcriptome import flair_transcriptome
        flair_transcriptome()
    elif module == 'quantify':
        from flair.flair_quantify import quantify
        quantify()
    elif module == 'combine':
        from flair import flair_combine
        flair_combine.combine()
    elif module == 'alleles':
        from flair.flair_allelotyping import getvariants
        getvariants()
    elif module == 'isoalleles':
        from flair.flair_isoalleles import getvariants
        getvariants()
    elif module == 'variantquant':
        from flair.flair_variantquant import quantvarpos
        quantvarpos()
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
