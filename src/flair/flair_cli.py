#! /usr/bin/env python3

import sys
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import time
from flair import VERSION, set_unix_path

def help():
    "temporary help until switched to argparse"
    help = f"""usage: flair module --help
module: align, correct, collapse, quantify

Version {VERSION}
"""
    print(help, file=sys.stderr)


def main():
    set_unix_path()

    path = '/'.join(os.path.realpath(__file__).split('/')[:-1])+'/'
    globals()['path'] = path
    if len(sys.argv) < 2:
        help()
        sys.exit(1)
    # remove mode from the sys arguments
    mode = sys.argv[1].lower()
    sys.argv.pop(1)

    if mode == '--help':
        help()
        sys.exit(0)

    if mode == '--version':
        sys.stderr.write(f'FLAIR {VERSION}\n')
        sys.exit(0)

    aligned_reads, corrected_reads, isoforms, isoform_sequences, counts_matrix = [0]*5
    start_time = time.time()
    if mode == 'align':
        from flair import flair_align
        aligned_reads = flair_align.align()

    elif mode == 'correct':
        from flair.flair_correct import correct
        if aligned_reads:
            corrected_reads = correct(aligned_reads=aligned_reads)
        else:
            corrected_reads = correct()

    elif mode == 'collapse':
        from flair.flair_collapse import collapse
        if corrected_reads:
            collapse(corrected_reads=corrected_reads)
        else:
            collapse()

    elif mode == 'quantify':
        from flair.flair_quantify import quantify
        if isoform_sequences:
            counts_matrix = quantify(isoform_sequences=isoform_sequences)
        else:
            counts_matrix = quantify()

    elif mode == 'combine':
        from flair import flair_combine
        flair_combine.combine()

    elif mode == 'collapse-range':
        sys.stderr.write('ERROR: This version of flair does not support collapse-range.\n')
        sys.stderr.write('We are in the process of making a better batch mode for flair.\n')
        sys.exit(1)

    elif mode == 'diffexp':
        from flair import flair_diffExp
        flair_diffExp.diffExp()
    elif mode == 'diffsplice':
        from flair import flair_diffSplice
        flair_diffSplice.diffSplice()
    else:
        print(f"Error: invalid flair mode`{mode}'", file=sys.stderr, flush=True)
        exit(1)

    end_time = time.time()
    print(f"Flair {mode} took {int((end_time - start_time)/60)} minutes and {int((end_time - start_time))%60} seconds and finished without issues.\n\nFLAIR HAS FINISHED", flush=True)


if __name__ == '__main__':
    main()
