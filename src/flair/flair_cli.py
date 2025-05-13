#! /usr/bin/env python3

import sys
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import time
from flair import VERSION, set_unix_path
from flair import flair_align
from flair.flair_correct import correct
from flair.flair_collapse import collapse
from flair.flair_quantify import quantify
from flair.flair_combine import combine
from flair.flair_straightfrombam import collapsefrombam

def help():
    "temporary help until switched to argparse"
    help = f"""usage: flair module --help
module: align, correct, collapse, quantify

Version {VERSION}
"""
    print(help, file=sys.stderr)


def main():
    set_unix_path()

    # FIXME: what is this for?
    path = '/'.join(os.path.realpath(__file__).split('/')[:-1])+'/'
    globals()['path'] = path
    if len(sys.argv) < 2:
        help()
        sys.exit(1)
    else:
        mode = sys.argv[1].lower()
    # remove mode from the sys arguments
    sys.argv.pop(1)

    if mode == '--help':
        help()
        sys.exit(0)

    if mode == '--version':
        sys.stderr.write(f'FLAIR {VERSION}\n')
        sys.exit(0)

    aligned_reads, corrected_reads, isoforms, isoform_sequences, counts_matrix = [0]*5
    cur_time = None
    start_time = time.time()
    last_time = start_time
    if mode == 'align':
        print(f"Starting align...", flush=True)
        aligned_reads = flair_align.align()
        cur_time = time.time()
        print(f"Flair align took {int((cur_time - last_time)/60)} minutes and {int((cur_time - last_time))%60} seconds", flush=True)
        last_time = cur_time

    if mode == 'correct':
        print(f"Starting correct...", flush=True)
        if aligned_reads:
            corrected_reads = correct(aligned_reads=aligned_reads)
        else:
            corrected_reads = correct()
        cur_time = time.time()
        print(f"Flair correct took {int((cur_time - last_time)/60)} minutes and {int((cur_time - last_time))%60} seconds", flush=True)
        last_time = cur_time

    if mode == 'collapse':
        print(f"Starting collapse...", flush=True)
        if corrected_reads:
            [isoforms, isoform_sequences] = collapse(corrected_reads=corrected_reads)
        else:
            [isoforms, isoform_sequences] = collapse()
        cur_time = time.time()
        print(f"Flair collapse took {int((cur_time - last_time)/60)} minutes and {int((cur_time - last_time))%60} seconds", flush=True)
        last_time = cur_time

    if mode == 'transcriptome':
        print(f"Starting transcriptome generation...", flush=True)
        collapsefrombam()
        cur_time = time.time()
        print(f"Flair transcriptome generation took {int((cur_time - last_time)/60)} minutes and {int((cur_time - last_time))%60} seconds", flush=True)
        last_time = cur_time

    if mode == 'quantify':
        print(f"Starting quantify...", flush=True)
        if isoform_sequences:
            counts_matrix = quantify(isoform_sequences=isoform_sequences)
        else:
            counts_matrix = quantify()
        cur_time = time.time()
        print(f"Flair quantify took {int((cur_time - last_time)/60)} minutes and {int((cur_time - last_time))%60} seconds", flush=True)
        last_time = cur_time

    if mode == 'combine':
        print(f"Starting combine...", flush=True)
        flair_combine.combine()
        cur_time = time.time()
        print(
                f"Flair combine took {int((cur_time - last_time) / 60)} minutes and {int((cur_time - last_time)) % 60} seconds",
                flush=True)

    if mode == 'collapse-range':
        sys.stderr.write('ERROR: This version of flair does not support collapse-range.\n')
        sys.stderr.write('We are in the process of making a better batch mode for flair.\n')
        sys.exit(1)

    if mode in ['diffexp', 'diffsplice', ]:
        sys.stderr.write('ERROR: This version of flair does not support diffExp and diffSplice.\n')
        sys.stderr.write('You can run flair_diffexp and flair_diffsplice as separate programs.\n')
        sys.exit(1)

    cur_time = time.time()
    print(f"\nFlair took {int((cur_time - start_time)/60)} minutes and {int((cur_time - start_time))%60} seconds and finished without issues.\n\nFLAIR HAS FINISHED\n\n", flush=True)


if __name__ == '__main__':
    main()
