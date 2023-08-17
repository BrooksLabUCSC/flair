#! /usr/bin/env python3

import sys
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import time
from multiprocessing import Pool
from flair_align import align
from flair_correct import correct
from flair_collapse import collapse
from flair_quantify import quantify

def main():
	path = '/'.join(os.path.realpath(__file__).split('/')[:-1])+'/'
	globals()['path'] = path
	if len(sys.argv) < 2:
		sys.stderr.write('usage: flair <mode> --help \n')
		sys.stderr.write('modes: align, correct, collapse, quantify, diffExp, diffSplice\n')
		sys.stderr.write('Multiple modules can be run when specified using numbers, e.g.:\n')
		sys.stderr.write('flair 1234 ...\n')
		sys.exit(1)
	else:
		mode = sys.argv[1].lower()

	aligned_reads, corrected_reads, isoforms, isoform_sequences, counts_matrix = [0]*5
	start_time = time.time()
	last_time = start_time
	if mode == 'align' or '1' in mode:
		print(f"Starting align...", flush=True)
		aligned_reads = align()
		cur_time = time.time()
		print(f"Flair align took {int((cur_time - last_time)/60)} minutes and {int((cur_time - last_time))%60} seconds", flush=True)
		last_time = cur_time

	if mode == 'correct' or '2' in mode:
		print(f"Starting correct...", flush=True)
		if aligned_reads:
			status = correct(aligned_reads=aligned_reads)
		else:
			status = correct()
		if status == 1:
			sys.exit(1)
		else:
			corrected_reads = status
		cur_time = time.time()
		print(f"Flair correct took {int((cur_time - last_time)/60)} minutes and {int((cur_time - last_time))%60} seconds", flush=True)
		last_time = cur_time

	if mode == 'collapse' or ('3' in mode and '3.5' not in mode):
		print(f"Starting collapse...", flush=True)
		if corrected_reads:
			status = collapse(corrected_reads=corrected_reads)
		else:
			status = collapse()
		if status == 1:
			sys.exit(1)
		else:
			isoforms, isoform_sequences = status
		cur_time = time.time()
		print(f"Flair collapse took {int((cur_time - last_time)/60)} minutes and {int((cur_time - last_time))%60} seconds", flush=True)
		last_time = cur_time


	if mode == 'quantify' or '4' in mode:
		print(f"Starting quantify...", flush=True)
		if isoform_sequences:
			status = quantify(isoform_sequences=isoform_sequences)
		else:
			status = quantify()
		if status == 1:
			sys.exit(1)
		else:
			counts_matrix = status
		cur_time = time.time()
		print(f"Flair quantify took {int((cur_time - last_time)/60)} minutes and {int((cur_time - last_time))%60} seconds", flush=True)
		last_time = cur_time

	
	if mode == '--version':
		sys.stderr.write('FLAIR v2.0.0\n')
		sys.exit(0)
	print(f"Flair took {int((cur_time - start_time)/60)} minutes and {int((cur_time - start_time))%60} seconds and finished without issues. If you do not see this message, please check your run.", flush=True)


if __name__ == '__main__':
    main()
