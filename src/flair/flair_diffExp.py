#! /usr/bin/env python3

import sys
import argparse
import subprocess
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'

# TODO: this is a complex interface to runDE, merge the two

def diffExp(counts_matrix=''):
	parser = argparse.ArgumentParser(description='flair-diffExp parse options',
		usage='flair diffExp -q counts_matrix.tsv --out_dir out_dir [options]')
#	parser.add_argument('diffExp')
	required = parser.add_argument_group('required named arguments')
	if not counts_matrix:
		required.add_argument('-q', '--counts_matrix', action='store', dest='q',
			type=str, required=True, help='Tab-delimited isoform count matrix from flair quantify module.')
	required.add_argument('-o', '--out_dir', action='store', dest='o',
		type=str, required=True, help='Output directory for tables and plots.')
	parser.add_argument('-t', '--threads', action='store', dest='t',
		type=int, required=False, default=4, help='Number of threads for parallel DRIMSeq.')
	parser.add_argument('-e', '--exp_thresh', action='store', dest='e', type=int, required=False,
		default=10, help='''Read count expression threshold. Isoforms in which
		both conditions contain fewer than E reads are filtered out (Default E=10)''')
	parser.add_argument('-of', '--out_dir_force', action='store_true', dest='of',
		required=False, help='''Specify this argument to force overwriting of files in
		an existing output directory''')
	args, unknown = parser.parse_known_args()
	if unknown:
		sys.stderr.write('DiffExp unrecognized arguments: {}\n'.format(' '.join(unknown)))
		if not counts_matrix:
			return 1
	if counts_matrix:
		args.q = counts_matrix
		args.o+'.diffExp'

	if not os.path.exists(args.q):
		sys.stderr.write('Counts matrix file path does not exist\n')
		return 1

	runDE = 'deFLAIR.py'
#	DEcommand = [sys.executable, '-W ignore', runDE, '--filter', str(args.e), '--threads',
	DEcommand = [runDE, '--filter', str(args.e), '--threads',
		str(args.t), '--outDir', args.o, '--matrix', args.q]
	if args.of:
		DEcommand += ['-of']

	if subprocess.call(DEcommand):
		return 1
	return

if __name__ == "__main__":
    diffExp()

