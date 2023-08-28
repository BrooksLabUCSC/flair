#! /usr/bin/env python3

import sys
import argparse
import os
import pipettor
os.environ['OPENBLAS_NUM_THREADS'] = '1'
from multiprocessing import Pool
from bam2Bed12 import bam2Bed12

# note to self: select interpreter for conda
def align():
	parser = argparse.ArgumentParser(description='flair-align parse options',
		usage='flair_align -g genome.fa -r <reads.fq>|<reads.fa> [options]')
	required = parser.add_argument_group('required named arguments')
	atleastone = parser.add_argument_group('Either one of the following arguments is required')
	required.add_argument('-r', '--reads', nargs='+', type=str, required=True, 
		help='FastA/FastQ files of raw reads')
	atleastone.add_argument('-g', '--genome', type=str, 
		help='FastA of reference genome, can be minimap2 indexed')
	atleastone.add_argument('--mm_index', type=str, default='', 
		help='minimap2 index .mmi file')
	parser.add_argument('-o', '--output', default='flair.aligned',
		help='output file name base (default: flair.aligned)')
	parser.add_argument('-t', '--threads', type=int, default=4,
		help='minimap2 number of threads (4)')
	parser.add_argument('--junction_bed', default='',
		help='annotated isoforms/junctions bed file for splice site-guided minimap2 genomic alignment')
	parser.add_argument('--nvrna', action='store_true', default=False,
		help='specify this flag to use native-RNA specific alignment parameters for minimap2')
	parser.add_argument('--quality', type=int, default=1,
		help='minimum MAPQ of read alignment to the genome (1)')
	parser.add_argument('-N', '--keep_supplementary', type=int, default=0,
        help='retain at most INT secondary alignments from minimap2 alignment (default 0)')
	parser.add_argument('--quiet', default=False, action='store_true', dest='quiet',
		help='''Suppress minimap progress statements from being printed''')
	
	no_arguments_passed = len(sys.argv) == 1
	if no_arguments_passed:
		parser.print_help()
		sys.exit(1)
	if 'align' in sys.argv:
		sys.argv.remove('align')
	args, unknown = parser.parse_known_args()
	if unknown and not args.quiet:
		sys.stderr.write('Align unrecognized arguments: {}\n'.format(' '.join(unknown)))

	# do we have multiple inputs?
	if ',' in args.reads[0]:
		args.reads = args.reads[0].split(',')
	for rfile in args.reads:
		if not os.path.exists(rfile):
			sys.stderr.write(f'Check that read file {rfile} exists\n')
			sys.exit(1)

	# define outputs
	bamout = args.output+'.bam'
	bedout = args.output+'.bed'

	# minimap
	mm2_cmd = ['minimap2', '-ax', 'splice', '-t', str(args.threads), args.genome]+args.reads
	if args.mm_index:
		mm2_cmd[5] = args.mm_index
	if args.nvrna:
		mm2_cmd[3:3] = ['-uf', '-k14']
	if args.junction_bed:
		mm2_cmd[3:3] = ['--junc-bed', args.junction_bed]
	if str(args.keep_supplementary) != '0':
		mm2_cmd[3:3] = ['-N', str(args.keep_supplementary)]
	else:
		mm2_cmd[3:3] = ['--secondary=no']
	mm2_cmd = tuple(mm2_cmd)
	
	# samtools; the dash at the end means STDIN
	samtools_filter_cmd = ('samtools', 'view', '-q', str(args.quality), '-h', '-')
	samtools_sort_cmd = ('samtools', 'sort', '-@', str(args.threads), '-o', bamout, '-')
	samtools_index_cmd = ('samtools', 'index', bamout)
	if not args.quiet:
		print('flair align minimap cmd:', mm2_cmd, file=sys.stderr)
		print('flair align samtools filter cmd:', samtools_filter_cmd, file=sys.stderr)
		print('flair align samtools sort cmd:', samtools_sort_cmd, file=sys.stderr)
		print('flair align samtools index cmd:', samtools_index_cmd, file=sys.stderr)
	pipettor.run([mm2_cmd, samtools_filter_cmd, samtools_sort_cmd])
	pipettor.run([samtools_index_cmd])
	
	bam2Bed12(bamout, bedout, args.keep_supplementary)
	return bedout

if __name__ == "__main__":
	align()

