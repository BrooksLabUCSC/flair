#! /usr/bin/env python3

import sys
import re
import argparse
import subprocess
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import tempfile
import glob
from multiprocessing import Pool


def align():
	parser = argparse.ArgumentParser(description='flair-align parse options',
		usage='flair align -g genome.fa -r <reads.fq>|<reads.fa> [options]')
	parser.add_argument('align')
	required = parser.add_argument_group('required named arguments')
	atleastone = parser.add_argument_group('Either one of the following arguments is required')
	required.add_argument('-r', '--reads', action='store', dest='r',
		nargs='+', type=str, required=True, help='FastA/FastQ files of raw reads')
	atleastone.add_argument('-g', '--genome', action='store', dest='g',
		type=str, help='FastA of reference genome, can be minimap2 indexed')
	atleastone.add_argument('--mm_index', action='store', dest='mm_index', type=str, default='',
		help='minimap2 index .mmi file')
	parser.add_argument('-o', '--output', action='store', dest='o', default='flair.aligned',
		help='output file name base (default: flair.aligned)')
	parser.add_argument('-t', '--threads', type=int,
		action='store', dest='t', default=4, help='minimap2 number of threads (4)')
	parser.add_argument('--junction_bed', action='store', dest='junction_bed', default='',
		help='annotated isoforms/junctions bed file for splice site-guided minimap2 genomic alignment')
	parser.add_argument('--nvrna', action='store_true', dest='n', default=False,
		help='specify this flag to use native-RNA specific alignment parameters for minimap2')
	parser.add_argument('--quality', type=int, action='store', dest='quality', default=1,
		help='minimum MAPQ of read alignment to the genome (1)')
	parser.add_argument('-N', type=int, action='store', dest='N', default=0,
		help='retain at most INT secondary alignments from minimap2 alignment (0)')
	parser.add_argument('--quiet', default=False, action='store_true', dest='quiet',
			help='''Suppress progress statements from being printed''')
	args, unknown = parser.parse_known_args()
	if unknown and not args.quiet:
		sys.stderr.write('Align unrecognized arguments: {}\n'.format(' '.join(unknown)))

	if ',' in args.r[0]:
		args.r = args.r[0].split(',')
	for rfile in args.r:
	#for i in range(len(args.r)):
		if not os.path.exists(rfile):
			sys.stderr.write(f'Check that read file {rfile} exists\n')
			return 1

	mm2_command = ['minimap2', '-ax', 'splice', '-t', str(args.t), args.g]+args.r
	if args.mm_index:
		mm2_command[5] = args.mm_index
	if args.n:
		mm2_command[3:3] = ['-uf', '-k14']
	if args.junction_bed:
		mm2_command[3:3] = ['--junc-bed', args.junction_bed]
	if str(args.N) != '0':
		mm2_command[3:3] = ['-N', str(args.N)]
	else:
		mm2_command[3:3] = ['--secondary=no']

	try:
		if args.quiet:
			if subprocess.call(mm2_command, stdout=open(args.o+'.sam', 'w'),
			stderr=open(args.o+'.mm2_stderr', 'w')):
				return 1
		elif subprocess.call(mm2_command, stdout=open(args.o+'.sam', 'w')):
			return 1
	except:
		sys.stderr.write('Possible minimap2 error, specify executable path with -m?\n')
		if args.quiet:
			sys.stderr.write('Check {}\n'.format(args.o+'.mm2_stderr'))
		return 1

	if args.quality != 0:
		if subprocess.call(['samtools', 'view', '-q', str(args.quality), '-h', '-S', args.o+'.sam'],
		stdout=open(args.o+'.q.sam', 'w'), stderr=open(args.o+'.samtools_stderr', 'w')):
			sys.stderr.write('Possible issue with samtools, see {}\n'.format(args.o+'.samtools_stderr'))
			return 1

		subprocess.check_call(['mv', args.o+'.q.sam', args.o+'.sam'])
		subprocess.check_call(['rm', args.o+'.samtools_stderr'])

	stderrfile = args.o + '.bam.stderr'
	if subprocess.call(['samtools', 'sort', '-@', str(args.t), args.o+'.sam', '-o', args.o+'.bam'],
		stderr=open(stderrfile, 'w')):
		sys.stderr.write(f'Samtools issue with sorting minimap2 sam, see {stderrfile}\n')
		return 1
	subprocess.check_call(['rm', stderrfile, args.o+'.sam'])

	subprocess.check_call(['samtools', 'index', args.o+'.bam'])

	bed_conversion_cmd = [sys.executable, path+'bam2Bed12.py', '-i', args.o+'.bam']
	if args.N != 0:
		bed_conversion_cmd += ['--keep_supplementary']
	subprocess.check_call(bed_conversion_cmd, stdout=open(args.o+'.bed', 'w'))

	return args.o+'.bed'


def correct(aligned_reads=''):
	parser = argparse.ArgumentParser(description='flair-correct parse options',
		usage='flair correct -q query.bed12 [-f annotation.gtf]v[-j introns.tab] -g genome.fa [options]')
	parser.add_argument('correct')
	required = parser.add_argument_group('required named arguments')
	atleastone = parser.add_argument_group('at least one of the following arguments is required')
	if not aligned_reads:
		required.add_argument('-q', '--query', type=str, default='', required=True,
			action='store', dest='q', help='uncorrected bed12 file')
	required.add_argument('-g', '--genome', action='store', dest='g',
		type=str, required=True, help='FastA of reference genome')
	atleastone.add_argument('-j', '--shortread', action='store', dest='j', type=str, default='',
		help='bed format splice junctions from short-read sequencing')
	atleastone.add_argument('-f', '--gtf', default='',
		action='store', dest='f', help='GTF annotation file')
	parser.add_argument('-o', '--output',
		action='store', dest='o', default='flair', help='output name base (default: flair)')
	parser.add_argument('-t', '--threads', type=int, action='store', dest='t', default=4,
		help='splice site correction script number of threads (4)')
	parser.add_argument('--nvrna', action='store_true', dest='n', default=False,
		help='''specify this flag to make the strand of a read consistent with the annotation during correction''')
	parser.add_argument('-w', '--ss_window', type=int, action='store', dest='w', default=15,
		help='window size for correcting splice sites (15)')
	parser.add_argument('--print_check',
		action='store_true', dest='p', default=False, help='Print err.txt with step checking.')
	args, unknown = parser.parse_known_args()
	if unknown:
		sys.stderr.write('Correct unrecognized arguments: {}\n'.format(' '.join(unknown)))
		if not aligned_reads:
			return 1

	if aligned_reads:
		args.q = aligned_reads
	if not os.path.exists(args.q):
		sys.stderr.write(f'File {args.q} does not appear to exist, please check.\n')
		return 1
	if not os.path.exists(args.g):
		sys.stderr.write(f'File {args.g} does not appear to exist, please check.\n')
		return 1
	if not args.j and not args.f:
		sys.stderr.write('Please specify at least one of the -f or -j arguments for correction\n')
		return 1
	if args.j and not os.path.exists(args.j):
		sys.stderr.write(f'File {args.j} does not appear to exist, please check.\n')
		return 1
	if args.f and not os.path.exists(args.f):
		sys.stderr.write(f'File {args.f} does not appear to exist, please check.\n')
		return 1

	correction_cmd = [sys.executable, path+'ssCorrect.py', '-i', args.q,
			'-w', str(args.w), '-p', str(args.t), '-o', args.o, '--progress', '-f', args.g]
	if not args.n:
		correction_cmd += ['--correctStrand']
	if args.j:
		correction_cmd += ['-j', args.j]
	if args.f:
		correction_cmd += ['-g', args.f]
	if args.p:
		correction_cmd += ['--print_check']

	if subprocess.call(correction_cmd):
		printcmd = ' '.join(correction_cmd)
		sys.stderr.write(f'Correction command did not exit with success status:\n{printcmd}\n\n')

	return args.o+'_all_corrected.bed'


def collapse_range(run_id='', corrected_reads='', aligned_reads=''):
	parser = argparse.ArgumentParser(description='flair-collapse parse options',
		usage='flair collapse-range -g genome.fa -r reads.bam -q <query.bed> [options]')
	parser.add_argument('collapse')
	required = parser.add_argument_group('required named arguments')
	required.add_argument('-r', '--reads', action='store', dest='r', nargs='+',
		type=str, required=True, help='bam file(s) of the aligned reads')
	if not corrected_reads:
		required.add_argument('-q', '--query', type=str, default='', required=True,
			action='store', dest='q', help='bed file of aligned/corrected reads')
	required.add_argument('-g', '--genome', action='store', dest='g',
		type=str, required=True, help='FastA of reference genome')
	parser.add_argument('-f', '--gtf', default='', action='store', dest='f',
		help='GTF annotation file, used for renaming FLAIR isoforms to annotated isoforms and adjusting TSS/TESs')
	parser.add_argument('--mmi', type=str, default='',
		action='store', dest='mmi', help='minimap2 indexed genome')
	parser.add_argument('-t', '--threads', type=int,
		action='store', dest='t', default=4, help='minimap2 number of threads (4)')
	parser.add_argument('-p', '--promoters', action='store', dest='p', default='',
		help='promoter regions bed file to identify full-length reads')
	parser.add_argument('-w', '--end_window', type=int, default=100, action='store', dest='w',
		help='window size for comparing TSS/TES (100)')
	parser.add_argument('-s', '--support', type=int, default=3, action='store', dest='s',
		help='minimum number of supporting reads for an isoform (3)')
	parser.add_argument('--stringent', default=False, action='store_true', dest='stringent',
		help='''specify if all supporting reads need to be full-length
		(80%% coverage and spanning 25 bp of the first and last exons)''')
	parser.add_argument('-n', '--no_redundant', default='none', action='store', dest='n',
		help='''For each unique splice junction chain, report options include:
		none--best TSSs/TESs chosen for each unique set of splice junctions;
		longest--single TSS/TES chosen to maximize length;
		best_only--single most supported TSS/TES used in conjunction chosen (none)''')
	parser.add_argument('-i', '--isoformtss', default=False, action='store_true', dest='i',
		help='''when specified, TSS/TES for each isoform will be determined from supporting reads
		for individual isoforms (default: not specified, determined at the gene level)''')
	parser.add_argument('--max_ends', type=int, default=2, action='store', dest='max_ends',
		help='maximum number of TSS/TES picked per isoform (2)')
	parser.add_argument('--trust_ends', default=False, action='store_true', dest='trust_ends',
		help='specify if reads are generated from a long read method with minimal fragmentation')
	parser.add_argument('--filter', default='default', action='store', dest='filter',
		help='''Report options include:
		nosubset--any isoforms that are a proper set of another isoform are removed;
		default--subset isoforms are removed based on support;
		comprehensive--default set + all subset isoforms;
		ginormous--comprehensive set + single exon subset isoforms''')
	parser.add_argument('--quality', type=int, action='store', dest='quality', default=1,
		help='minimum MAPQ of read assignment to an isoform (1)')
	parser.add_argument('--keep_intermediate', default=False, action='store_true', dest='keep_intermediate',
		help='''specify if intermediate and temporary files are to be kept for debugging.
		Intermediate files include: promoter-supported reads file,
		read assignments to firstpass isoforms''')
	parser.add_argument('--generate_map', default=False, action='store_true', dest='generate_map',
		help='''specify this argument to generate a txt file of which reads are assigned to each isoform.
		(default: not specified)''')
	parser.add_argument('--quiet', default=False, action='store_true', dest='quiet',
		help='''Suppress progress statements from being printed''')
	parser.add_argument('--temp_dir', default='', action='store', dest='temp_dir',
		help='directory to put temporary files. use "./" to indicate current directory (default: python tempfile directory)')
	parser.add_argument('-o', '--output', default='flair.collapse',
		action='store', dest='o', help='output file name base for FLAIR isoforms (default: flair.collapse)')
	args, unknown = parser.parse_known_args()
	if unknown and not args.quiet:
		sys.stderr.write('Collapse-range unrecognized arguments: {}\n'.format(' '.join(unknown)))

	if corrected_reads:
		args.q = corrected_reads
		args.r = [aligned_reads[:-3]+'bam']

	if args.r[0][-3:] != 'bam':
		sys.stderr.write('Must provide genome alignment BAM with -r if range is specified\n')
		return 1

	if args.temp_dir == '':
		args.temp_dir = tempfile.NamedTemporaryFile().name+'/'
		if not args.quiet:
			sys.stderr.write('Writing temporary files to {}\n'.format(args.temp_dir))
	if not os.path.isdir(args.temp_dir): # make temporary directory
		if subprocess.call(['mkdir', args.temp_dir]):
			sys.stderr.write('Could not make temporary directory {}\n'.format(args.temp_dir))
			return 1
	if args.temp_dir[-1] != '/':
		args.temp_dir += '/'

	if args.q.endswith('psl'):
		sys.stderr.write('** Error. Flair no longer accepts PSL input. Please use psl_to_bed first.\n')
		return 1

	# partition the bed file into independent regions
	subprocess.check_call(['sort', '-k1,1', '-k2,2n', '--parallel='+str(args.t), args.q],
		stdout=open(args.temp_dir+run_id+'.sorted.bed', 'w'))
	if subprocess.call(['bedPartition', '-parallel='+str(args.t), args.temp_dir+run_id+'.sorted.bed', args.o+'.ranges.bed']):
		sys.stderr.write('''Make sure bedPartition (http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/)
		is an executable in your $PATH\n''')
		return 1
	ranges = []
	for line in open(args.o+'.ranges.bed'):
		line = line.rstrip().split('\t')
		ranges += [line[0]+':'+line[1]+'-'+line[2]]

	# index the bed file
	subprocess.check_call(['bgzip', args.temp_dir+run_id+'.sorted.bed'])
	if subprocess.call(['tabix', '-f', '--preset', 'bed', '--zero-based', args.temp_dir+run_id+'.sorted.bed.gz']):
		return 1

	# call collapse on all the ranges
	p = Pool(args.t)
	if 1 in p.map(collapse, ranges): # if a process failed
		return 1
	p.terminate()

	# consolidate all the isoforms from all the ranges
	subprocess.check_call([sys.executable, path+'consolidate_isoforms.py', args.temp_dir, run_id, args.o,
		'--generate_map', str(args.generate_map)])
	subprocess.check_call([sys.executable, path+'bed_to_sequence.py', args.o+'.isoforms.bed', args.f,
		args.o+'.isoforms.fa'])
	subprocess.check_call([sys.executable, path+'bed_to_gtf.py', args.o+'.isoforms.bed'],
		stdout=open(args.o+'.isoforms.gtf', 'w'))
	return args.o+'.isoforms.bed', args.o+'.isoforms.fa'


def collapse(genomic_range='', corrected_reads=''):
	parser = argparse.ArgumentParser(description='flair-collapse parse options',
		usage='''flair collapse -g genome.fa -q query.bed
		-r <reads.fq>/<reads.fa> [options]''')
	parser.add_argument('collapse')
	required = parser.add_argument_group('required named arguments')
	if not corrected_reads:
		required.add_argument('-q', '--query', type=str, default='', required=True,
			action='store', dest='q', help='bed file of aligned/corrected reads')
	required.add_argument('-g', '--genome', action='store', dest='g',
		type=str, required=True, help='FastA of reference genome')
	required.add_argument('-r', '--reads', action='store', dest='r', nargs='+',
		type=str, required=True, help='FastA/FastQ files of raw reads, can specify multiple files')
	parser.add_argument('-o', '--output', default='flair.collapse',
		action='store', dest='o', help='output file name base for FLAIR isoforms (default: flair.collapse)')
	parser.add_argument('-t', '--threads', type=int,
		action='store', dest='t', default=4, help='minimap2 number of threads (4)')
	parser.add_argument('-f', '--gtf', default='', action='store', dest='f',
		help='GTF annotation file, used for renaming FLAIR isoforms to annotated isoforms and adjusting TSS/TESs')
	parser.add_argument('--generate_map', default=False, action='store_true', dest='generate_map',
		help='''specify this argument to generate a txt file of read-isoform assignments
		(default: not specified)''')
	parser.add_argument('--annotation_reliant', default=False, action='store', dest='annotation_reliant',
		help='''specify transcript fasta that corresponds to transcripts in the gtf to run annotation-
		reliant flair collapse; to ask flair to make transcript sequences given the gtf and genome fa,
		type --annotation_reliant generate''')
	# supporting read assignment options
	parser.add_argument('-s', '--support', type=float, default=3.0, action='store', dest='s',
		help='''minimum number of supporting reads for an isoform;
		if s < 1, it will be treated as a percentage of expression of the gene (3)''')
	parser.add_argument('--stringent', default=False, action='store_true', dest='stringent',
		help='''specify if all supporting reads need to be full-length
		(80%% coverage and spanning 25 bp of the first and last exons)''')
	parser.add_argument('--check_splice', default=False, action='store_true', dest='check_splice',
		help='''enforce coverage of 4 out of 6 bp around each splice site and no
		insertions greater than 3 bp at the splice site''')
	parser.add_argument('--trust_ends', default=False, action='store_true', dest='trust_ends',
		help='specify if reads are generated from a long read method with minimal fragmentation')
	parser.add_argument('--quality', type=int, action='store', dest='quality', default=1,
		help='minimum MAPQ of read assignment to an isoform (1)')
	# variant options
	parser.add_argument('--longshot_bam', action='store', dest='longshot_bam', default='',
		help='bam from longshot containing haplotype information for each read')
	parser.add_argument('--longshot_vcf', action='store', dest='longshot_vcf', default='',
		help='vcf from longshot')
	# ends options
	parser.add_argument('-w', '--end_window', type=int, default=100, action='store', dest='w',
		help='window size for comparing TSS/TES (100)')
	parser.add_argument('-p', '--promoters', action='store', dest='p', default='',
		help='promoter regions bed file to identify full-length reads')
	parser.add_argument('--3prime_regions', action='store', dest='threeprime', default='',
		help='TES regions bed file to identify full-length reads')
	parser.add_argument('-n', '--no_redundant', default='none', action='store', dest='n',
		help='''For each unique splice junction chain, report options include:
		none--best TSSs/TESs chosen for each unique set of splice junctions;
		longest--single TSS/TES chosen to maximize length;
		best_only--single most supported TSS/TES used in conjunction chosen (none)''')
	parser.add_argument('-i', '--isoformtss', default=False, action='store_true', dest='i',
		help='''when specified, TSS/TES for each isoform will be determined from supporting reads
		for individual isoforms (default: not specified, determined at the gene level)''')
	parser.add_argument('--no_gtf_end_adjustment', default=False, action='store_true',
		dest='no_end_adjustment',
		help='''when specified, TSS/TES from the gtf provided with -f will not be used to adjust
		isoform TSSs/TESs each isoform will be determined from supporting reads''')
	parser.add_argument('--max_ends', type=int, default=2, action='store', dest='max_ends',
		help='maximum number of TSS/TES picked per isoform (2)')
	parser.add_argument('--filter', default='default', action='store', dest='filter',
		help='''Report options include:
		nosubset--any isoforms that are a proper set of another isoform are removed;
		default--subset isoforms are removed based on support;
		comprehensive--default set + all subset isoforms;
		ginormous--comprehensive set + single exon subset isoforms''')
	# other
	parser.add_argument('--temp_dir', default='', action='store', dest='temp_dir',
		help='directory for temporary files. use "./" to indicate current directory (default: python tempfile directory)')
	parser.add_argument('--keep_intermediate', default=False, action='store_true', dest='keep_intermediate',
		help='''specify if intermediate and temporary files are to be kept for debugging.
		Intermediate files include: promoter-supported reads file,
		read assignments to firstpass isoforms''')
	parser.add_argument('--fusion_dist', default=0, type=int, action='store', dest='fusion_dist',
			help='''minimium distance between separate read alignments on the same chromosome to be
			considered a fusion, otherwise no reads will be assumed to be fusions''')
	parser.add_argument('--mm2_args', action='store', dest='mm2_args',
		type=str, default=[], help='''additional minimap2 arguments when aligning reads first-pass transcripts;
		separate args by commas, e.g. --mm2_args=-I8g,--MD ''')
	parser.add_argument('--quiet', default=False, action='store_true', dest='quiet',
			help='''Suppress progress statements from being printed''')
	parser.add_argument('--annotated_bed', default=False, action='store', dest='annotated_bed',
		help='''annotation_reliant also requires a bedfile of annotated isoforms; if this isn't provided,
		flair collapse will generate the bedfile from the gtf. eventually this argument will be removed''')
	parser.add_argument('--range', default='', action='store', dest='range',
		help='''interval for which to collapse isoforms, formatted chromosome:coord1-coord2 or
		tab-delimited; if a range is specified, then the aligned reads bam must be specified with -r
		and the query must be a sorted, bgzip-ed bed file''')
	args, unknown = parser.parse_known_args()
	if unknown and not args.quiet:
		sys.stderr.write('Collapse unrecognized arguments: {}\n'.format(' '.join(unknown)))
		if not corrected_reads:
			return 1
	if corrected_reads:
		args.q = corrected_reads

	# housekeeping stuff
	tempfile_dir = tempfile.NamedTemporaryFile().name
	tempfile_name = tempfile_dir[tempfile_dir.rfind('/')+1:]+'.'
	if args.temp_dir == '':
		args.temp_dir = tempfile_dir+'/'
		if not args.quiet:
			sys.stderr.write('Writing temporary files to {}\t\n'.format(args.temp_dir))
	if not os.path.isdir(args.temp_dir): # make temporary directory
		if subprocess.call(['mkdir', args.temp_dir]):
			sys.stderr.write('Could not make temporary directory {}\n'.format(args.temp_dir))
			return 1
	if args.temp_dir[-1] != '/':
		args.temp_dir += '/'

	if genomic_range: # this module was called internally from collapse_range
		args.range = genomic_range
		args.o = args.temp_dir+run_id
		args.q = args.temp_dir+run_id+'.sorted.bed.gz'
		args.quiet = True

	args.quality = '0' if args.trust_ends else args.quality
	args.o += '.'
	min_reads = float(args.s) if float(args.s) >= 1 else 3
	if args.fusion_dist:
		args.trust_ends = True

	if ',' in args.r[0]:
		args.r = args.r[0].split(',')
	for rfile in args.r:
		if not os.path.exists(rfile):
			sys.stderr.write(f'Read file path does not exist: {rfile}\n')
			return 1
	if not os.path.exists(args.q):
		sys.stderr.write('Query file path does not exist: {}\n'.format(args.q))
		return 1
	if not os.path.exists(args.g):
		sys.stderr.write('Genome file path does not exist: {}\n'.format(args.q))
		return 1
	if os.stat(args.q).st_size == 0:
		sys.stderr.write('Query file is empty\n')
		return 1
	if float(args.s) < 1 and not args.f:
		sys.stderr.write('Provide gtf for gene grouping if -s is percentage of total gene expression\n')
		return 1

	if args.range:
		# subset out the read sequences and corrected reads corresponding to the specified range
		if '\t' in args.range:
			args.range = args.range.split('\t')
			args.range = args.range[0]+':'+args.range[1]+'-'+args.range[2]
		args.o += args.range+'.'
		if args.r[0][-3:] != 'bam':
			sys.stderr.write('Must provide genome alignment BAM with -r if range is specified\n')
			return 1
		bams = []
		for i in range(len(args.r)): # subset bam file for alignments within range
			bams += [args.temp_dir+tempfile_name+args.range+str(i)+'.bam']
			if subprocess.call(['samtools', 'view', '-h', args.r[i], args.range],
				stdout=open(bams[-1], 'w')):
				return 1
		args.r = []
		for i in range(len(bams)): # read sequences of the alignments within range
			args.r += [bams[i][:-3]+'fasta']
			subprocess.check_call(['samtools', 'fasta', bams[i]],
				stdout=open(args.r[-1], 'w'),
				stderr=open(args.temp_dir+tempfile_name+'bam2fq_stderr', 'w'))
		subprocess.check_call(['rm'] + bams)
		chrom = args.range[:args.range.find(':')]
		coord1 = args.range[args.range.find(':')+1:args.range.find('-')]
		coord2 = args.range[args.range.find('-')+1:]
		precollapse = args.temp_dir+tempfile_name+args.range+'.bed' # name of subsetted query file
		coordfile = open(args.temp_dir+tempfile_name+args.range+'.range.bed', 'wt') # write range to a bed file
		coordfile.write('\t'.join([chrom, coord1, coord2]))
		coordfile.close()
		if subprocess.call(['tabix', '-R', args.temp_dir+tempfile_name+args.range+'.range.bed', args.q],
			stdout=open(precollapse, 'w')):
			sys.stderr.write('Query file needs to be a sorted, bgzip-ed, tabix-indexed bed file if range is specified\n')
			return 1
	else:
		if args.q.endswith('psl'):
			sys.stderr.write('** Error. Flair no longer accepts PSL input. Please use psl_to_bed first.\n')
			return 1
		precollapse = args.q # query file unchanged
		args.r = args.r[0].split(',') if ',' in args.r[0] else args.r # read sequences
		for r in args.r:
			if not os.path.exists(args.q):
				sys.stderr.write('Check that read file {} exists\n'.format(r))
				return 1

	intermediate = []
	if args.p:
		if not args.quiet:
			sys.stderr.write('Filtering out reads without promoter-supported TSS\n')
		if subprocess.call([sys.executable, path+'pull_starts.py', args.q, args.temp_dir+tempfile_name+'tss.bed']):
			return 1
		if subprocess.call(['bedtools', 'intersect', '-a', args.temp_dir+tempfile_name+'tss.bed', '-b', args.p],
			stdout=open(args.temp_dir+tempfile_name+'promoter_intersect.bed', 'w')):
			return 1
		precollapse = args.o+'promoter_supported.bed' # filename of promoter-supported, corrected reads
		subprocess.check_call([sys.executable, path+'select_from_bed.py', args.temp_dir+tempfile_name+'promoter_intersect.bed',
			args.q, precollapse])
		intermediate += [args.temp_dir+tempfile_name+'tss.bed', precollapse]

	if args.threeprime:
		if not args.quiet:
			sys.stderr.write('Filtering out reads without TES support\n')
		if subprocess.call([sys.executable, path+'pull_starts.py', precollapse, args.temp_dir+tempfile_name+'tes.bed', 'reverse']):
			return 1
		if subprocess.call(['bedtools', 'intersect', '-a', args.temp_dir+tempfile_name+'tes.bed', '-b', args.threeprime],
			stdout=open(args.temp_dir+tempfile_name+'tes_intersect.bed', 'w')):
			return 1
		precollapse = args.o+'tes_supported.bed' # filename of 3' end-supported, corrected reads
		subprocess.check_call([sys.executable, path+'select_from_bed.py', args.temp_dir+tempfile_name+'tes_intersect.bed',
			args.q, precollapse])
		intermediate += [args.temp_dir+tempfile_name+'tes.bed', precollapse]

	if args.range: # for collapse_range, make sure the minimum number of reads is present
		count = len(open(precollapse).readlines())
		if count < min_reads:
			return 1

	if args.annotation_reliant:
		if not args.generate_map:
			args.generate_map = True
		if args.annotation_reliant == 'generate' or not args.annotated_bed:
			if not os.path.exists(args.f):
				if not args.f:
					sys.stderr.write('Please specify annotated gtf with -f for --annotation_reliant generate\n')
				else:
					sys.stderr.write('GTF file path does not exist\n')
				return 1

			if not args.quiet:
				sys.stderr.write('Making transcript fasta using annotated gtf and genome sequence\n')
			args.annotated_bed = args.o+'annotated_transcripts.bed'
			subprocess.check_call([sys.executable, path+'gtf_to_bed.py', args.f, args.annotated_bed, '--include_gene'])

		if args.annotation_reliant == 'generate':
			args.annotation_reliant = args.o+'annotated_transcripts.fa'
			subprocess.check_call([sys.executable, path+'bed_to_sequence.py', args.o+'annotated_transcripts.bed', args.g, args.annotation_reliant])

		if not args.quiet:
			sys.stderr.write('Aligning reads to reference transcripts\n')
		if subprocess.call(['minimap2', '-a', '-t', str(args.t), '-N', '4', args.annotation_reliant] + args.r,
			stdout=open(args.o+'annotated_transcripts.alignment.sam', 'w'),
			stderr=open(args.o+'annotated_transcripts.alignment.mm2_stderr', 'w')):
			sys.stderr.write('Minimap2 issue, check stderr file\n')
			return 1
		intermediate += [args.o+'annotated_transcripts.alignment.sam', args.o+'annotated_transcripts.alignment.mm2_stderr']

		if not args.quiet:
			sys.stderr.write('Counting supporting reads for annotated transcripts\n')
		count_cmd = [sys.executable, path+'count_sam_transcripts.py', '-s', args.o+'annotated_transcripts.alignment.sam',
			'-o', args.o+'annotated_transcripts.alignment.counts', '-t', str(args.t), '--quality', str(args.quality),
			'-w', str(args.w), '--generate_map', args.o+'annotated_transcripts.isoform.read.map.txt']
		if args.stringent:
			count_cmd += ['--stringent']
		if args.check_splice:
			count_cmd += ['--check_splice']
		if args.check_splice or args.stringent:
			count_cmd += ['-i', args.annotated_bed] # annotated isoform bed file
		if args.trust_ends:
			count_cmd += ['--trust_ends']
		if subprocess.call(count_cmd):
			sys.stderr.write('Failed at counting step for annotated isoform read support\n')
			return 1

		if not args.quiet:
			sys.stderr.write('Setting up unassigned reads for flair-collapse novel isoform detection\n')
		subprocess.check_call([sys.executable, path+'match_counts.py', args.o+'annotated_transcripts.alignment.counts',
			args.annotated_bed, str(min_reads), args.o+'annotated_transcripts.supported.bed'])

		subset_reads = args.o+'unassigned.fasta'
		subprocess.check_call([sys.executable, path+'subset_unassigned_reads.py', args.o+'annotated_transcripts.isoform.read.map.txt',
			precollapse, str(min_reads), args.o+'unassigned.bed']+args.r, stdout=open(subset_reads, 'w'))
		precollapse = args.o+'unassigned.bed'
		args.r = [subset_reads]
		intermediate += [subset_reads, precollapse]

	collapse_cmd = [sys.executable, path+'collapse_isoforms_precise.py', '-q', precollapse, '-t', str(args.t),
			'-m', str(args.max_ends), '-w', str(args.w), '-n', args.n, '-o', args.o+'firstpass.unfiltered.bed']
	if args.f and not args.no_end_adjustment:
		collapse_cmd += ['-f', args.f]
	if args.i:
		collapse_cmd += ['-i']
	if args.quiet:
		collapse_cmd += ['--quiet']

	if subprocess.call(collapse_cmd):
		return 1

	# filtering out subset isoforms with insuficient support
	filter_cmd = [sys.executable, path+'filter_collapsed_isoforms.py',
		args.o+'firstpass.unfiltered.bed', args.filter, args.o+'firstpass.bed', str(args.w)]
	if float(args.s) < 1:
		filter_cmd += ['keep']

	if subprocess.call(filter_cmd):
		return 1
	intermediate += [args.o+'firstpass.unfiltered.bed']

	# rename first-pass isoforms to annotated transcript IDs if they match
	if args.f:
		if not args.quiet:
			sys.stderr.write('Renaming isoforms using gtf\n')
		renaming_cmd = [sys.executable, path+'identify_gene_isoform.py',
			args.o+'firstpass.bed', args.f, args.o+'firstpass.named.bed']
		if args.annotation_reliant:
			renaming_cmd += ['--annotation_reliant']
		if subprocess.call(renaming_cmd):
			return 1
		subprocess.check_call(['mv', args.o+'firstpass.named.bed', args.o+'firstpass.bed'])
		if float(args.s) < 1:
			subprocess.check_call([sys.executable, path+'filter_isoforms_by_proportion_of_gene_expr.py',
				args.o+'firstpass.bed', str(args.s), args.o+'firstpass.bed'])

	if subprocess.call([sys.executable, path+'bed_to_sequence.py', args.o+'firstpass.bed',
		args.g, args.o+'firstpass.fa']):
		return 1

	# reassign reads to first-pass isoforms
	if not args.quiet:
		sys.stderr.write('Aligning reads to first-pass isoform reference\n')
	align_files = []
	alignout = args.temp_dir + tempfile_name + 'firstpass.'

	if args.mm2_args:
		args.mm2_args = args.mm2_args.split(',')
	# note that we cannot get a return value at this point because there's a pipe to samtools or count_sam_transcripts.py
	ps = subprocess.Popen(['minimap2', '-a', '-t', str(args.t), '-N', '4'] + args.mm2_args + [args.o+'firstpass.fa'] + args.r,
		stdout=subprocess.PIPE)

	# count the number of supporting reads for each first-pass isoform
	count_file = args.o+'firstpass.q.counts'
	count_cmd = [sys.executable, path+'count_sam_transcripts.py', '-s', '-',
		'-o', count_file, '-t', str(args.t), '--quality', str(args.quality), '-w', str(args.w)]
	if args.stringent:
		count_cmd += ['--stringent']
	if args.check_splice:
		count_cmd += ['--check_splice']
	if args.check_splice or args.stringent:
		count_cmd += ['-i', args.o+'firstpass.bed']
	if args.trust_ends:
		count_cmd += ['--trust_ends']
	if args.generate_map:
		count_cmd += ['--generate_map', args.o+'isoform.read.map.txt']
	if args.fusion_dist:
		count_cmd += ['--fusion_dist', str(args.fusion_dist)]
	if subprocess.call(count_cmd, stdin=ps.stdout):
		sys.stderr.write('Failed at counting step for isoform read support\n')
		return 1

	if not args.quiet:
		sys.stderr.write('Filtering isoforms by read coverage\n')
	match_count_cmd = [sys.executable, path+'match_counts.py', count_file,
		args.o+'firstpass.bed', str(min_reads), args.o+'isoforms.bed']
	if args.generate_map or args.annotation_reliant:
		match_count_cmd += ['--generate_map', args.o+'isoform.read.map.txt']
	subprocess.check_call(match_count_cmd)

	if args.annotation_reliant:
		subprocess.check_call([sys.executable, path+'filter_collapsed_isoforms_from_annotation.py', '-s', str(min_reads),
			'-i', args.o+'isoforms.bed', '--map_i', args.o+'isoform.read.map.txt',
			'-a', args.o+'annotated_transcripts.supported.bed', '--map_a', args.o+'annotated_transcripts.isoform.read.map.txt',
			'-o', args.o+'isoforms.bed', '--new_map', args.o+'combined.isoform.read.map.txt'])

	if not args.range: # also write .fa and .gtf files
		if args.longshot_bam:
			subprocess.check_call([sys.executable, path+'get_phase_sets.py', '-i', args.o+'isoforms.bed',
				'-b', args.longshot_bam, '-m', args.o+'combined.isoform.read.map.txt',
				'-o', args.o+'phase_sets.txt', '--out_iso', args.o+'isoforms.bed'])

		to_sequence_cmd = [sys.executable, path+'bed_to_sequence.py', args.o+'isoforms.bed',
			args.g, args.o+'isoforms.fa']
		if args.longshot_bam:
			to_sequence_cmd += ['--vcf', args.longshot_vcf, '--isoform_haplotypes', args.o+'phase_sets.txt', '--vcf_out', args.o+'flair.vcf']
		subprocess.check_call(to_sequence_cmd)
		if args.f:
			subprocess.check_call([sys.executable, path+'bed_to_gtf.py', args.o+'isoforms.bed'],
				stdout=open(args.o+'isoforms.gtf', 'w'))

	subprocess.check_call(['rm', '-rf', args.o+'firstpass.fa', alignout+'q.counts'])
	if not args.keep_intermediate:
		subprocess.check_call(['rm', args.o+'firstpass.q.counts', args.o+'firstpass.bed'])
		subprocess.check_call(['rm', '-rf'] + glob.glob(args.temp_dir+'*'+tempfile_name+'*') + align_files + intermediate)
	return args.o+'isoforms.bed', args.o+'isoforms.fa'


def quantify(isoform_sequences=''):
	parser = argparse.ArgumentParser(description='flair-quantify parse options',
		usage='flair quantify -r reads_manifest.tsv -i isoforms.fa [options]')
	parser.add_argument('quantify')
	required = parser.add_argument_group('required named arguments')
	if not isoform_sequences:
		required.add_argument('-r', '--reads_manifest', action='store', dest='r', type=str,
			required=True, help='Tab delimited file containing sample id, condition, batch, reads.fq')
		required.add_argument('-i', '--isoforms', action='store', dest='i',
			type=str, required=True, help='FastA of FLAIR collapsed isoforms')
	else:
		required.add_argument('--reads_manifest', action='store', dest='r', type=str,
			required=True, help='Tab delimited file containing sample id, condition, batch, reads.fq')
	parser.add_argument('-o', '--output', type=str, action='store', dest='o', default='flair.quantify',
		help='''output file name base for FLAIR quantify (default: flair.quantify)''')
	parser.add_argument('-t', '--threads', type=int,
		action='store', dest='t', default=4, help='minimap2 number of threads (4)')
	parser.add_argument('--temp_dir', default='', action='store', dest='temp_dir',
		help='''directory to put temporary files. use './" to indicate current directory
		(default: python tempfile directory)''')
	parser.add_argument('--sample_id_only', default=False, action='store_true', dest='sample_id_only',
		help='''only use sample id in output header''')
	parser.add_argument('--tpm', action='store_true', dest='tpm', default=False,
		help='Convert counts matrix to transcripts per million and output as a separate file named <output>.tpm.tsv')
	parser.add_argument('--quality', type=int, action='store', dest='quality', default=1,
		help='''minimum MAPQ of read assignment to an isoform (1)''')
	parser.add_argument('--trust_ends', default=False, action='store_true', dest='trust_ends',
		help='specify if reads are generated from a long read method with minimal fragmentation')
	parser.add_argument('--generate_map', default=False, action='store_true', dest='generate_map',
		help='''create read-to-isoform assignment files for each sample (default: not specified)''')
	parser.add_argument('--isoform_bed', '--isoformbed', default='', type=str, action='store', dest='isoforms',
		help='''isoform .bed file, must be specified if --stringent or check_splice is specified''')
	parser.add_argument('--stringent', default=False, action='store_true', dest='stringent',
		help='''Supporting reads must cover 80 percent of their isoform and extend at least 25 nt into the
                first and last exons. If those exons are themselves shorter than 25 nt, the requirement becomes
                'must start within 4 nt from the start" or "must end within 4 nt from the end" ''')
	parser.add_argument('--check_splice', default=False, action='store_true', dest='check_splice',
		help='''enforce coverage of 4 out of 6 bp around each splice site and no
		insertions greater than 3 bp at the splice site''')
	args, unknown = parser.parse_known_args()
	if unknown:
		sys.stderr.write('Quantify unrecognized arguments: {}\n'.format(' '.join(unknown)))
		if not isoform_sequences:
			return 1

	if isoform_sequences:
		args.i = isoform_sequences
		args.o += '.counts_matrix.tsv'
	if (args.stringent or args.check_splice):
		if not args.isoforms:
			sys.stderr.write('Please specify isoform models as .bed file using --isoform_bed\n')
			return 1
		elif not os.path.exists(args.isoforms):
			sys.stderr.write('Isoform models bed file path does not exist\n')
			return 1
		elif args.isoforms.endswith('.psl'):
			sys.stderr.write('** Error. Flair no longer accepts PSL input. Please use psl_to_bed first.\n')
			return 1
	if not os.path.exists(args.i):
		sys.stderr.write('Isoform sequences fasta file path does not exist\n')
		return 1

	try:
		import numpy as np
		import codecs
	except:
		sys.stderr.write('Numpy import error. Please pip install numpy. Exiting.\n')
		return 1

	samData = list()
	with codecs.open(args.r, 'r', encoding='utf-8', errors='ignore') as lines:
		for line in lines:
			cols = line.rstrip().split('\t')
			if len(cols) < 4:
				sys.stderr.write(f'Expected 4 columns in tab-delimited manifest.tsv, got len(cols). Exiting.\n')
				return 1

			sample, group, batch, readFile = cols
			if args.sample_id_only is False:
				if '_' in sample or '_' in group or '_' in batch:
					sys.stderr.write(f'Please do not use underscores in the id, condition, or batch fields of {args.r}. Exiting. \n')
					return 1

			readFileRoot = tempfile.NamedTemporaryFile().name
			if args.temp_dir != '':
				if not os.path.isdir(args.temp_dir):
					if subprocess.call(['mkdir', args.temp_dir]):
						sys.stderr.write('Could not make temporary directory {}\n'.format(args.temp_dir))
						return 1
				readFileRoot = args.temp_dir + '/' + readFileRoot[readFileRoot.rfind('/')+1:]

			if not os.path.exists(cols[3]):
				sys.stderr.write('Query file path does not exist: {}\n'.format(cols[3]))
				return 1
			samData.append(cols + [readFileRoot + '.sam'])
	sys.stderr.write('Writing temporary files with prefixes similar to {}\n'.format(readFileRoot))

	for num, sample in enumerate(samData, 0):
		sys.stderr.write('Step 1/3. Aligning sample %s_%s, %s/%s \r' % (sample[0], sample[2], num+1, len(samData)))
		mm2_command = ['minimap2', '-a', '-N', '4', '-t', str(args.t), args.i, sample[-2]]

		try:
			if subprocess.call(mm2_command, stdout=open(sample[-1], 'w'),
				stderr=open(sample[-1]+'.mm2_stderr.txt', 'w')):
				sys.stderr.write('Check {} file\n'.format(sample[-1]+'.mm2_stderr.txt'))
				return 1
		except:
			sys.stderr.write('''Possible minimap2 error, please check that all file, directory,
				and executable paths exist\n''')
			return 1
		subprocess.check_call(['rm', sample[-1]+'.mm2_stderr.txt'])
		sys.stderr.flush()

	countData = dict()
	for num, data in enumerate(samData):
		sample, group, batch, readFile, samOut = data
		sys.stderr.write('Step 2/3. Quantifying isoforms for sample %s_%s: %s/%s \r' % (sample, batch, num+1, len(samData)))

		count_cmd = [sys.executable, path+'count_sam_transcripts.py', '-s', samOut,
			'-o', samOut+'.counts.txt', '-t', str(args.t), '--quality', str(args.quality)]
		if args.trust_ends:
			count_cmd += ['--trust_ends']
		if args.stringent:
			count_cmd += ['--stringent']
		if args.check_splice:
			count_cmd += ['--check_splice']
		if args.check_splice or args.stringent:
			count_cmd += ['-i', args.isoforms]
		if args.generate_map:
			count_cmd += ['--generate_map', args.o+'.'+sample+'.'+group+'.isoform.read.map.txt']

		if subprocess.call(count_cmd):
			return 1
		for line in open(samOut+'.counts.txt'):
			line = line.rstrip().split('\t')
			iso, numreads = line[0], line[1]
			if iso not in countData:
				countData[iso] = np.zeros(len(samData))
			countData[iso][num] = numreads
		sys.stderr.flush()
		subprocess.check_call(['rm', samOut])

	sys.stderr.write('Step 3/3. Writing counts to {} \r'.format(args.o+'.counts.tsv'))
	countMatrix = open(args.o+'.counts.tsv', 'w')

	if args.sample_id_only:
		countMatrix.write('\t'.join(['ID']+[x[0] for x in samData])+'\n')
	else:
		countMatrix.write('ids\t%s\n' % '\t'.join(['_'.join(x[:3]) for x in samData]))

	features = sorted(list(countData.keys()))
	for f in features:
		countMatrix.write('%s\t%s\n' % (f, '\t'.join(str(x) for x in countData[f])))

	countMatrix.close()
	sys.stderr.flush()
	sys.stderr.write('\n')

	if args.tpm:
		subprocess.check_call([sys.executable, path+'counts_to_tpm.py', args.o+'.counts.tsv', args.o+'.tpm.tsv'])
	return args.o+'.counts.tsv'


def diffExp(counts_matrix=''):
	parser = argparse.ArgumentParser(description='flair-diffExp parse options',
		usage='flair diffExp -q counts_matrix.tsv --out_dir out_dir [options]')
	parser.add_argument('diffExp')
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

	scriptsBin = path
	runDE = scriptsBin + 'deFLAIR.py'
	DEcommand = [sys.executable, '-W ignore', runDE, '--filter', str(args.e), '--threads',
		str(args.t), '--outDir', args.o, '--matrix', args.q]
	if args.of:
		DEcommand += ['-of']

	if subprocess.call(DEcommand):
		return 1
	return


def emptyMatrix(infile):
	'''Returns true if file has only a header line'''
	with open(infile, 'r') as inf:
		if len(inf.readlines()) <= 1:
			return True
	return False


def diffSplice(isoforms='', counts_matrix=''):
	parser = argparse.ArgumentParser(description='flair-diffSplice parse options',
		usage='flair diffSplice -i isoforms.bed -q counts_matrix.tsv [options]')
	parser.add_argument('diffSplice')
	required = parser.add_argument_group('required named arguments')
	if not isoforms:
		required.add_argument('-i', '--isoforms', action='store', dest='i', required=True,
			type=str, help='isoforms in bed format')
		required.add_argument('-q', '--counts_matrix', action='store', dest='q',
			type=str, required=True, help='tab-delimited isoform count matrix from flair quantify module')
	required.add_argument('-o', '--out_dir', action='store', dest='o',
		type=str, required=True, help='Output directory for tables and plots.')
#	parser.add_argument('-o', '--output', action='store', dest='o', default='flair.diffsplice', type=str,
#		required=False, help='output file name base for FLAIR isoforms (default: flair.diffsplice)')
	parser.add_argument('-t', '--threads', action='store', dest='t',
		type=int, required=False, default=4, help='Number of threads for parallel DRIMSeq (4)')
	parser.add_argument('--test', action='store_true', dest='test',
		required=False, default=False, help='Run DRIMSeq statistical testing')
	parser.add_argument('--drim1', action='store', dest='drim1', type=int, required=False, default=6,
		help='''The minimum number of samples that have coverage over an AS event inclusion/exclusion
		for DRIMSeq testing; events with too few samples are filtered out and not tested (6)''')
	parser.add_argument('--drim2', action='store', dest='drim2', type=int, required=False, default=3,
		help='''The minimum number of samples expressing the inclusion of an AS event;
		 events with too few samples are filtered out and not tested (3)''')
	parser.add_argument('--drim3', action='store', dest='drim3', type=int, required=False, default=15,
		help='''The minimum number of reads covering an AS event inclusion/exclusion for DRIMSeq testing,
		 events with too few samples are filtered out and not tested (15)''')
	parser.add_argument('--drim4', action='store', dest='drim4', type=int, required=False, default=5,
		help='''The minimum number of reads covering an AS event inclusion for DRIMSeq testing,
		 events with too few samples are filtered out and not tested (5)''')
	parser.add_argument('--batch', action='store_true', dest='batch', required=False, default=False,
		help='''If specified with --test, DRIMSeq will perform batch correction''')
	parser.add_argument('--conditionA', action='store', dest='conditionA', required=False, default='',
		help='''Implies --test. Specify one condition corresponding to samples in the counts_matrix to be compared against
		condition2; by default, the first two unique conditions are used''')
	parser.add_argument('--conditionB', action='store', dest='conditionB', required=False, default='',
		help='''Specify another condition corresponding to samples in the counts_matrix to be compared against
		conditionA''')
	parser.add_argument('-of', '--out_dir_force', action='store_true', dest='of',
		required=False, help='''Specify this argument to force overwriting of files in
		an existing output directory''')
	args, unknown = parser.parse_known_args()
	if unknown:
		sys.stderr.write('DiffSplice unrecognized arguments: {}\n'.format(' '.join(unknown)))
		if not isoforms:
			return 1

	if isoforms:
		args.i = isoforms
		args.q = counts_matrix

	if not os.path.exists(args.q):
		sys.stderr.write('Counts matrix file path does not exist\n')
		return 1
	if not os.path.exists(args.i):
		sys.stderr.write('Isoform bed file path does not exist\n')
		return 1

	# Create output directory including a working directory for intermediate files.
	workdir = os.path.join(args.o, 'workdir')
	if args.of:
		if not os.path.exists(workdir):
			os.makedirs(workdir)
		pass
	elif not os.path.exists(args.o):
		try:
			os.makedirs(workdir, 0o700)
		except OSError as e:
			if e.errno != errno.EEXIST:
				raise
	else:
		sys.stderr.write(f'** Error. Name {args.o} already exists. Choose another name for out_dir\n')
		return 1
	if args.i.endswith('psl'):
		sys.stderr.write('** Error. Flair no longer accepts PSL input. Please use psl_to_bed first.\n')
		return 1

	filebase = os.path.join(args.o, 'diffsplice')
	subprocess.check_call([sys.executable, path+'call_diffsplice_events.py', args.i, filebase, args.q])
	subprocess.check_call([sys.executable, path+'es_as.py', args.i], stdout=open(filebase+'.es.events.tsv', 'w'))
	subprocess.check_call([sys.executable, path+'es_as_inc_excl_to_counts.py', args.q, filebase+'.es.events.tsv'],
		stdout=open(filebase+'.es.events.quant.tsv', 'w'))
	subprocess.check_call(['rm', filebase+'.es.events.tsv'])

	if args.test or args.conditionA:
		sys.stderr.write('DRIMSeq testing for each AS event type\n')
		drim1, drim2, drim3, drim4 = [str(x) for x in [args.drim1, args.drim2, args.drim3, args.drim4]]
		ds_command = [sys.executable, path+'runDS.py', '--threads', str(args.t), '--outDir', args.o,
			'--drim1', drim1, '--drim2', drim2, '--drim3', drim3, '--drim4', drim4]
		if args.batch:
			ds_command += ['--batch']
		if args.conditionA:
			if not args.conditionB:
				sys.stderr.write('Both conditionA and conditionB must be specified, or both left unspecified\n')
				return 1
			ds_command += ['--conditionA', args.conditionA, '--conditionB', args.conditionB]

		with open(workdir+'/ds.stderr.txt', 'w') as ds_stderr:
			for event in ['es', 'alt5', 'alt3', 'ir']:
				matrixfile = f'{filebase}.{event}.events.quant.tsv'
				if emptyMatrix(matrixfile):
					sys.stderr.write(f'{event} event matrix file empty, not running DRIMSeq\n')
					continue
				cur_command = ds_command + ['--matrix', matrixfile, '--prefix', event]
				if subprocess.call(cur_command, stderr=ds_stderr):
					print(f'\nDRIMSeq failed on {event} event')
					print(f'Check {workdir}/ds.stderr.txt for details.\nCommand was:')
					print(' '.join([x for x in cur_command]))
	else:
		# workdir only necessary for drimseq output
		os.rmdir(workdir)
	return


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
	if mode == 'align' or '1' in mode:
		status = align()
		if status == 1:
			sys.exit(1)
		else:
			aligned_reads = status

	if mode == 'correct' or '2' in mode:
		if aligned_reads:
			status = correct(aligned_reads=aligned_reads)
		else:
			status = correct()
		if status == 1:
			sys.exit(1)
		else:
			corrected_reads = status

	if mode == 'collapse' or ('3' in mode and '3.5' not in mode):
		if corrected_reads:
			status = collapse(corrected_reads=corrected_reads)
		else:
			status = collapse()
		if status == 1:
			sys.exit(1)
		else:
			isoforms, isoform_sequences = status

	if mode == 'collapse-range' or '3.5' in mode:
		tempfile_name = tempfile.NamedTemporaryFile().name
		run_id = tempfile_name[tempfile_name.rfind('/')+1:]

		if corrected_reads and not aligned_reads:
			sys.stderr.write('''Collapse 3.5 run consecutively without align module; will assume {}
			 to be the name of the aligned reads bam file\n'''.format(corrected_reads[:-18]+'.bam'))
			status = collapse_range(run_id=run_id,corrected_reads=corrected_reads,
				aligned_reads=corrected_reads[:-18]+'.bam')
		elif corrected_reads and aligned_reads:
			status = collapse_range(run_id=run_id,corrected_reads=corrected_reads, aligned_reads=aligned_reads)
		elif not corrected_reads and aligned_reads:
			sys.stderr.write('Correct module not run...\n')
			status = collapse_range(run_id=run_id,corrected_reads=aligned_reads, aligned_reads=aligned_reads)
		else:
			status = collapse_range(run_id=run_id,)
		if status == 1:
			sys.exit(1)
		else:
			isoforms, isoform_sequences = status
		mode = mode.replace('3.5', 'x')

	if mode == 'quantify' or '4' in mode:
		if isoform_sequences:
			status = quantify(isoform_sequences=isoform_sequences)
		else:
			status = quantify()
		if status == 1:
			sys.exit(1)
		else:
			counts_matrix = status

	if mode == 'diffexp' or '5' in mode:
		if counts_matrix:
			status = diffExp(counts_matrix=counts_matrix)
		else:
			status = diffExp()
		if status == 1:
			sys.exit(1)

	if mode == 'diffsplice' or '6' in mode:
		if counts_matrix and isoforms:
			status = diffSplice(isoforms=isoforms, counts_matrix=counts_matrix)
		elif not isoforms and counts_matrix:
			sys.stderr.write('DiffSplice run consecutively without collapse module, exiting\n')
			sys.exit(1)
		else:
			status = diffSplice()
		if status == 1:
			sys.exit(1)

	if mode == '--version':
		sys.stderr.write('FLAIR v1.7.0\n')


if __name__ == '__main__':
    main()
