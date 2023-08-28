#! /usr/bin/env python3

import sys
import argparse
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import tempfile
import glob
import pipettor
from multiprocessing import Pool
# TODO: put all of these in utils.py
from pull_starts import pull_starts
from select_from_bed import select_from_bed
from bed_to_sequence import bed_to_sequence
from gtf_to_bed import gtf_to_bed
from match_counts import match_counts
from filter_collapsed_isoforms import filter_collapsed_isoforms
from identify_gene_isoform import identify_gene_isoform
from filter_collapsed_isoforms_from_annotation import filter_collapsed_isoforms_from_annotation
from get_phase_sets import get_phase_sets
from bed_to_gtf import bed_to_gtf
from subset_unassigned_reads import subset_unassigned_reads
from filter_isoforms_by_proportion_of_gene_expr import filter_isoforms_by_proportion_of_gene_expr

run_id = 'removeme'
# TODO: do not redefine args variables, it breaks your head.

def collapse(genomic_range='', corrected_reads=''):
	parser = argparse.ArgumentParser(description='flair-collapse parse options',
		usage='''flair collapse -g genome.fa -q query.bed
		-r <reads.fq>/<reads.fa> [options]''')
	required = parser.add_argument_group('required named arguments')
	if not corrected_reads:
		required.add_argument('-q', '--query', type=str, default='', required=True,
			  help='bed file of aligned/corrected reads')
	required.add_argument('-g', '--genome',  
		type=str, required=True, help='FastA of reference genome')
	required.add_argument('-r', '--reads',   nargs='+',
		type=str, required=True, help='FastA/FastQ files of raw reads, can specify multiple files')
	parser.add_argument('-o', '--output', default='flair.collapse',
		  help='output file name base for FLAIR isoforms (default: flair.collapse)')
	parser.add_argument('-t', '--threads', type=int, default=4, 
		help='minimap2 number of threads (4)')
	parser.add_argument('-f', '--gtf', default='',  
		help='GTF annotation file, used for renaming FLAIR isoforms to annotated isoforms and adjusting TSS/TESs')
	parser.add_argument('--generate_map', default=False, action='store_true', 
		help='''specify this argument to generate a txt file of read-isoform assignments
		(default: not specified)''')
	parser.add_argument('--annotation_reliant', default=False,  
		help='''specify transcript fasta that corresponds to transcripts in the gtf to run annotation-
		reliant flair collapse; to ask flair to make transcript sequences given the gtf and genome fa,
		type --annotation_reliant generate''')
	# supporting read assignment options
	parser.add_argument('-s', '--support', type=float, default=3.0,  
		help='''minimum number of supporting reads for an isoform;
		if s < 1, it will be treated as a percentage of expression of the gene (3)''')
	parser.add_argument('--stringent', default=False, action='store_true', 
		help='''specify if all supporting reads need to be full-length
		(80%% coverage and spanning 25 bp of the first and last exons)''')
	parser.add_argument('--check_splice', default=False, action='store_true', 
		help='''enforce coverage of 4 out of 6 bp around each splice site and no
		insertions greater than 3 bp at the splice site''')
	parser.add_argument('--trust_ends', default=False, action='store_true', 
		help='specify if reads are generated from a long read method with minimal fragmentation')
	parser.add_argument('--quality', type=int,  default=1,
		help='minimum MAPQ of read assignment to an isoform (1)')
	# variant options
	parser.add_argument('--longshot_bam',  default='',
		help='bam from longshot containing haplotype information for each read')
	parser.add_argument('--longshot_vcf',  default='',
		help='vcf from longshot')
	# ends options
	parser.add_argument('-w', '--end_window', type=int, default=100,  
		help='window size for comparing TSS/TES (100)')
	parser.add_argument('-p', '--promoters',   default='',
		help='promoter regions bed file to identify full-length reads')
	parser.add_argument('--3prime_regions',  dest='threeprime', default='',
		help='TES regions bed file to identify full-length reads')
	parser.add_argument('-n', '--no_redundant', default='none',  
		help='''For each unique splice junction chain, report options include:
		none--best TSSs/TESs chosen for each unique set of splice junctions;
		longest--single TSS/TES chosen to maximize length;
		best_only--single most supported TSS/TES used in conjunction chosen (none)''')
	parser.add_argument('-i', '--isoformtss', default=False, action='store_true', 
		help='''when specified, TSS/TES for each isoform will be determined from supporting reads
		for individual isoforms (default: not specified, determined at the gene level)''')
	parser.add_argument('--no_gtf_end_adjustment', default=False, action='store_true',
		dest='no_end_adjustment',
		help='''when specified, TSS/TES from the gtf provided with -f will not be used to adjust
		isoform TSSs/TESs each isoform will be determined from supporting reads''')
	parser.add_argument('--max_ends', type=int, default=2,  
		help='maximum number of TSS/TES picked per isoform (2)')
	parser.add_argument('--filter', default='default',  
		help='''Report options include:
		nosubset--any isoforms that are a proper set of another isoform are removed;
		default--subset isoforms are removed based on support;
		comprehensive--default set + all subset isoforms;
		ginormous--comprehensive set + single exon subset isoforms''')
	# other
	parser.add_argument('--temp_dir', default='',  
		help='directory for temporary files. use "./" to indicate current directory (default: python tempfile directory)')
	parser.add_argument('--keep_intermediate', default=False, action='store_true', 
		help='''specify if intermediate and temporary files are to be kept for debugging.
		Intermediate files include: promoter-supported reads file,
		read assignments to firstpass isoforms''')
	parser.add_argument('--mm2_args', type=str, default=[], 
		help='''additional minimap2 arguments when aligning reads first-pass transcripts;
		separate args by commas, e.g. --mm2_args=-I8g,--MD ''')
	parser.add_argument('--quiet', default=False, action='store_true', 
			help='''Suppress progress statements from being printed''')
	parser.add_argument('--annotated_bed', default=False,  
		help='''annotation_reliant also requires a bedfile of annotated isoforms; if this isn't provided,
		flair collapse will generate the bedfile from the gtf. eventually this argument will be removed''')
	parser.add_argument('--range', default='',  
		help='''interval for which to collapse isoforms, formatted chromosome:coord1-coord2 or
		tab-delimited; if a range is specified, then the aligned reads bam must be specified with -r
		and the query must be a sorted, bgzip-ed bed file''')
	
	no_arguments_passed = len(sys.argv) == 1
	if no_arguments_passed:
		parser.print_help()
		sys.exit(1)
#	if 'collapse' in sys.argv:
#		sys.argv.remove('collapse')
	args, unknown = parser.parse_known_args()
	#if unknown and not args.quiet:
	if unknown:
		sys.stderr.write('Collapse unrecognized arguments: {}\n'.format(' '.join(unknown)))

	if corrected_reads:
		query = corrected_reads
	else:
		query = args.query

	# housekeeping stuff
	tempfile_dir = tempfile.NamedTemporaryFile().name
	tempfile_name = tempfile_dir[tempfile_dir.rfind('/')+1:]+'.'
	if args.temp_dir == '':
		args.temp_dir = tempfile_dir+'/'
		if not args.quiet:
			sys.stderr.write('Writing temporary files to {}\t\n'.format(args.temp_dir))
	if not os.path.isdir(args.temp_dir): # make temporary directory
		pipettor.run([('mkdir', args.temp_dir)])

	if args.temp_dir[-1] != '/':
		args.temp_dir += '/'

	if genomic_range: # this module was called internally from collapse_range
		args.range = genomic_range
		args.output = args.temp_dir+run_id
		query = args.temp_dir+run_id+'.sorted.bed.gz'
		args.quiet = True

	args.quality = '0' if args.trust_ends else args.quality
	args.output += '.'
	min_reads = float(args.support) if float(args.support) >= 1 else 3

	if ',' in args.reads[0]:
		args.reads = args.reads[0].split(',')
	for rfile in args.reads:
		if not os.path.exists(rfile):
			sys.stderr.write(f'Read file path does not exist: {rfile}\n')
			sys.exit(1)
	if not os.path.exists(query):
		sys.stderr.write('Query file path does not exist: {}\n'.format(query))
		sys.exit(1)
	if not os.path.exists(args.genome):
		sys.stderr.write('Genome file path does not exist: {}\n'.format(args.genome))
		sys.exit(1)
	if os.stat(query).st_size == 0:
		sys.stderr.write('Query file is empty\n')
		sys.exit(1)
	if float(args.support) < 1 and not args.gtf:
		sys.stderr.write('Provide gtf for gene grouping if -s is percentage of total gene expression\n')
		sys.exit(1)

	if args.range:
		# subset out the read sequences and corrected reads corresponding to the specified range
		if '\t' in args.range:
			args.range = args.range.split('\t')
			args.range = args.range[0]+':'+args.range[1]+'-'+args.range[2]
		args.output += args.range+'.'
		if args.reads[0][-3:] != 'bam':
			sys.stderr.write('Must provide genome alignment BAM with -r if range is specified\n')
			sys.exit(1)
		bams = []
		for i in range(len(args.reads)): # subset bam file for alignments within range
			bams += [args.temp_dir+tempfile_name+args.range+str(i)+'.bam']
			samtools_cmd = ('samtools', 'view', '-h', args.reads[i], args.range)
			pipettor.Popen([samtools_cmd], 'w', stdout=bams[-1])

		args.reads = []
		for i in range(len(bams)): # read sequences of the alignments within range
			args.reads += [bams[i][:-3]+'fasta']
			pipettor.Popen([('samtools', 'fasta', bams[i])], 'w', stdout=args.reads[-1]) # TODO add stderr
#			subprocess.check_call(['samtools', 'fasta', bams[i]],
#				stdout=open(args.reads[-1], 'w'),
#				stderr=open(args.temp_dir+tempfile_name+'bam2fq_stderr', 'w'))
		pipettor.run([('rm', bams)])  # TODO: does this work, or need join?
		chrom = args.range[:args.range.find(':')]
		coord1 = args.range[args.range.find(':')+1:args.range.find('-')]
		coord2 = args.range[args.range.find('-')+1:]
		precollapse = args.temp_dir+tempfile_name+args.range+'.bed' # name of subsetted query file
		coordfile = open(args.temp_dir+tempfile_name+args.range+'.range.bed', 'wt') # write range to a bed file
		coordfile.write('\t'.join([chrom, coord1, coord2]))
		coordfile.close()
		try:
			tabix_cmd = ('tabix', '-R', args.temp_dir+tempfile_name+args.range+'.range.bed', query)
			pipettor.Popen([tabix_cmd], 'w', stdout=precollapse)
#		if subprocess.call(['tabix', '-R', args.temp_dir+tempfile_name+args.range+'.range.bed', query],
#			stdout=open(precollapse, 'w')):
#			
		except Exception as e:
			sys.stderr.write('Query file needs to be a sorted, bgzip-ed, tabix-indexed bed file if range is specified\n')
			sys.exit(1)
	else:
		if query.endswith('psl'):
			sys.stderr.write('** Error. Flair no longer accepts PSL input. Please use psl_to_bed first.\n')
			sys.exit(1)
		precollapse = query # query file unchanged
		args.reads = args.reads[0].split(',') if ',' in args.reads[0] else args.reads # read sequences
		for r in args.reads:
			if not os.path.exists(query):
				sys.stderr.write('Check that read file {} exists\n'.format(r))
				sys.exit(1)

	intermediate = []
	if args.promoters:
		if not args.quiet:
			sys.stderr.write('Filtering out reads without promoter-supported TSS\n')
		tss_bedfile = args.temp_dir+tempfile_name+'tss.bed'
		pull_starts(query, tss_bedfile)
		bedtools_output = args.temp_dir+tempfile_name+'promoter_intersect.bed'
		bedtools_cmd = ('bedtools', 'intersect', '-a', tss_bedfile, '-b', args.promoters)
		pipettor.Popen([bedtools_cmd], 'w', stdout=bedtools_output)		
		# if subprocess.call(['bedtools', 'intersect', '-a', tss_bedfile, '-b', args.promoters],
		# 	stdout=open(args.temp_dir+tempfile_name+'promoter_intersect.bed', 'w')):
		# 	sys.exit(1)
		precollapse = args.output+'promoter_supported.bed' # filename of promoter-supported, corrected reads
		select_from_bed(bedtools_output, query, precollapse)
#		subprocess.check_call([sys.executable, path+'select_from_bed.py', args.temp_dir+tempfile_name+'promoter_intersect.bed',
#			args.query, precollapse])
		intermediate += [tss_bedfile, precollapse]

	if args.threeprime:
		if not args.quiet:
			sys.stderr.write('Filtering out reads without TES support\n')
		tes_bedfile = args.temp_dir+tempfile_name+'tes.bed'
		pull_starts(precollapse, tes_bedfile, reverse=True)
#		if subprocess.call([sys.executable, path+'pull_starts.py', precollapse, tes_bedfile, 'reverse']):
#			sys.exit(1)
		bedtools_output = args.temp_dir+tempfile_name+'tes_intersect.bed'
		bedtools_cmd = ('bedtools', 'intersect', '-a', tes_bedfile, '-b', args.threeprime)
		pipettor.Popen([bedtools_cmd], 'w', stdout=bedtools_output)	
#		if subprocess.call(['bedtools', 'intersect', '-a', tes_bedfile, '-b', args.threeprime],
#			stdout=open(args.temp_dir+tempfile_name+'tes_intersect.bed', 'w')):
#			sys.exit(1)
		precollapse = args.output+'tes_supported.bed' # filename of 3' end-supported, corrected reads
		select_from_bed(bedtools_output, query, precollapse)
#		subprocess.check_call([sys.executable, path+'select_from_bed.py', args.temp_dir+tempfile_name+'tes_intersect.bed',
#			args.query, precollapse])
		intermediate += [tes_bedfile, precollapse]

	if args.range: # for collapse_range, make sure the minimum number of reads is present
		count = len(open(precollapse).readlines())
		if count < min_reads:
			sys.stderr.write('ERROR: not enough reads in', precollapse)
			sys.exit(1)

	if args.annotation_reliant:
		if not args.generate_map:
			args.generate_map = True
		if args.annotation_reliant == 'generate' or not args.annotated_bed:
			if not os.path.exists(args.gtf):
				if not args.gtf:
					sys.stderr.write('Please specify annotated gtf with -f for --annotation_reliant generate\n')
				else:
					sys.stderr.write('GTF file path does not exist\n')
				sys.exit(1)

			if not args.quiet:
				sys.stderr.write('Making transcript fasta using annotated gtf and genome sequence\n')
			args.annotated_bed = args.output+'annotated_transcripts.bed'
			gtf_to_bed(args.annotated_bed, args.gtf, include_gene=True, chrom_sizes=False)
#			subprocess.check_call([sys.executable, path+'gtf_to_bed.py', args.gtf, args.annotated_bed, '--include_gene'])

		if args.annotation_reliant == 'generate':
			args.annotation_reliant = args.output+'annotated_transcripts.fa'
			bed_to_sequence(query=args.output+'annotated_transcripts.bed', genome=args.genome, 
				outfilename=args.annotation_reliant)

#			subprocess.check_call([sys.executable, path+'bed_to_sequence.py', 
#args.output+'annotated_transcripts.bed', args.genome, args.annotation_reliant])

		# minimap (results are piped into count_sam_transcripts.py)
		mm2_cmd = ['minimap2', '-a', '-t', str(args.threads), '-N', '4', args.annotation_reliant] + args.reads
		mm2_cmd = tuple(mm2_cmd)
		
#		if subprocess.call(['minimap2', '-a', '-t', str(args.threads), '-N', '4', args.annotation_reliant] + args.reads,
#			stdout=open(args.output+'annotated_transcripts.alignment.sam', 'w'),
#			stderr=open(args.output+'annotated_transcripts.alignment.mm2_stderr', 'w')):
#			sys.stderr.write('Minimap2 issue, check stderr file\n')
#			sys.exit(1)
#		intermediate += [args.output+'annotated_transcripts.alignment.sam', args.output+'annotated_transcripts.alignment.mm2_stderr']
		# count sam transcripts ; the dash at the end means STDIN
		count_cmd = ['count_sam_transcripts.py', '--sam', '-',
			'-o', args.output+'annotated_transcripts.alignment.counts', '-t', str(args.threads), 
			'--quality', str(args.quality), '-w', str(args.end_window), 
			'--generate_map', args.output+'annotated_transcripts.isoform.read.map.txt']
		if args.stringent:
			count_cmd += ['--stringent']
		if args.check_splice:
			count_cmd += ['--check_splice']
		if args.check_splice or args.stringent:
			count_cmd += ['-i', args.annotated_bed] # annotated isoform bed file
		if args.trust_ends:
			count_cmd += ['--trust_ends']
#		if subprocess.call(count_cmd):
#			sys.stderr.write('Failed at counting step for annotated isoform read support\n')
#			sys.exit(1)
		count_cmd = tuple(count_cmd)

		if not args.quiet:
			sys.stderr.write('Aligning reads to reference transcripts\n')
			sys.stderr.write('Counting supporting reads for annotated transcripts\n')
		pipettor.run([mm2_cmd, count_cmd]) 

		if not args.quiet:
			sys.stderr.write('Setting up unassigned reads for flair-collapse novel isoform detection\n')
		counts_file = args.output+'annotated_transcripts.alignment.counts'
		supported_bed = args.output+'annotated_transcripts.supported.bed'
		match_counts(counts_file=counts_file, output_file=supported_bed, psl=args.annotated_bed, min_reads=args.support)
# NOTE: this won't work as is anymore because match_counts.py now requires named input args
#		subprocess.check_call([sys.executable, path+'match_counts.py', 
# args.output+'annotated_transcripts.alignment.counts',
#			args.annotated_bed, str(min_reads), args.output+'annotated_transcripts.supported.bed'])

		subset_reads = args.output+'unassigned.fasta'
		# NOTE: min_reads is passed but not used
		subset_unassigned_reads(readmap=args.output+'annotated_transcripts.isoform.read.map.txt',
                        query=precollapse, support=str(min_reads), output=args.output+'unassigned.bed', 
			fastx=args.reads, outfa=subset_reads)
#		subprocess.check_call([sys.executable, path+'subset_unassigned_reads.py', args.output+'annotated_transcripts.isoform.read.map.txt',
#			precollapse, str(min_reads), args.output+'unassigned.bed']+args.reads, stdout=open(subset_reads, 'w'))
		# TODO: Get rid of this args renaming!
		precollapse = args.output+'unassigned.bed'
		args.reads = [subset_reads]
		intermediate += [subset_reads, precollapse]


	gtfname = None
	if args.gtf and not args.no_end_adjustment:
		gtfname=args.gtf

	# TODO: This program uses pool and map, which makes it difficult to capture in a function
#	collapse_isoforms_precise(queryfile=precollapse, threads=args.threads, 
#			   max_results=args.max_ends, window=args.end_window, no_redundant=args.no_redundant, 
#			   outputfname=args.output+'firstpass.unfiltered.bed', gtfname=gtfname, isoformtss=args.isoformtss, 
#			   quiet=args.quiet)
	collapse_cmd = ['collapse_isoforms_precise.py', '-q', precollapse, '-t', str(args.threads),
		'-m', str(args.max_ends), '-w', str(args.end_window), '-n', args.no_redundant, 
		'-o', args.output+'firstpass.unfiltered.bed']
	if args.gtf and not args.no_end_adjustment:
		collapse_cmd += ['-f', args.gtf]
	if args.isoformtss:
		collapse_cmd += ['-i']
	if args.quiet:
		collapse_cmd += ['--quiet']
	collapse_cmd = tuple(collapse_cmd)
	pipettor.run([collapse_cmd])


	# filtering out subset isoforms with insufficient support
	keep_extra_column = False
	if float(args.support) < 1:
		keep_extra_column = True
	filter_collapsed_isoforms(queryfile=args.output+'firstpass.unfiltered.bed', mode=args.filter, 
		tol=args.end_window, outfile=args.output+'firstpass.bed', keep_extra_column=keep_extra_column)
#	filter_cmd = [sys.executable, path+'filter_collapsed_isoforms.py',
#		args.output+'firstpass.unfiltered.bed', args.filter, args.output+'firstpass.bed', str(args.end_window)]
#	if float(args.support) < 1:
#		filter_cmd += ['keep']
#	if subprocess.call(filter_cmd):
#		sys.exit(1)
	intermediate += [args.output+'firstpass.unfiltered.bed']

	# rename first-pass isoforms to annotated transcript IDs if they match
	if args.gtf:
		if not args.quiet:
			sys.stderr.write('Renaming isoforms using gtf\n')
		identify_gene_isoform(query=args.output+'firstpass.bed', gtf=args.gtf, 
			outfilename=args.output+'firstpass.named.bed', annotation_reliant=args.annotation_reliant)
#		renaming_cmd = [sys.executable, path+'identify_gene_isoform.py',
#			args.output+'firstpass.bed', args.gtf, args.output+'firstpass.named.bed']
#		if args.annotation_reliant:
#			renaming_cmd += ['--annotation_reliant']
#		if subprocess.call(renaming_cmd):
#			sys.exit(1)

		# TODO: get rid of this renaming (maybe by copying, or otherwise by having a currentbed variable)
		os.rename(args.output+'firstpass.named.bed', args.output+'firstpass.bed')

		if float(args.support) < 1:
			filter_isoforms_by_proportion_of_gene_expr(isoforms=args.output+'firstpass.bed', 
				outfilename=args.output+'firstpass.filtered.bed', support=args.support)
			os.rename(args.output+'firstpass.filtered.bed', args.output+'firstpass.bed')
#			subprocess.check_call(['mv', args.output+'firstpass.filtered.bed', args.output+'firstpass.bed'])

# please note: in the original, input and output filenames were identical. Don't do this
#			subprocess.check_call([sys.executable, path+'filter_isoforms_by_proportion_of_gene_expr.py',
#				args.output+'firstpass.bed', str(args.support), args.output+'firstpass.bed'])

	bed_to_sequence(query=args.output+'firstpass.bed', genome=args.genome, outfilename=args.output+'firstpass.fa')
	# if subprocess.call([sys.executable, path+'bed_to_sequence.py', args.output+'firstpass.bed',
	# 	args.genome, args.output+'firstpass.fa']):
	# 	sys.exit(1)

	# reassign reads to first-pass isoforms
	if not args.quiet:
		sys.stderr.write('Aligning reads to first-pass isoform reference\n')
#	align_files = []  # TODO: this doesn't get filled
	alignout = args.temp_dir + tempfile_name + 'firstpass.'

	# minimap
	if args.mm2_args:
		args.mm2_args = args.mm2_args.split(',')
	mm2_cmd = ['minimap2', '-a', '-t', str(args.threads), '-N', '4'] + args.mm2_args + [args.output+'firstpass.fa'] + args.reads
#	ps = subprocess.Popen(['minimap2', '-a', '-t', str(args.threads), '-N', '4'] + args.mm2_args + [args.output+'firstpass.fa'] + args.reads,
#		stdout=subprocess.PIPE)

	# count the number of supporting reads for each first-pass isoform
	count_file = args.output+'firstpass.q.counts'
	count_cmd = ['count_sam_transcripts.py', '--sam', '-',
		'-o', count_file, '-t', str(args.threads), '--quality', str(args.quality), '-w', str(args.end_window)]
	if args.stringent:
		count_cmd += ['--stringent']
	if args.check_splice:
		count_cmd += ['--check_splice']
	if args.check_splice or args.stringent:
		count_cmd += ['-i', args.output+'firstpass.bed']
	if args.trust_ends:
		count_cmd += ['--trust_ends']
	if args.generate_map:
		count_cmd += ['--generate_map', args.output+'isoform.read.map.txt']
	count_cmd = tuple(count_cmd)

	if not args.quiet:
		sys.stderr.write('Aligning reads to firstpass transcripts\n')
		sys.stderr.write('Counting supporting reads for firstpass transcripts\n')
	pipettor.run([mm2_cmd, count_cmd]) 
#	if subprocess.call(count_cmd, stdin=ps.stdout):
#		sys.stderr.write('Failed at counting step for isoform read support\n')
#		sys.exit(1)

	if not args.quiet:
		sys.stderr.write('Filtering isoforms by read coverage\n')
	
	mc_output = args.output+'isoforms.bed'
	firstpassbed = args.output+'firstpass.bed'
	isoform_file = False
#	match_count_cmd = [sys.executable, path+'match_counts.py', count_file,
#		args.output+'firstpass.bed', str(min_reads), args.output+'isoforms.bed']
	# TODO: make a test that shows this has different output
	if args.generate_map or args.annotation_reliant:
		isoform_file = args.output+'isoform.read.map.txt'
#		match_count_cmd += ['--generate_map', args.output+'isoform.read.map.txt']
#	subprocess.check_call(match_count_cmd)
	match_counts(counts_file=count_file, output_file=mc_output, psl=firstpassbed, min_reads=min_reads, 
	      isoform_file = isoform_file)

	if args.annotation_reliant:
		filter_output = args.output+'isoforms.filter.bed'
		filter_collapsed_isoforms_from_annotation(support=min_reads, 
			queryfile=mc_output,
			map_i=args.output+'isoform.read.map.txt', 
			annotation=args.output+'annotated_transcripts.supported.bed', 
			map_a=args.output+'annotated_transcripts.isoform.read.map.txt',
			outputfile=filter_output,
			new_map=args.output+'combined.isoform.read.map.txt')
		os.rename(filter_output, mc_output)
#		subprocess.check_call([sys.executable, path+'filter_collapsed_isoforms_from_annotation.py', '-s', str(min_reads),
#			'-i', args.output+'isoforms.bed', '--map_i', args.output+'isoform.read.map.txt',
#			'-a', args.output+'annotated_transcripts.supported.bed', '--map_a', args.output+'annotated_transcripts.isoform.read.map.txt',
#			'-o', args.output+'isoforms.bed', '--new_map', args.output+'combined.isoform.read.map.txt'])

	# TODO: Test this section, longshot has not been tested at all
	if not args.range: # also write .fa and .gtf files
		if args.longshot_bam:
			get_phase_sets(isoforms=args.output+'isoforms.bed', 
				bam=args.longshot_bam, 
				isoform_reads_map=args.output+'combined.isoform.read.map.txt', 
				output=args.output+'phase_sets.txt', 
				outiso=args.output+'isoforms.ls.bed')
			os.rename(args.output+'isoforms.ls.bed', args.output+'isoforms.bed')

#			subprocess.check_call([sys.executable, path+'get_phase_sets.py', '-i', args.output+'isoforms.bed',
#				'-b', args.longshot_bam, '-m', args.output+'combined.isoform.read.map.txt',
#				'-o', args.output+'phase_sets.txt', '--out_iso', args.output+'isoforms.bed'])

			bed_to_sequence(query=args.output+'isoforms.bed', genome=args.genome, 
				outfilename=args.output+'isoforms.fa', 
				vcfinput=args.longshot_vcf, isoform_haplotypes=args.output+'phase_sets.txt', 
				vcf_out=args.output+'flair.vcf')
		else:
			bed_to_sequence(query=args.output+'isoforms.bed', genome=args.genome, 
				outfilename=args.output+'isoforms.fa')
			
		# to_sequence_cmd = [sys.executable, path+'bed_to_sequence.py', args.output+'isoforms.bed',
		# 	args.genome, args.output+'isoforms.fa']
		# if args.longshot_bam:
		# 	to_sequence_cmd += ['--vcf', args.longshot_vcf, '--isoform_haplotypes', args.output+'phase_sets.txt', '--vcf_out', args.output+'flair.vcf']
		# subprocess.check_call(to_sequence_cmd)
		if args.gtf:
			bed_to_gtf(query=args.output+'isoforms.bed', outputfile=args.output+'isoforms.gtf') 
		# 	subprocess.check_call([sys.executable, path+'bed_to_gtf.py', args.output+'isoforms.bed'],
		# 		stdout=open(args.output+'isoforms.gtf', 'w'))

	files_to_remove = [args.output+'firstpass.fa', alignout+'q.counts']
	#subprocess.check_call(['rm', '-rf', args.output+'firstpass.fa', alignout+'q.counts'])
	if not args.keep_intermediate:
		files_to_remove += [args.output+'firstpass.q.counts', args.output+'firstpass.bed'] + intermediate
		files_to_remove += glob.glob(args.temp_dir+'*'+tempfile_name+'*') # TODO: CHECK
		#subprocess.check_call(['rm', args.output+'firstpass.q.counts', args.output+'firstpass.bed'])
		#subprocess.check_call(['rm', '-rf'] + glob.glob(args.temp_dir+'*'+tempfile_name+'*') + align_files + intermediate)
	return [args.output+'isoforms.bed', args.output+'isoforms.fa']

if __name__ == "__main__":
    collapse()
