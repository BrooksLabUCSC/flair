#! /usr/bin/env python3

import sys
import argparse
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import tempfile
import glob
import pipettor
# TODO: put all of these in utils.py
from flair.pull_starts import pull_starts
from flair.select_from_bed import select_from_bed
from flair.bed_to_sequence import bed_to_sequence
from flair.gtf_to_bed import gtf_to_bed
from flair.match_counts import match_counts
from flair.filter_collapsed_isoforms import filter_collapsed_isoforms
from flair.identify_gene_isoform import identify_gene_isoform
from flair.filter_collapsed_isoforms_from_annotation import filter_collapsed_isoforms_from_annotation
from flair.get_phase_sets import get_phase_sets
from flair.bed_to_gtf import bed_to_gtf
from flair.subset_unassigned_reads import subset_unassigned_reads
from flair.filter_isoforms_by_proportion_of_gene_expr import filter_isoforms_by_proportion_of_gene_expr

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
	parser.add_argument('--remove_internal_priming', default=False, action='store_true',
						help='specify if want to remove reads with internal priming')
	parser.add_argument('--intprimingthreshold', type=int, default=12,
						help='number of bases that are at leas 75%% As required to call read as internal priming')
	parser.add_argument('--intprimingfracAs', type=float, default=0.6,
						help='number of bases that are at leas 75%% As required to call read as internal priming')

	no_arguments_passed = len(sys.argv) == 1
	if no_arguments_passed:
		parser.print_help()
		sys.exit(1)

	args = parser.parse_args()
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
		pipettor.run(['rm'] + bams)

		chrom = args.range[:args.range.find(':')]
		coord1 = args.range[args.range.find(':')+1:args.range.find('-')]
		coord2 = args.range[args.range.find('-')+1:]
		precollapse = args.temp_dir+tempfile_name+args.range+'.bed' # name of subsetted query file
		coordfile = open(args.temp_dir+tempfile_name+args.range+'.range.bed', 'wt') # write range to a bed file
		coordfile.write('\t'.join([chrom, coord1, coord2]))
		coordfile.close()

		# input bedfile must be tabix compressed
		try:
			tabix_cmd = ('tabix', '-R', args.temp_dir+tempfile_name+args.range+'.range.bed', query)
			pipettor.Popen([tabix_cmd], 'w', stdout=precollapse)
		except Exception as ex:
			raise Exception("** tabix FAILED for %s. File needs to be a sorted, bgzip-ed, tabix-indexed bed file if range is specified") from ex
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
		# get the starts
		tss_bedfile = args.temp_dir+tempfile_name+'tss.bed'
		pull_starts(query, tss_bedfile)

		# intersect with the promoters
		bedtools_output = args.temp_dir+tempfile_name+'promoter_intersect.bed'
		bedtools_cmd = ('bedtools', 'intersect', '-a', tss_bedfile, '-b', args.promoters)
		pipettor.Popen([bedtools_cmd], 'w', stdout=bedtools_output)

		# select reads that contain promoters
		precollapse = args.output+'promoter_supported.bed' # filename of promoter-supported, corrected reads
		select_from_bed(bedtools_output, query, precollapse)
		intermediate += [tss_bedfile, precollapse]

	if args.threeprime:
		if not args.quiet:
			sys.stderr.write('Filtering out reads without TES support\n')

		# get the ends
		tes_bedfile = args.temp_dir+tempfile_name+'tes.bed'
		pull_starts(precollapse, tes_bedfile, reverse=True)

		# intersect with known gene ends
		bedtools_output = args.temp_dir+tempfile_name+'tes_intersect.bed'
		bedtools_cmd = ('bedtools', 'intersect', '-a', tes_bedfile, '-b', args.threeprime)
		pipettor.Popen([bedtools_cmd], 'w', stdout=bedtools_output)

		# select reads that contain ends
		precollapse = args.output+'tes_supported.bed' # filename of 3' end-supported, corrected reads
		select_from_bed(bedtools_output, query, precollapse)
		intermediate += [tes_bedfile, precollapse]

	if args.range: # for collapse_range, make sure the minimum number of reads is present
		count = len(open(precollapse).readlines())
		if count < min_reads:
			sys.stderr.write('ERROR: not enough reads in', precollapse)
			sys.exit(1)

	if args.annotation_reliant or args.remove_internal_priming:
		if not args.generate_map:
			args.generate_map = True
		if (not args.annotation_reliant and args.remove_internal_priming) or args.annotation_reliant == 'generate' or not args.annotated_bed:
			if not os.path.exists(args.gtf):
				if not args.gtf:
					sys.stderr.write('Please specify annotated gtf with -f for --annotation_reliant generate\n')
				else:
					sys.stderr.write('GTF file path does not exist\n')
				sys.exit(1)

			if not args.quiet:
				sys.stderr.write('Making transcript fasta using annotated gtf and genome sequence\n')
			args.annotated_bed = args.output+'annotated_transcripts.bed'
			# gtf to bed
			gtf_to_bed(args.annotated_bed, args.gtf, include_gene=True)

		if (not args.annotation_reliant and args.remove_internal_priming) or args.annotation_reliant == 'generate':
			# get transcript sequences
			args.annotation_reliant = args.output+'annotated_transcripts.fa'
			bed_to_sequence(query=args.output+'annotated_transcripts.bed', genome=args.genome,
				outfilename=args.annotation_reliant)


		# minimap (results are piped into count_sam_transcripts.py)
		mm2_cmd = ['minimap2', '-a', '-t', str(args.threads), '-N', '4', '--MD', args.annotation_reliant] + args.reads #'--split-prefix', 'minimap2transcriptomeindex', ##doesn't work with MD tag

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
		if args.remove_internal_priming:
			count_cmd += ['--remove_internal_priming', '--intprimingthreshold', str(args.intprimingthreshold),
						  '--intprimingfracAs', str(args.intprimingfracAs), '--transcriptomefasta',
						  args.annotation_reliant]

		if not args.quiet:
			sys.stderr.write('Aligning reads to reference transcripts\n')
			sys.stderr.write('Counting supporting reads for annotated transcripts\n')
			sys.stderr.write(' '.join(mm2_cmd) + '\n')
			sys.stderr.write(' '.join(count_cmd) + '\n')
		pipettor.run([mm2_cmd, count_cmd])

		if not args.quiet:
			sys.stderr.write('Setting up unassigned reads for flair-collapse novel isoform detection\n')
		# match counts
		counts_file = args.output+'annotated_transcripts.alignment.counts'
		supported_bed = args.output+'annotated_transcripts.supported.bed'
		match_counts(counts_file=counts_file, output_file=supported_bed, bed=args.annotated_bed, min_reads=args.support)

		# get the unassigned reads separately
		subset_reads = args.output+'unassigned.fasta'
		# NOTE: min_reads is passed but not used
		subset_unassigned_reads(readmap=args.output+'annotated_transcripts.isoform.read.map.txt',
                        query=precollapse, support=str(min_reads), output=args.output+'unassigned.bed',
			fastx=args.reads, outfa=subset_reads)

		# TODO: Get rid of this args renaming!
		precollapse = args.output+'unassigned.bed'
		args.reads = [subset_reads]
		intermediate += [subset_reads, precollapse]


	# TODO: collapse_isoforms_precise uses pool and map, which makes it difficult to capture in a function
	collapse_cmd = ['collapse_isoforms_precise.py', '-q', precollapse, '-t', str(args.threads),
		'-m', str(args.max_ends), '-w', str(args.end_window), '-n', args.no_redundant,
		'-o', args.output+'firstpass.unfiltered.bed']
	if args.gtf and not args.no_end_adjustment:
		collapse_cmd += ['-f', args.gtf]
	if args.isoformtss:
		collapse_cmd += ['-i']
	if args.support < 1:
		collapse_cmd += ['-s', str(args.support)]
	if args.quiet:
		collapse_cmd += ['--quiet']
	print(" ".join((str(a) for a in collapse_cmd)), file=sys.stderr)
	pipettor.run(collapse_cmd)


	# filtering out subset isoforms with insufficient support
	keep_extra_column = False
	if float(args.support) < 1:
		keep_extra_column = True
	filter_collapsed_isoforms(queryfile=args.output+'firstpass.unfiltered.bed', mode=args.filter,
		tol=args.end_window, outfile=args.output+'firstpass.bed', keep_extra_column=keep_extra_column)
	intermediate += [args.output+'firstpass.unfiltered.bed']

	# rename first-pass isoforms to annotated transcript IDs if they match
	# if args.gtf:
	if not args.quiet:
		sys.stderr.write('Renaming isoforms using gtf\n')
	identify_gene_isoform(query=args.output+'firstpass.bed', gtf=args.gtf,
		outfilename=args.output+'firstpass.named.bed', annotation_reliant=args.annotation_reliant)

	# TODO: get rid of this renaming (maybe by copying, or otherwise by having a currentbed variable)
	os.rename(args.output+'firstpass.named.bed', args.output+'firstpass.bed')

	# if we want a certain read fraction to support the isoform (instead of a number of reads)
	if float(args.support) < 1:
		filter_isoforms_by_proportion_of_gene_expr(isoforms=args.output+'firstpass.bed',
			outfilename=args.output+'firstpass.filtered.bed', support=args.support)
		os.rename(args.output+'firstpass.filtered.bed', args.output+'firstpass.bed')

	# get the isoform sequences for the first pass we just did
	bed_to_sequence(query=args.output+'firstpass.bed', genome=args.genome, outfilename=args.output+'firstpass.fa')

	# reassign reads to first-pass isoforms
	if not args.quiet:
		sys.stderr.write('Aligning reads to first-pass isoform reference\n')
#	align_files = []  # TODO: this doesn't get filled
	alignout = args.temp_dir + tempfile_name + 'firstpass.'

	# minimap
	if args.mm2_args:
		args.mm2_args = args.mm2_args.split(',')
	mm2_cmd = ['minimap2', '-a', '-t', str(args.threads), '-N', '4', '--MD'] + args.mm2_args + [args.output+'firstpass.fa'] + args.reads

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
	if args.remove_internal_priming:
		count_cmd += ['--remove_internal_priming', '--intprimingthreshold', str(args.intprimingthreshold),
					  '--intprimingfracAs', str(args.intprimingfracAs), '--transcriptomefasta', args.annotation_reliant]
	if not args.quiet:
		sys.stderr.write('Aligning reads to firstpass transcripts\n')
		sys.stderr.write('Counting supporting reads for firstpass transcripts\n')
		sys.stderr.write(' '.join(mm2_cmd) + '\n')
		sys.stderr.write(' '.join(count_cmd) + '\n')
	pipettor.run([mm2_cmd, count_cmd])

	if not args.quiet:
		sys.stderr.write('Filtering isoforms by read coverage\n')

	# match counts
	mc_output = args.output+'isoforms.bed'
	firstpassbed = args.output+'firstpass.bed'
	isoform_file = False
	# TODO: make a test that shows this has different output
	if args.generate_map or args.annotation_reliant:
		isoform_file = args.output+'isoform.read.map.txt'
	match_counts(counts_file=count_file, output_file=mc_output, bed=firstpassbed, min_reads=min_reads,
	      isoform_file = isoform_file)

	if args.annotation_reliant:
		# filter the collapsed isoforms from the annotation
		filter_output = args.output+'isoforms.filter.bed'
		filter_collapsed_isoforms_from_annotation(support=min_reads,
			queryfile=mc_output,
			map_i=args.output+'isoform.read.map.txt',
			annotation=args.output+'annotated_transcripts.supported.bed',
			map_a=args.output+'annotated_transcripts.isoform.read.map.txt',
			outputfile=filter_output,
			new_map=args.output+'combined.isoform.read.map.txt')
		os.rename(filter_output, mc_output)

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

		# convert isoform bed to gtf
		# if args.gtf:
		bed_to_gtf(query=args.output+'isoforms.bed', outputfile=args.output+'isoforms.gtf')

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
