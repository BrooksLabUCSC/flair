""" ADT, CMS """

import sys, argparse, subprocess, os, tempfile

if len(sys.argv) > 1 and (sys.argv[1].lower() == 'align' or sys.argv[1] == '1'):
	mode = 'align'
elif len(sys.argv) > 1 and (sys.argv[1].lower() == 'correct' or sys.argv[1] == '2'):
	mode = 'correct'
elif len(sys.argv) > 1 and (sys.argv[1].lower() == 'collapse' or sys.argv[1] == '3'):
	mode = 'collapse'
elif len(sys.argv) > 1 and (sys.argv[1].lower() == 'quantify' or sys.argv[1] == '4'):
	mode = 'quantify'
elif len(sys.argv) > 1 and (sys.argv[1].lower() == 'diffexp' or sys.argv[1] == '5'):
	mode = 'diffExp'
elif len(sys.argv) > 1 and (sys.argv[1].lower() == 'diffsplice' or sys.argv[1] == '6'):
	mode = 'diffSplice'
elif len(sys.argv) > 1 and sys.argv[1] == '--version':
	sys.stderr.write('FLAIR v1.4.0\n')
	sys.exit(1)
else:
	sys.stderr.write('usage: python flair.py <mode> --help \n')
	sys.stderr.write('modes: align, correct, collapse, quantify, diffExp, diffSplice\n')
	sys.exit(1)

path = '/'.join(os.path.realpath(__file__).split("/")[:-1])+'/'
if mode == 'align':
	parser = argparse.ArgumentParser(description='flair-align parse options', \
		usage='python flair.py align -g genome.fa -r <reads.fq>|<reads.fa> [options]')
	parser.add_argument('align')
	required = parser.add_argument_group('required named arguments')
	required.add_argument('-r', '--reads', action='store', dest='r', \
		nargs='+', type=str, required=True, help='FastA/FastQ files of raw reads')
	required.add_argument('-g', '--genome', action='store', dest='g', \
		type=str, required=True, help='FastA of reference genome, can be minimap2 indexed')
	parser.add_argument('-m', '--minimap2', type=str, default='minimap2', \
		action='store', dest='m', help='path to minimap2 if not in $PATH')
	parser.add_argument('-o', '--output', \
		action='store', dest='o', default='flair.aligned', help='output file name base (default: flair.aligned)')
	parser.add_argument('-t', '--threads', type=str, \
		action='store', dest='t', default='4', help='minimap2 number of threads (4)')
	parser.add_argument('-sam', '--samtools', action='store', dest='sam', default='samtools', \
		help='samtools executable path if not in $PATH')
	parser.add_argument('-c', '--chromsizes', type=str, action='store', dest='c', default='', \
		help='''chromosome sizes tab-separated file, used for converting sam to genome-browser
		compatible psl file''')
	parser.add_argument('-n', '--nvrna', action='store_true', dest='n', default=False, \
		help='specify this flag to use native-RNA specific alignment parameters for minimap2')
	parser.add_argument('-p', '--psl', action='store_true', dest='p', \
		help='also output sam-converted psl')
	parser.add_argument('-v1.3', '--version1.3', action='store_true', dest='v', \
		help='specify if samtools version 1.3+')
	parser.add_argument('--quality', type=int, action='store', dest='quality', default=1, \
		help='minimum MAPQ of read alignment to the genome (1)')
	args = parser.parse_args()

	args.quality = str(args.quality)
	if args.m[-8:] != 'minimap2':
		if args.m[-1] == '/':
			args.m += 'minimap2'
		else:
			args.m += '/minimap2'

	try:
		if args.n:
			if subprocess.call([args.m, '-ax', 'splice', '-uf', '-k14', '-t', args.t, \
			'--secondary=no', args.g] + args.r, stdout=open(args.o+'.sam', 'w')):
				sys.exit(1)
		elif subprocess.call([args.m, '-ax', 'splice', '-t', args.t, '--secondary=no', args.g]+args.r, \
				stdout=open(args.o+'.sam', 'w')):
			sys.exit(1)
	except:
		sys.stderr.write('Possible minimap2 error, specify executable path with -m\n')
		sys.exit(1)

	if args.p and subprocess.call([sys.executable, path+'bin/sam_to_psl.py', args.o+'.sam', \
		args.o+'.psl', args.c]):
		sys.exit(1)

	sys.stderr.write('Converting sam output to bed\n')
	if args.quality != '0':
		if subprocess.call([args.sam, 'view', '-q', args.quality, '-h', '-S', args.o+'.sam'], \
		stdout=open(args.o+'.q.sam', 'w')):
			sys.stderr.write('Possible issue with samtools executable\n')
			sys.exit(1)
		subprocess.call(['mv', args.o+'.q.sam', args.o+'.sam'])

	if subprocess.call([args.sam, 'view', '-h', '-Sb', '-@', args.t, args.o+'.sam'], \
		stdout=open(args.o+'.unsorted.bam', 'w')):  # calls samtools view, exit if an error code that != 0 results
		sys.stderr.write('Possible issue with samtools executable\n')
		sys.exit(1)

	if not args.v:  # samtools version is < 1.3 or unspecified --> detect version
		ver = subprocess.Popen([args.sam], stderr=subprocess.PIPE, universal_newlines=True)
		for line in ver.stderr:
			if 'Version:' in line:
				v = line.rstrip()[line.find('Version:')+9:line.find('Version:')+12]
				try:
					if float(v) >= 1.3:
						sys.stderr.write('Samtools version >= 1.3 detected\n')
						args.v = True
						break
				except:
					sys.stderr.write('Could not detect samtools version, will assume < 1.3\n')

	if args.v:  # samtools verison 1.3+
		subprocess.call([args.sam, 'sort', '-@', args.t, args.o+'.unsorted.bam', '-o', args.o+'.bam'])
	elif subprocess.call([args.sam, 'sort', '-@', args.t, args.o+'.unsorted.bam', args.o]):
		sys.stderr.write('If using samtools v1.3+, please specify -v1.3 argument\n')
		sys.exit(1)

	subprocess.call([args.sam, 'index', args.o+'.bam'])
	subprocess.call([sys.executable, path+'bin/bam2Bed12.py', '-i', args.o+'.bam'], stdout=open(args.o+'.bed', 'w'))
	subprocess.call(['rm', args.o+'.unsorted.bam'])

elif mode == 'correct':
	parser = argparse.ArgumentParser(description='flair-correct parse options', \
		usage='python flair.py correct -q query.bed12 [-f annotation.gtf]v[-j introns.tab] -g genome.fa [options]')
	parser.add_argument('correct')
	required = parser.add_argument_group('required named arguments')
	atleastone = parser.add_argument_group('at least one of the following arguments is required')
	required.add_argument('-q', '--query', type=str, default='', required=True, \
		action='store', dest='q', help='uncorrected bed12 file')
	required.add_argument('-g', '--genome', action='store', dest='g', \
		type=str, required=True, help='FastA of reference genome')
	atleastone.add_argument('-j', '--shortread', action='store', dest='j', type=str, default='', \
		help='bed format splice junctions from short-read sequencing')
	atleastone.add_argument('-f', '--gtf', default='', \
		action='store', dest='f', help='GTF annotation file')
	parser.add_argument('-c', '--chromsizes', type=str, \
		action='store', dest='c', default='', help='chromosome sizes tab-separated file')
	parser.add_argument('-n', '--nvrna', action='store_true', dest='n', default=False, help='specify this flag to keep \
		the strand of a read consistent after correction')
	parser.add_argument('-t', '--threads', type=str, action='store', dest='t', default='4', \
		help='splice site correction script number of threads (4)')
	parser.add_argument('-w', '--window', action='store', dest='w', default='10', \
		help='window size for correcting splice sites (W=10)')
	parser.add_argument('-o', '--output', \
		action='store', dest='o', default='flair', help='output name base (default: flair)')
	parser.add_argument('--print_check', \
		action='store_true', dest='p', default=False, help='Print err.txt with step checking.')

	args = parser.parse_args()
	if not args.j and not args.f:
		sys.stderr.write('Please specify at least one of the -f or -j arguments for correction\n')
		sys.exit(1)
	correction_cmd = [sys.executable, path+'bin/ssCorrect.py', '-i', args.q, \
			'-w', args.w, '-p', args.t, '-o', args.o, '--progress', '-f', args.g]
	if not args.n:
		correction_cmd += ['--correctStrand']
	if args.j:
		correction_cmd += ['-j', args.j]
	if args.f:
		correction_cmd += ['-g', args.f]
	if args.p:
		correction_cmd += ['--print_check']

	if subprocess.call(correction_cmd):
		sys.stderr.write('Correction command did not exit with success status\n')
		sys.exit(1)

	if args.c and subprocess.call([sys.executable, path+'bin/bed_to_psl.py', args.c, args.o+'_all_corrected.bed', \
		args.o+'_all_corrected.psl']):
		sys.exit(1)

elif mode == 'collapse':
	parser = argparse.ArgumentParser(description='flair-collapse parse options', \
		usage='python flair.py collapse -g genome.fa -r <reads.fq>/<reads.fa> -q <query.psl>|<query.bed> [options]')
	parser.add_argument('collapse')
	required = parser.add_argument_group('required named arguments')
	required.add_argument('-r', '--reads', action='store', dest='r', nargs='+', \
		type=str, required=True, help='FastA/FastQ files of raw reads')
	required.add_argument('-q', '--query', type=str, default='', required=True, \
		action='store', dest='q', help='PSL file of aligned/corrected reads')
	required.add_argument('-g', '--genome', action='store', dest='g', \
		type=str, required=True, help='FastA of reference genome')
	parser.add_argument('-f', '--gtf', default='', action='store', dest='f', \
		help='GTF annotation file, used for renaming FLAIR isoforms to annotated isoforms and adjusting TSS/TESs')
	parser.add_argument('-m', '--minimap2', type=str, default='minimap2', \
		action='store', dest='m', help='path to minimap2 if not in $PATH')
	parser.add_argument('-t', '--threads', type=int, \
		action='store', dest='t', default=4, help='minimap2 number of threads (4)')
	parser.add_argument('-p', '--promoters', action='store', dest='p', default='', \
		help='promoter regions bed file to identify full-length reads')
	parser.add_argument('-b', '--bedtools', action='store', dest='b', default='bedtools', \
		help='bedtools executable path, provide if promoter regions specified and bedtools is not in $PATH')
	parser.add_argument('-sam', '--samtools', action='store', dest='sam', default='samtools', \
		help='samtools executable path if not in $PATH')
	parser.add_argument('-w', '--window', default='100', action='store', dest='w', \
		help='window size for comparing TSS/TES (100)')
	parser.add_argument('-s', '--support', default='3', action='store', dest='s', \
		help='minimum number of supporting reads for an isoform (3)')
	parser.add_argument('--stringent', default=False, action='store_true', dest='stringent', \
		help='''specify if all supporting reads need to be full-length \
		(80%% coverage and spanning 25 bp of the first and last exons)''')
	parser.add_argument('-n', '--no_redundant', default='none', action='store', dest='n', \
		help='''For each unique splice junction chain, report options include:
		none--best TSSs/TESs chosen for each unique set of splice junctions;
		longest--single TSS/TES chosen to maximize length;
		best_only--single most supported TSS/TES used in conjunction chosen (none)''')
	parser.add_argument('-i', '--isoformtss', default=False, action='store_true', dest='i', \
		help='when specified, TSS/TES for each isoform will be determined from supporting reads \
		for individual isoforms (default: not specified, determined at the gene level)')
	parser.add_argument('--max_ends', default=2, action='store', dest='max_ends', \
		help='maximum number of TSS/TES picked per isoform (2)')
	parser.add_argument('--trust_ends', default=False, action='store_true', dest='trust_ends', \
		help='specify if reads are generated from a long read method with minimal fragmentation')
	parser.add_argument('-e', '--filter', default='default', action='store', dest='e', \
		help='''Report options include:
		nosubset--any isoforms that are a proper set of another isoform are removed;
		default--subset isoforms are removed based on support;
		comprehensive--default set + all subset isoforms;
		ginormous--comprehensive set + single exon subset isoforms''')
	parser.add_argument('--quality', type=int, action='store', dest='quality', default=1, \
		help='minimum MAPQ of read assignment to an isoform (1)')
	parser.add_argument('--keep_intermediate', default=False, action='store_true', dest='keep_intermediate', \
		help='''specify if intermediate files are to be kept for debugging.
		Intermediate files include: promoter-supported reads file,
		read assignments to firstpass isoforms (default: not specified)''')
	parser.add_argument('--generate_map', default=False, action='store_true', dest='generate_map', \
		help='''specify this argument to generate a txt file of which reads are assigned to each isoform. 
		note: only works if the quantification method is not using salmon (default: not specified)''')
	parser.add_argument('--salmon', type=str, action='store', dest='salmon', \
		default='', help='Path to salmon executable, specify if salmon quantification is desired')
	parser.add_argument('--temp_dir', default='', action='store', dest='temp_dir', \
		help='directory to put temporary files. use "./" to indicate current directory (default: python tempfile directory)')
	parser.add_argument('-o', '--output', default='flair.collapse', \
		action='store', dest='o', help='output file name base for FLAIR isoforms (default: flair.collapse)')
	args = parser.parse_args()

	if args.m[-8:] != 'minimap2':
		if args.m[-1] == '/':
			args.m += 'minimap2'
		else:
			args.m += '/minimap2'

	args.t, args.quality = str(args.t), str(args.quality)
	if args.trust_ends:
		args.quality = '0'

	if not os.path.exists(args.q):
		sys.stderr.write('Query file path does not exist\n')
		sys.exit(1)
	if os.stat(args.q).st_size == 0:
		sys.stderr.write('Query file is empty\n')
		sys.exit(1)
	suffix = args.q[-3:]  # bed or psl

	intermediate, temporary = [], []
	precollapse = args.q  # filename
	if args.p:
		sys.stderr.write('Filtering out reads without promoter-supported TSS\n')
		if subprocess.call([sys.executable, path+'bin/pull_starts.py', args.q, args.o+'.tss.bed']):
			sys.exit(1)
		if subprocess.call([args.b, 'intersect', '-a', args.o+'.tss.bed', '-b', args.p], \
			stdout=open(args.o+'.promoter_intersect.bed', 'w')):
			sys.exit(1)
		subprocess.call([sys.executable, path+'bin/psl_reads_from_bed.py', args.o+'.promoter_intersect.bed', \
			args.q, args.o+'.promotersupported.'+suffix])
		precollapse = args.o+'.promotersupported.'+suffix  # filename of promoter-supported, corrected reads
		intermediate += [args.o+'.tss.bed', args.o+'.promotersupported.'+suffix, args.o+'.promoter_intersect.bed']

	collapse_cmd = [sys.executable, path+'bin/collapse_isoforms_precise.py', '-q', precollapse, \
			'-m', str(args.max_ends), '-w', args.w, '-n', args.n, '-o', args.o+'.firstpass.unfiltered.'+suffix]
	if args.f:
		collapse_cmd += ['-f', args.f]
	if args.i:
		collapse_cmd += ['-i']
	if subprocess.call(collapse_cmd):
		sys.exit(1)

	sys.stderr.write('Filtering isoforms\n')  # filtering out subset isoforms with fewer reads
	if subprocess.call([sys.executable, path+'bin/filter_collapsed_isoforms.py', \
		args.o+'.firstpass.unfiltered.'+suffix, args.e, args.o+'.firstpass.'+suffix, args.w]):
		sys.exit(1)

	intermediate += [args.o+'.firstpass.unfiltered.'+suffix]

	if args.f:
		sys.stderr.write('Renaming isoforms\n')
		if subprocess.call([sys.executable, path+'bin/identify_gene_isoform.py', \
			args.o+'.firstpass.'+suffix, args.f, args.o+'.firstpass.named.'+suffix]):
			sys.exit(1)
		subprocess.call(['mv', args.o+'.firstpass.named.'+suffix, args.o+'.firstpass.'+suffix])

	if subprocess.call([sys.executable, path+'bin/psl_to_sequence.py', args.o+'.firstpass.'+suffix, \
		args.g, args.o+'.firstpass.fa']):
		sys.exit(1)

	sys.stderr.write('Aligning reads to first-pass isoform reference\n')
	if ',' in args.r[0]:
		args.r = args.r[0].split(',')
	count_files, align_files = [], []

	# align reads to first-pass isoforms
	tempfile_name = tempfile.NamedTemporaryFile().name
	if args.temp_dir == '':
		alignout = tempfile_name+'.firstpass'
	else:
		if not os.path.isdir(args.temp_dir):
			if subprocess.call(['mkdir', args.temp_dir]):
				sys.stderr.write('Could not make temporary directory {}\n'.format(args.temp_dir))
				sys.exit(1)
		alignout = args.temp_dir + '/' + tempfile_name[tempfile_name.rfind('/'):]+'.firstpass'

	try:
		if subprocess.call([args.m, '-a', '-t', args.t, '-N', '4', args.o+'.firstpass.fa'] + args.r, \
			stdout=open(alignout+'.sam', "w")):
			sys.exit(1)
	except Exception as e:
		sys.stderr.write(str(e)+'\n\n\nMinimap2 error, please check that all file, directory, and executable paths exist\n')
		sys.exit(1)

	if args.salmon:
		if subprocess.call([args.sam, 'view', '-F', '4', '-h', '-S', alignout+'.sam'], \
			stdout=open(alignout+'.mapped.sam', 'w')):
			sys.exit(1)
		subprocess.call(['mv', alignout+'.mapped.sam', alignout+'.sam'])
		subprocess.call([args.salmon, 'quant', '-t', args.o+'.firstpass.fa', '-o',  alignout+'.salmon', \
			'-p', args.t, '-l', 'U', '-a', alignout+'.sam'], stderr=open('salmon_stderr.txt', 'w'))
		count_file = alignout+'.salmon/quant.sf'
		align_files += [alignout+'.sam', alignout+'.salmon/quant.sf']
		temporary += [alignout+'.salmon', 'salmon_stderr.txt']
	else:
		if args.quality != '0':
			subprocess.call([args.sam, 'view', '-q', args.quality, '-h', '-S', alignout+'.sam'], \
				stdout=open(alignout+'.q.sam', 'w'))
			align_files += [alignout+'.sam']
		else:
			subprocess.call(['mv', alignout+'.sam', alignout+'.q.sam'])
		count_cmd = [sys.executable, path+'bin/count_sam_transcripts.py', '-s', alignout+'.q.sam', \
			'-o', alignout+'.q.counts', '-t', args.t, '--quality', args.quality]
		if args.stringent:
			count_cmd += ['--stringent', '-i', args.o+'.firstpass.'+suffix]
		if args.trust_ends:
			count_cmd += ['--trust_ends']
		if args.generate_map:
			count_cmd += ['--generate_map', args.o+'.isoform.read.map.txt']
		if subprocess.call(count_cmd):
			sys.exit(1)
		count_file = alignout+'.q.counts'
		align_files += [alignout+'.q.sam']

	subprocess.call([sys.executable, path+'bin/combine_counts.py', count_file, args.o+'.firstpass.q.counts'])

	sys.stderr.write('Filtering isoforms by read coverage\n')
	subprocess.call([sys.executable, path+'bin/match_counts.py', args.o+'.firstpass.q.counts', \
		args.o+'.firstpass.'+suffix, args.s, args.o+'.isoforms.'+suffix])
	subprocess.call([sys.executable, path+'bin/psl_to_sequence.py', args.o+'.isoforms.'+suffix, \
		args.g, args.o+'.isoforms.fa'])
	if args.f:
		subprocess.call([sys.executable, path+'bin/psl_to_gtf.py', args.o+'.isoforms.'+suffix], \
			stdout=open(args.o+'.isoforms.gtf', 'w'))

	subprocess.call(['rm', '-rf'] + temporary + [args.o+'.firstpass.fa', alignout+'.q.counts'])
	if not args.keep_intermediate:
		sys.stderr.write('Removing intermediate files/done\n')
		subprocess.call(['rm', args.o+'.firstpass.q.counts', args.o+'.firstpass.'+suffix])
		subprocess.call(['rm', '-rf'] + align_files + intermediate)

elif mode == 'quantify':
	parser = argparse.ArgumentParser(description='flair-quantify parse options', \
		usage='python flair.py quantify -r reads_manifest.tsv -i isoforms.fa [options]')
	parser.add_argument('quantify')
	required = parser.add_argument_group('required named arguments')
	required.add_argument('-r', '--reads_manifest', action='store', dest='r', type=str, \
		required=True, help='Tab delimited file containing sample id, condition, batch, reads.fq')
	required.add_argument('-i', '--isoforms', action='store', dest='i', \
		type=str, required=True, help='FastA of FLAIR collapsed isoforms')
	parser.add_argument('-m', '--minimap2', type=str, default='minimap2', \
		action='store', dest='m', help='path to minimap2 if not in $PATH')
	parser.add_argument('-t', '--threads', type=int, \
		action='store', dest='t', default=4, help='minimap2 number of threads (4)')
	parser.add_argument('-sam', '--samtools', action='store', dest='sam', default='samtools', \
		help='specify a samtools executable path if not in $PATH if --quality is also used')
	parser.add_argument('--quality', type=int, action='store', dest='quality', default=1, \
		help='''minimum MAPQ of read assignment to an isoform. If using salmon, all alignments are
		used (1)''')
	parser.add_argument('-o', '--output', type=str, action='store', dest='o', \
		default='counts_matrix.tsv', help='Counts matrix output name (counts_matrix.tsv)')
	parser.add_argument('--salmon', type=str, action='store', dest='salmon', \
		default='', help='Path to salmon executable, specify if salmon quantification is desired')
	parser.add_argument('--tpm', action='store_true', dest='tpm', default=False, \
		help='specify this flag to output additional file with expression in TPM')
	parser.add_argument('--trust_ends', default=False, action='store_true', dest='trust_ends', \
		help='specify if reads are generated from a long read method with minimal fragmentation')
	parser.add_argument('--temp_dir', default='', action='store', dest='temp_dir', \
		help='''directory to put temporary files. use "./" to indicate current directory
		(default: python tempfile directory)''')
	args = parser.parse_args()

	try:
		import numpy as np
		import codecs
	except:
		sys.stderr.write('Numpy import error. Please pip install numpy. Exiting.\n')
		sys.exit(1)

	if args.m[-8:] != 'minimap2':
		if args.m[-1] == '/':
			args.m += 'minimap2'
		else:
			args.m += '/minimap2'
	args.t, args.quality = str(args.t), str(args.quality)

	samData = list()
	with codecs.open(args.r, "r", encoding='utf-8', errors='ignore') as lines:
		for line in lines:

			cols = line.rstrip().split('\t')
			if len(cols)<4:
				sys.stderr.write('Expected 4 columns in manifest.tsv, got %s. Exiting.\n' % len(cols))
				sys.exit(1)
			sample, group, batch, readFile = cols

			readFileRoot = tempfile.NamedTemporaryFile().name
			if args.temp_dir != '':
				if not os.path.isdir(args.temp_dir):
					subprocess.call(['mkdir', args.temp_dir])
				readFileRoot = args.temp_dir + '/' + readFileRoot[readFileRoot.rfind('/')+1:]

			# readFileRoot = readFile[readFile.rfind('/')+1:]
			samData.append(cols + [readFileRoot + '.sam'])

		for num,sample in enumerate(samData,0):
			sys.stderr.write("Step 1/3. Aligning sample %s_%s: %s/%s \r" % (sample[0],sample[2],num+1,len(samData)))
			mm2_command = [args.m, '-a', '-N', '4', '-t', args.t, args.i, sample[-2]]
			try:
				if subprocess.call(mm2_command, stdout=open(sample[-1], 'w'), \
					stderr=open(sample[-1]+'.mm2_Stderr.txt', 'w')):
					sys.stderr.write('Check {} file\n'.format(sample[-1]+'.mm2_Stderr.txt'))
					sys.exit(1)
			except:
				sys.stderr.write('''Possible minimap2 error, please check that all file, directory,
					and executable paths exist\n''')
				sys.exit(1)
			subprocess.call(['rm', sample[-1]+'.mm2_Stderr.txt'])
			sys.stderr.flush()

			if args.quality != '0' and not args.trust_ends and not args.salmon:
				if subprocess.call([args.sam, 'view', '-q', args.quality, '-h', '-S', sample[-1]], \
					stdout=open(sample[-1]+'.qual.sam', 'w')):
					sys.exit(1)
				subprocess.call(['mv', sample[-1]+'.qual.sam', sample[-1]])

	countData = dict()
	for num,data in enumerate(samData):
		sample, group, batch, readFile, samOut = data
		sys.stderr.write("Step 2/3. Quantifying isoforms for sample %s_%s: %s/%s \r" % (sample,batch,num+1,len(samData)))

		if not args.salmon:
			count_cmd = [sys.executable, path+'bin/count_sam_transcripts.py', '-s', samOut, \
				'-o', samOut+'.counts.txt', '-t', args.t, '--quality', args.quality]
			if args.trust_ends:
				count_cmd += ['--trust_ends']
			subprocess.call(count_cmd)
			for line in open(samOut+'.counts.txt'):
				line = line.rstrip().split('\t')
				iso, numreads = line[0], line[1]
				if iso not in countData: countData[iso] = np.zeros(len(samData))
				countData[iso][num] = numreads
			# subprocess.call(['rm', samOut+'.counts.txt'])
		else:
			subprocess.call([args.salmon, 'quant', '-t', args.i, '-o', samOut[:-4]+'.salmon', \
				'-p', args.t, '-l', 'U', '-a', samOut], stderr=open('salmon_stderr.txt', 'w'))
			salmonOut = open(samOut[:-4]+'.salmon/quant.sf')
			salmonOut.readline()  # header
			for line in salmonOut:
				line = line.rstrip().split('\t')
				iso, tpm, numreads = line[0], line[3], line[4]
				if iso not in countData: countData[iso] = np.zeros(len(samData))
				if args.tpm:
					countData[iso][num] = tpm
				else:
					countData[iso][num] = numreads
			subprocess.call(['rm', 'salmon_stderr.txt'])
			subprocess.call(['rm', '-r', samOut[:-4]+'.salmon/'])
		sys.stderr.flush()
		subprocess.call(['rm', samOut])

	sys.stderr.write("Step 3/3. Writing counts to {} \r".format(args.o))
	countMatrix = open(args.o,'w')

	countMatrix.write("ids\t%s\n" % "\t".join(["_".join(x[:3]) for x in samData]))

	features = sorted(list(countData.keys()))
	for f in features:
		countMatrix.write("%s\t%s\n" % (f,"\t".join(str(x) for x in countData[f])))

	countMatrix.close()
	sys.stderr.flush()
	sys.stderr.write("\n")

	if args.tpm and not args.salmon:
		subprocess.call([sys.executable, path+'bin/counts_to_tpm.py', args.o, args.o+'.tpm'])

elif mode == 'diffExp':
	parser = argparse.ArgumentParser(description='flair-diffExp parse options', \
		usage='python flair.py diffExp -q counts_matrix.tsv --out_dir out_dir [options]')
	parser.add_argument('diffExp')
	required = parser.add_argument_group('required named arguments')

	required.add_argument('-q', '--counts_matrix', action='store', dest='q', \
		type=str, required=True, help='Tab-delimited isoform count matrix from flair quantify module.')
	required.add_argument('-o', '--out_dir', action='store', dest='o', \
		type=str, required=True, help='Output directory for tables and plots.')
	parser.add_argument('-t', '--threads', action='store', dest='t', \
		type=int, required=False, default=4, help='Number of threads for parallel DRIMSeq.')
	parser.add_argument('-e', '--exp_thresh', action='store', dest='e', type=int, required=False, \
		default=10, help='Read count expression threshold. Isoforms in which \
		both conditions contain fewer than E reads are filtered out (Default E=10)')
	required.add_argument('-of', '--out_dir_force', action='store_true', dest='of', \
		required=False, help='''Specify this argument to force overwriting of files in
		an existing output directory''')

	args = parser.parse_args()

	scriptsBin	= path + "bin/"
	runDE		= scriptsBin + "deFLAIR.py"
	DEcommand	= [sys.executable, '-W ignore', runDE, '--filter', str(args.e), '--threads', \
		str(args.t), '--outDir', args.o, '--matrix', args.q]

	if args.of:
		DEcommand += ['-of']

	subprocess.call(DEcommand)


elif mode == 'diffSplice':
	parser = argparse.ArgumentParser(description='flair-diffSplice parse options', \
		usage='python flair.py diffSplice -i isoforms.bed|isoforms.psl -q counts_matrix.tsv [options]')
	parser.add_argument('diffExp')
	required = parser.add_argument_group('required named arguments')

	required.add_argument('-i', '--isoforms', action='store', dest='i', required=True, \
		type=str, help='isoforms in bed or psl format')
	required.add_argument('-q', '--counts_matrix', action='store', dest='q', \
		type=str, required=True, help='tab-delimited isoform count matrix from flair quantify module')
	parser.add_argument('-o', '--output', action='store', dest='o', default='flair.diffsplice', type=str, \
		required=False, help='output file name base for FLAIR isoforms (default: flair.diffsplice)')
	parser.add_argument('--test', action='store_true', dest='test', \
		required=False, default=False, help='Run DRIMSeq statistical testing')
	parser.add_argument('-t', '--threads', action='store', dest='t', \
		type=int, required=False, default=1, help='Number of threads DRIMSeq (1)')
	parser.add_argument('--drim1', action='store', dest='drim1', type=int, required=False, default=6, \
		help='''The minimum number of samples that have coverage over an AS event inclusion/exclusion
		for DRIMSeq testing; events with too few samples are filtered out and not tested (6)''')
	parser.add_argument('--drim2', action='store', dest='drim2', type=int, required=False, default=3, \
		help='''The minimum number of samples expressing the inclusion of an AS event;
		 events with too few samples are filtered out and not tested (3)''')
	parser.add_argument('--drim3', action='store', dest='drim3', type=int, required=False, default=15, \
		help='''The minimum number of reads covering an AS event inclusion/exclusion for DRIMSeq testing,
		 events with too few samples are filtered out and not tested (15)''')
	parser.add_argument('--drim4', action='store', dest='drim4', type=int, required=False, default=5, \
		help='''The minimum number of reads covering an AS event inclusion for DRIMSeq testing,
		 events with too few samples are filtered out and not tested (5)''')
	parser.add_argument('--batch', action='store_true', dest='batch', required=False, default=False, \
		help='''If specified with --test, DRIMSeq will perform batch correction''')
	parser.add_argument('--conditionA', action='store', dest='conditionA', required=False, default='', \
		help='''Specify one condition corresponding to samples in the counts_matrix to be compared against
		condition2; by default, the first two unique conditions are used''')
	parser.add_argument('--conditionB', action='store', dest='conditionB', required=False, default='', \
		help='''Specify another condition corresponding to samples in the counts_matrix to be compared against
		conditionA''')
	args = parser.parse_args()

	if args.i[-3:].lower() == 'psl':
		subprocess.call([sys.executable, path+'bin/psl_to_bed.py', args.i, args.i+'.bed'])
		args.i = args.i+'.bed'

	subprocess.call([sys.executable, path+'bin/call_diffsplice_events.py', args.i, args.o, args.q])
	subprocess.call([sys.executable, path+'bin/es_as.py', args.i], stdout=open(args.o+'.es.events.tsv','w'))
	subprocess.call([sys.executable, path+'bin/es_as_inc_excl_to_counts.py', args.q, args.o+'.es.events.tsv'], \
		stdout=open(args.o+'.es.events.quant.tsv','w'))
	subprocess.call(['rm', args.o+'.es.events.tsv'])

	if args.test or args.drim1 or args.drim2 or args.drim4 or args.drim4:
		sys.stderr.write('DRIMSeq testing for each AS event type\n')
		drim1, drim2, drim3, drim4 = [str(x) for x in [args.drim1, args.drim2, args.drim3, args.drim4]]
		ds_command = [sys.executable, path+'bin/runDS.py', '--threads', str(args.t), \
			'--drim1', drim1, '--drim2', drim2, '--drim3', drim3, '--drim4', drim4]
		if args.batch:
			ds_command += ['--batch']
		if args.conditionA:
			if not args.conditionB:
				sys.stderr.write('Both conditionA and conditionB must be specified, or both left unspecified\n')
				sys.exit(1)
			ds_command += ['--conditionA', args.conditionA, '--conditionB',  args.conditionB]

		with open(args.o+'.stderr.txt', 'w') as ds_stderr:
			subprocess.call(ds_command + ['--matrix', args.o+'.es.events.quant.tsv', '--prefix', args.o+'.es'], stderr=ds_stderr)
			subprocess.call(ds_command + ['--matrix', args.o+'.alt5.events.quant.tsv', '--prefix', args.o+'.alt5'], stderr=ds_stderr)
			subprocess.call(ds_command + ['--matrix', args.o+'.alt3.events.quant.tsv', '--prefix', args.o+'.alt3'], stderr=ds_stderr)
			subprocess.call(ds_command + ['--matrix', args.o+'.ir.events.quant.tsv', '--prefix', args.o+'.ir'], stderr=ds_stderr)

