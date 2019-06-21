""" ADT, CMS """

import sys, argparse, subprocess, os, tempfile

if len(sys.argv) > 1 and sys.argv[1] == 'align':
	mode = 'align'
elif len(sys.argv) > 1 and sys.argv[1] == 'correct':
	mode = 'correct'
elif len(sys.argv) > 1 and sys.argv[1] == 'collapse':
	mode = 'collapse'
elif len(sys.argv) > 1 and sys.argv[1] == 'quantify':
	mode = 'quantify'
elif len(sys.argv) > 1 and sys.argv[1] == 'diffExp':
	mode = 'diffExp'
elif len(sys.argv) > 1 and sys.argv[1] == '--version':
	sys.stderr.write('FLAIR v1.4.0\n')
	sys.exit(1)
else:
	sys.stderr.write('usage: python flair.py <mode> --help \n')
	sys.stderr.write('modes: align, correct, collapse, quantify, diffExp\n')
	sys.exit(1)

path = sys.argv[0][:sys.argv[0].rfind('/')+1] if '/' in sys.argv[0] else ''
if mode == 'align':
	parser = argparse.ArgumentParser(description='flair-align parse options', \
		usage='python flair.py align -g genome.fa -r <reads.fq>|<reads.fa> [options]')
	parser.add_argument('align')
	required = parser.add_argument_group('required named arguments')
	required.add_argument('-r', '--reads', action='store', dest='r', \
		type=str, required=True, help='FastA/FastQ files of raw reads')
	required.add_argument('-g', '--genome', action='store', dest='g', \
		type=str, required=True, help='FastA of reference genome')
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
	args = parser.parse_args()

	if args.m[-8:] != 'minimap2':
		if args.m[-1] == '/':
			args.m += 'minimap2'
		else:
			args.m += '/minimap2'

	try:
		if args.n and subprocess.call([args.m, '-ax', 'splice', '-uf', '-k14', '-t', args.t, \
			'--secondary=no', args.g, args.r], stdout=open(args.o+'.sam', 'w')):
			sys.exit(1)
		elif subprocess.call([args.m, '-ax', 'splice', '-t', args.t, '--secondary=no', args.g, args.r], \
				stdout=open(args.o+'.sam', 'w')):
			sys.exit(1)
	except:
		sys.stderr.write('Possible minimap2 error, specify executable path with -m\n')
		sys.exit(1)

	if args.p and subprocess.call([sys.executable, path+'bin/sam_to_psl.py', args.o+'.sam', \
		args.o+'.psl', args.c]):
		sys.exit(1)
	
	sys.stderr.write('Converting sam output to bed\n')
	if subprocess.call([args.sam, 'view', '-h', '-Sb', '-@', args.t, args.o+'.sam'], \
		stdout=open(args.o+'.unsorted.bam', 'w')):  # calls samtools view, exit if an error code that != 0 results
		sys.stderr.write('Possible issue with samtools executable\n')
		sys.exit(1)

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
		usage='python flair.py correct -c chromsizes.tsv -q query.bed12 [-f annotation.gtf]v[-j introns.tab] [options]')
	parser.add_argument('correct')
	required = parser.add_argument_group('required named arguments')
	atleastone = parser.add_argument_group('at least one of the following arguments is required')
	required.add_argument('-q', '--query', type=str, default='', required=True, \
		action='store', dest='q', help='uncorrected bed12 file')
	required.add_argument('-c', '--chromsizes', type=str, required=True, \
		action='store', dest='c', default='', help='chromosome sizes tab-separated file')
	required.add_argument('-g', '--genome', action='store', dest='g', \
		type=str, required=True, help='FastA of reference genome')
	atleastone.add_argument('-j', '--shortread', action='store', dest='j', type=str, default='', \
		help='bed format splice junctions from short-read sequencing')
	atleastone.add_argument('-f', '--gtf', default='', \
		action='store', dest='f', help='GTF annotation file')
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

	if subprocess.call([sys.executable, path+'bin/bed_to_psl.py', args.c, args.o+'_all_corrected.bed', \
		args.o+'_all_corrected.psl']):
		sys.exit(1)

elif mode == 'collapse':
	parser = argparse.ArgumentParser(description='flair-collapse parse options', \
		usage='python flair.py collapse -g genome.fa -r <reads.fq>/<reads.fa> -q <query.psl>|<query.bed12> [options]')
	parser.add_argument('collapse')
	required = parser.add_argument_group('required named arguments')
	required.add_argument('-r', '--reads', action='store', dest='r', \
		type=str, required=True, help='FastA/FastQ files of raw reads')
	required.add_argument('-q', '--query', type=str, default='', required=True, \
		action='store', dest='q', help='PSL file of aligned/corrected reads')
	required.add_argument('-g', '--genome', action='store', dest='g', \
		type=str, required=True, help='FastA of reference genome')
	parser.add_argument('-f', '--gtf', default='', action='store', dest='f', \
		help='GTF annotation file, used for identifying annotated isoforms')
	parser.add_argument('-m', '--minimap2', type=str, default='minimap2', \
		action='store', dest='m', help='path to minimap2 if not in $PATH')
	parser.add_argument('-t', '--threads', type=int, \
		action='store', dest='t', default=4, help='minimap2 number of threads (4)')
	parser.add_argument('-p', '--promoters', action='store', dest='p', default='', \
		help='promoter regions bed file to identify full-length reads')
	parser.add_argument('-b', '--bedtools', action='store', dest='b', default='bedtools', \
		help='bedtools executable path, provide if promoter regions specified')
	parser.add_argument('-sam', '--samtools', action='store', dest='sam', default='samtools', \
		help='samtools executable path if not in $PATH')
	parser.add_argument('-w', '--window', default='100', action='store', dest='w', \
		help='window size for comparing TSS/TES (100)')
	parser.add_argument('-s', '--support', default='3', action='store', dest='s', \
		help='minimum number of supporting reads for an isoform (3)')
	parser.add_argument('--quality', type=int, action='store', dest='quality', default=1, \
		help='minimum MAPQ of read assignment to an isoform. If using salmon, all alignments are used (1)')
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
	parser.add_argument('-e', '--filter', default='default', action='store', dest='e', \
		help='''Report options include:
		default--full-length isoforms only;
		comprehensive--default set + subset isoforms;
		ginormous--comprehensive set + single exon subset isoforms''')
	parser.add_argument('--keep_intermediate', default=False, action='store_true', dest='keep_intermediate', \
		help='''specify if intermediate files are to be kept for debugging. 
		Intermediate files include: promoter-supported reads file, 
		read assignments to firstpass isoforms (default: not specified)''')
	parser.add_argument('--salmon', type=str, action='store', dest='salmon', \
		default='', help='Path to salmon executable, specify if salmon quantification is desired')
	parser.add_argument('-o', '--output', default='flair.collapse', \
		action='store', dest='o', help='output file name base for FLAIR isoforms (default: flair.collapse)')
	args = parser.parse_args()

	if args.m[-8:] != 'minimap2':
		if args.m[-1] == '/':
			args.m += 'minimap2'
		else:
			args.m += '/minimap2'

	args.t, args.quality = str(args.t), str(args.quality)

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
			args.q, args.o+'.promotersupported.psl'])
		precollapse = args.o+'.promotersupported.psl'  # filename of promoter-supported, corrected reads
		intermediate += [args.o+'.tss.bed', args.o+'.promotersupported.psl', args.o+'.promoter_intersect.bed']

	collapse_cmd = [sys.executable, path+'bin/collapse_isoforms_precise.py', '-q', precollapse, \
			'-m', str(args.max_ends), '-w', args.w, '-s', '1', '-n', args.n, '-o', args.o+'.firstpass.unfiltered.psl']
	if args.f:
		collapse_cmd += ['-f', args.f]
	if args.i:
		collapse_cmd += ['-i']
	if subprocess.call(collapse_cmd):
		sys.exit(1)

	sys.stderr.write('Filtering isoforms\n')  # filter more
	if subprocess.call([sys.executable, path+'bin/filter_collapsed_isoforms.py', args.o+'.firstpass.unfiltered.psl', \
		args.e, args.o+'.firstpass.psl', args.w]):
		sys.exit(1)

	if args.f:
		sys.stderr.write('Renaming isoforms\n')
		if subprocess.call([sys.executable, path+'bin/identify_gene_isoform.py', \
			args.o+'.firstpass.psl', args.f, args.o+'.firstpass.named.psl']):
			sys.exit(1)
		subprocess.call(['mv', args.o+'.firstpass.named.psl', args.o+'.firstpass.psl'])

	if subprocess.call([sys.executable, path+'bin/psl_to_sequence.py', args.o+'.firstpass.psl', \
		args.g, args.o+'.firstpass.fa']):
		sys.exit(1)

	sys.stderr.write('Aligning reads to first-pass isoform reference\n')
	reads_files = args.r.split(',')
	count_files, align_files = [], []
	try:
		for r in reads_files:  # align reads to first-pass isoforms
			alignout = tempfile.NamedTemporaryFile().name+'.firstpass'
			if args.salmon:
				if subprocess.call([args.m, '-a', '-t', args.t, args.o+'.firstpass.fa', r], \
					stdout=open(alignout+'.sam', "w")):
					sys.exit(1)
				if subprocess.call([args.sam, 'view', '-F', '4', '-h', '-S', alignout+'.sam'], \
					stdout=open(alignout+'.mapped.sam')):
					sys.exit(1)
				subprocess.call(['mv', alignout+'.mapped.sam', alignout+'.sam'])
				subprocess.call([args.salmon, 'quant', '-t', args.o+'.firstpass.fa', '-o',  alignout+'.salmon', \
					'-p', args.t, '-l', 'U', '-a', alignout+'.sam'], stderr=open('salmon_stderr.txt', 'w'))
				
				count_files += [alignout+'.salmon/quant.sf']
				subprocess.call(['rm', 'salmon_stderr.txt'])
				align_files += [alignout+'.sam', alignout+'.salmon/quant.sf']
				temporary += [alignout+'.salmon']
			else:
				if subprocess.call([args.m, '-a', '-t', args.t, '--secondary=no', \
					args.o+'.firstpass.fa', r], stdout=open(alignout+'.sam', "w")):
					sys.exit(1)
				subprocess.call([args.sam, 'view', '-q', args.quality, '-h', '-S', alignout+'.sam'], \
					stdout=open(alignout+'.q1.sam', 'w'))
				subprocess.call([sys.executable, path+'bin/count_sam_genes.py', alignout+'.q1.sam', \
					alignout+'.q1.counts'])
				count_files += [alignout+'.q1.counts']
				align_files = [alignout+'.sam', alignout+'.q1.sam']
	except:
		sys.stderr.write('Possible minimap2/samtools error, specify paths or make sure they are in $PATH\n')
		sys.exit(1)

	subprocess.call([sys.executable, path+'bin/combine_counts.py'] + count_files + [args.o+'.firstpass.q1.counts'])
	subprocess.call(['rm'] + count_files)

	sys.stderr.write('Filtering isoforms by read coverage\n')
	subprocess.call([sys.executable, path+'bin/match_counts.py', args.o+'.firstpass.q1.counts', \
		args.o+'.firstpass.psl', args.s, args.o+'.isoforms.psl'])
	subprocess.call([sys.executable, path+'bin/psl_to_sequence.py', args.o+'.isoforms.psl', \
		args.g, args.o+'.isoforms.fa'])
	subprocess.call([sys.executable, path+'bin/psl_to_gtf.py', args.o+'.isoforms.psl'], \
		stdout=open(args.o+'.isoforms.gtf', 'w'))

	if args.stringent:  # filtering for isoforms from the high-confidence set with full-length supporting reads
		sys.stderr.write('Applying stringent filtering\n')
		map_files, alignment_psls = [], []
		for r in reads_files:  # align reads to high-confidence isoforms
			alignout = tempfile.NamedTemporaryFile().name+'.stringent'
			subprocess.call([args.m, '-a', '-t', args.t, '--secondary=no', \
				args.o+'.isoforms.fa', r], stdout=open(alignout+'.sam', 'w'))
			subprocess.call([args.sam, 'view', '-q', '1', '-h', '-S', alignout+'.sam'], \
				stdout=open(alignout+'.q1.sam', "w"))
			subprocess.call([sys.executable, path+'bin/sam_to_psl.py', alignout+'.q1.sam', alignout+'.q1.sam.psl'])
			subprocess.call([sys.executable, path+'bin/sam_to_map.py', alignout+'.q1.sam', alignout+'.q1.map'])
			map_files += [alignout+'.q1.map']
			alignment_psls += [alignout+'.q1.sam.psl']
			align_files += [alignout+'.sam', alignout+'.q1.sam']

		subprocess.call(['cat'] + map_files, stdout=open(args.o+'.stringent.map', 'w'))
		subprocess.call(['cat'] + alignment_psls, stdout=open(args.o+'.stringent.q1.sam.psl', 'w'))
		intermediate += [args.o+'.stringent.q1.sam.psl', args.o+'.stringent.map']
		subprocess.call(['rm'] + map_files + alignment_psls)
		subprocess.call([sys.executable, path+'bin/filter_stringent_support.py', args.o+'.isoforms.psl', \
			args.o+'.stringent.q1.sam.psl', args.s, args.o+'.isoforms.stringent.psl',  \
			args.o+'.isoforms.stringent.map'])
		subprocess.call([sys.executable, path+'bin/psl_to_sequence.py', args.o+'.isoforms.stringent.psl', \
			args.g, args.o+'.isoforms.stringent.fa'])
		subprocess.call([sys.executable, path+'bin/psl_to_gtf.py', args.o+'.isoforms.stringent.psl'], \
			stdout=open(args.o+'.isoforms.stringent.gtf', 'w'))

	subprocess.call(['rm', '-rf'] + temporary)
	if not args.keep_intermediate:
		sys.stderr.write('Removing intermediate files/done\n')
		subprocess.call(['rm', args.o+'.firstpass.unfiltered.psl'])
		subprocess.call(['rm', args.o+'.firstpass.fa', args.o+'.firstpass.q1.counts', args.o+'.firstpass.psl'])
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
	parser.add_argument('-t', '--threads', type=str, \
		action='store', dest='t', default='4', help='minimap2 number of threads (4)')
	parser.add_argument('-sam', '--samtools', action='store', dest='sam', default='samtools', \
		help='specify a samtools executable path if not in $PATH if --quality is also used')
	parser.add_argument('--quality', type=str, action='store', dest='quality', default='0', \
		help='minimum MAPQ of read assignment to an isoform. If using salmon, all alignments are used (0)')
	parser.add_argument('-o', '--output', type=str, action='store', dest='o', \
		default='counts_matrix.tsv', help='Counts matrix output name (counts_matrix.tsv)')
	parser.add_argument('--salmon', type=str, action='store', dest='salmon', \
		default='', help='Path to salmon executable, specify if salmon quantification is desired')
	parser.add_argument('--tpm', action='store_true', dest='tpm', default=False, \
		help='specify this flag to output additional file with expression in TPM')
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

	samData = list()
	with codecs.open(args.r, "r", encoding='utf-8', errors='ignore') as lines:
		for line in lines:

			cols = line.rstrip().split('\t')
			if len(cols)<4:
				sys.stderr.write('Expected 4 columns in manifest.tsv, got %s. Exiting.\n' % len(cols))
				sys.exit(1)
			sample, group, batch, readFile = cols
			readFileRoot = tempfile.NamedTemporaryFile().name
			# readFileRoot = readFile[readFile.rfind('/')+1:]
			samData.append(cols + [readFileRoot + '.sam'])
		
		samData.sort(key=lambda x: x[1], reverse=True)
		for num,sample in enumerate(samData,0):
			sys.stderr.write("Step 1/3. Aligning sample %s_%s: %s/%s \r" % (sample[0],sample[2],num+1,len(samData)))
			mm2_command = [args.m, '-a', '-t', args.t, args.i, sample[-2]]
			if not args.salmon:
				mm2_command += ['--secondary=no']
			try:
				if subprocess.call(mm2_command, stdout=open(sample[-1], 'w'), \
					stderr=open(sample[-1]+'.mm2_Stderr.txt', 'w')):
					sys.stderr.write('Check stderr files\n')
					sys.exit(1)
			except:
				sys.stderr.write('Possible minimap2 error, specify executable path with -m\n')
				sys.exit(1)
			sys.stderr.flush()
			if args.quality != '0':
				if subprocess.call([args.sam, 'view', '-q', args.quality, '-h', '-S', sample[-1]], \
					stdout=open(sample[-1]+'.qual.sam', 'w')):
					sys.exit(1)
				subprocess.call(['mv', sample[-1]+'.qual.sam', sample[-1]])
			subprocess.call([sys.executable, path+'bin/sam_to_map.py', sample[-1], sample[-1]+'.map'])
	countData = dict()
	for num,data in enumerate(samData):
		sample, group, batch, readFile, samOut = data
		sys.stderr.write("Step 2/3. Quantifying isoforms for sample %s_%s: %s/%s \r" % (sample,batch,num+1,len(samData)))

		if not args.salmon:
			with open(samOut) as lines:
				for line in lines:
					if line[0] == "@": continue
					cols = line.split()
					aFlag,iso,mapq = cols[1],cols[2],cols[4]

					if aFlag == "16" or aFlag == "0": pass
					else: continue
					if iso not in countData: countData[iso] = np.zeros(len(samData))
					countData[iso][num] += 1 
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
			subprocess.call(['rm', '-rf', samOut[:-4]+'.salmon'])
		subprocess.call(['rm', samOut])
		sys.stderr.flush()

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
		subprocess.call([sys.executable, path+'bin/fasta_seq_lengths.py', args.i, args.i+'.sizes'])
		subprocess.call([sys.executable, path+'bin/counts_to_tpm.py', args.i+'.sizes', args.o, args.o+'.tpm'])

elif mode == 'diffExp':
	parser = argparse.ArgumentParser(description='flair-diffExp parse options', \
		usage='python flair.py diffExp -q count_matrix.tsv --out_dir out_dir [options]')
	parser.add_argument('diffExp')
	required = parser.add_argument_group('required named arguments')
	
	required.add_argument('-q', '--count_matrix', action='store', dest='q', \
		type=str, required=True, help='Tab-delimited isoform count matrix from flair quantify module.')
	required.add_argument('-o', '--out_dir', action='store', dest='o', \
		type=str, required=True, help='Output directory for tables and plots.')
	required.add_argument('-t', '--threads', action='store', dest='t', \
		type=int, required=False, default=4, help='Number of threads for parallel DRIM-Seq.')
	parser.add_argument('-e', '--exp_thresh', action='store', dest='e', type=int, required=False, \
		default=10, help='Read count expression threshold. Isoforms in which \
		both conditions contain fewer than E reads are filtered out (Default E=10)')

	args = parser.parse_args()

	scriptsBin = "/".join(os.path.realpath(__file__).split("/")[:-1])  + "/bin/"
	runDE      = scriptsBin + "deFLAIR.py"


	subprocess.call([sys.executable, '-W ignore', runDE, '--filter', str(args.e), '--threads', str(args.t), '--outDir', args.o, '--matrix', args.q])




