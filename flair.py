""" AD Tang 2018 """

import sys, argparse, subprocess, csv
# from subprocess import Popen, PIPE

if len(sys.argv) > 1 and sys.argv[1] == 'align':
	mode = 'align'
elif len(sys.argv) > 1 and sys.argv[1] == 'correct':
	mode = 'correct'
elif len(sys.argv) > 1 and sys.argv[1] == 'collapse':
	mode = 'collapse'
else:
	sys.stderr.write('usage: python flair.py <mode> --help \n')
	sys.stderr.write('modes: align, correct, collapse\n')
	sys.exit()

path = sys.argv[0][:sys.argv[0].rfind('/')+1] if '/' in sys.argv[0] else ''
if mode == 'align':
	parser = argparse.ArgumentParser(description='flair-align parse options', \
				usage='python flair.py align -r <reads.fq>/<reads.fa> -g genome.fa [options]')
	parser.add_argument('align')
	required = parser.add_argument_group('required named arguments')
	required.add_argument('-r', '--reads', action='store', dest='r', \
		type=str, required=True, \
		help='FastA/FastQ files of raw reads')
	required.add_argument('-g', '--genome', action='store', dest='g', \
		type=str, required=True, help='FastA of reference genome')
	parser.add_argument('-m', '--minimap2', type=str, default='minimap2', \
		action='store', dest='m', help='path to minimap2')
	parser.add_argument('-o', '--output', \
		action='store', dest='o', default='', help='output file name')
	parser.add_argument('-t', '--threads', type=str, \
		action='store', dest='t', default='4', help='Minimap2 number of threads (4)')
	parser.add_argument('-c', '--chromsizes', type=str, \
		action='store', dest='c', default='', help='Chromosome sizes tab-separated file')
	args = parser.parse_args()

	args.o = args.o if args.o else args.r[:args.r.rfind('.')+1]+'sam'
	if args.m[-8:] != 'minimap2':
		if args.m[-1] == '/':
			args.m += 'minimap2'
		else:
			args.m += '/minimap2'

	sys.stderr.write('Aligning to the genome with minimap2\n')
	subprocess.call([args.m, '-a', '-t', args.t, '--secondary=no', args.g, args.r], stdout=open(args.o, 'w'))

	if args.c:
		subprocess.call(['python', path+'bin/sam_to_psl.py', args.o, args.o[:-3]+'psl', args.c])
	else:
		subprocess.call(['python', path+'bin/sam_to_psl.py', args.o, args.o[:-3]+'psl'])

elif mode == 'correct':
	parser = argparse.ArgumentParser(description='flair-correct parse options', \
				usage='python flair.py correct -a annotated.gp -g genome.fa -q query.psl [options]')
	parser.add_argument('correct')
	required = parser.add_argument_group('required named arguments')
	required.add_argument('-a', '--annotation', action='store', dest='a', \
		type=str, required=True, \
		help='genepred format splice junctions from annotation')
	parser.add_argument('-s', '--shortread', action='store', dest='s', \
		type=str, default='', \
		help='genepred format splice junctions from short-read sequencing, consider merging with annotated junctions')
	required.add_argument('-g', '--genome', action='store', dest='g', \
		type=str, required=True, help='fasta of reference genome')
	required.add_argument('-q', '--query', type=str, default='', required=True, \
		action='store', dest='q', help='uncorrected PSL file')
	parser.add_argument('-w', '--window', \
		action='store', dest='w', default='10', \
		help='window size for correcting splice sites (10)')
	parser.add_argument('-m', '--mergesize', \
		action='store', dest='m', default='30', help='merge size for gaps within exons (30)')
	parser.add_argument('-o', '--output', \
		action='store', dest='o', default='', help='output file name (correction_output)')
	args = parser.parse_args()
	args.s = args.s if args.s else args.a
	stranded_psl = args.q[:args.q.rfind('.')] + '_strand.psl'
	args.o = args.o if args.o else 'correction_output'

	sys.stderr.write('Inferring strand for reads in PSL\n')
	subprocess.call(['python', path+'bin/infer_strand_for_psl.py', args.q, args.a, stranded_psl, 'gp'])

	sys.stderr.write('Correcting reads in {}\n'.format(stranded_psl))
	subprocess.call(['python', path+'bin/correctSplice.py', '-a', args.a, '-j', args.s, \
	'-g', args.g, '-q', stranded_psl, '-w', args.w, '-m', args.m, '-o', args.o, '-mode', 'closest'])

	sys.stderr.write('Converting output gp to PSL\n')
	subprocess.call(['python', path+'bin/genePredToPSL.py', stranded_psl, args.o+'/corrected.gp',\
		args.o+'/'+stranded_psl[:-4]+'_corrected_novels.psl'])

	sys.stderr.write('Filtering out reads containing unsupported junctions\n')
	subprocess.call(['python', path+'bin/remove_novel.py', args.a, \
		args.o+'/'+stranded_psl[:-4]+'_corrected_novels.psl', \
		args.o+'/'+stranded_psl[:-4]+'_corrected.psl'])
	subprocess.call(['rm', args.o+'/'+stranded_psl[:-4]+'_corrected_novels.psl'])

elif mode == 'collapse':
	parser = argparse.ArgumentParser(description='flair-collapse parse options', \
				usage='python flair.py correct -r <reads.fq>/<reads.fa> -q <query.psl>/<query.bed12> -g genome.fa [options]')
	parser.add_argument('collapse')
	required = parser.add_argument_group('required named arguments')
	required.add_argument('-r', '--reads', action='store', dest='r', \
		type=str, required=True, \
		help='FastA/FastQ files of raw reads')
	required.add_argument('-q', '--query', type=str, default='', required=True, \
		action='store', dest='q', help='PSL file of aligned/corrected reads')
	required.add_argument('-g', '--genome', action='store', dest='g', \
		type=str, required=True, help='fasta of reference genome')
	parser.add_argument('-m', '--minimap2', type=str, default='minimap2', \
		action='store', dest='m', help='path to minimap2')
	parser.add_argument('-t', '--threads', type=str, \
		action='store', dest='t', default='4', help='Minimap2 number of threads (4)')
	parser.add_argument('-w', '--window', default='20', \
		action='store', dest='w', help='Window size for comparing TSS/TES (20)')
	parser.add_argument('-s', '--support', default='3', \
		action='store', dest='s', help='Minimum number of supporting reads for an isoform (3)')
	parser.add_argument('-f', '--gtf', default='', \
		action='store', dest='f', help='GTF annotation file, used for renaming annotated isoforms')
	parser.add_argument('-e', '--report', default='default', \
		action='store', dest='e', help='Report options include\
		default: only full-length isoforms\
		comprehensive: default set + partial isoforms\
		ginormous: comprehensive + single exon isoforms')
	args = parser.parse_args()

	if args.m[-8:] != 'minimap2':
		if args.m[-1] == '/':
			args.m += 'minimap2'
		else:
			args.m += '/minimap2'

	sys.stderr.write('Collapsing isoforms\n')
	subprocess.call(['python', path+'bin/collapse_isoforms_precise.py', '-q', args.q, \
		'-w', args.w, '-s', '1', '-o', args.q[:-3]+'collapse1.psl'])
	sys.stderr.write('Filtering isoforms\n')
	subprocess.call(['python', path+'bin/filter_collapsed_isoforms.py', args.q[:-3]+'collapse1.psl', \
		args.e, args.q[:-3]+'collapse1.filtered.psl'])
	subprocess.call(['mv', args.q[:-3]+'collapse1.filtered.psl', args.q[:-3]+'collapse1.psl'])

	if args.f:
		sys.stderr.write('Renaming isoforms\n')
		subprocess.call(['python', path+'bin/identify_similar_annotated_isoform.py', \
			args.q[:-3]+'collapse1.psl', args.f, args.q[:-3]+'collapse1.renamed.psl'])
		subprocess.call(['mv', args.q[:-3]+'collapse1.renamed.psl', args.q[:-3]+'collapse1.psl'])

	subprocess.call(['python', path+'bin/psl_to_sequence.py', args.q[:-3]+'collapse1.psl', \
		args.g, args.q[:-3]+'collapse1.fa'])
	
	sys.stderr.write('Aligning reads to first-pass isoform reference\n')
	try:
		subprocess.call([args.m, '-a', '-t', args.t, '--secondary=no', \
			args.q[:-3]+'collapse1.fa', args.r], stdout=open(args.q[:-3]+'collapse1.sam', "w"))
	except:
		sys.stderr.write('Minimap2 error, specify minimap2 path with -m or make sure minimap2 is in your path\n')
		sys.exit()
		
	sys.stderr.write('Counting isoform expression\n')
	subprocess.call(['python', path+'bin/count_sam_genes.py', args.q[:-3]+'collapse1.sam', \
		args.q[:-3]+'collapse1.counts'])
	
	sys.stderr.write('Filtering isoforms by read coverage\n')
	subprocess.call(['python', path+'bin/match_counts.py', args.q[:-3]+'collapse1.counts', \
		args.q[:-3]+'collapse1.psl', args.s, args.q[:-3]+'isoforms.psl'])
	subprocess.call(['python', path+'bin/psl_to_sequence.py', args.q[:-3]+'isoforms.psl', \
		args.g, args.q[:-3]+'isoforms.fa'])
	if args.q[-3:] != 'psl':
		subprocess.call(['pslToBed', args.q[:-3]+'isoforms.psl', args.q[:-3]+'isoforms.bed'])
	
	sys.stderr.write('Removing intermediate files/done!\n')
	subprocess.call(['rm', args.q[:-3]+'collapse1.psl'])
	subprocess.call(['rm', args.q[:-3]+'collapse1.fa'])
	subprocess.call(['rm', args.q[:-3]+'collapse1.sam'])
	subprocess.call(['rm', args.q[:-3]+'collapse1.counts'])
