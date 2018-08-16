""" AD Tang 2018 """

import sys, argparse, subprocess

if len(sys.argv) > 1 and sys.argv[1] == 'align':
	mode = 'align'
elif len(sys.argv) > 1 and sys.argv[1] == 'correct':
	mode = 'correct'
elif len(sys.argv) > 1 and sys.argv[1] == 'collapse':
	mode = 'collapse'
else:
	sys.stderr.write('usage: python flair.py mode --help \n')
	sys.stderr.write('modes: correct, collapse\n')
	sys.exit()

path = sys.argv[0][:sys.argv[0].rfind('/')+1] if '/' in sys.argv[0] else ''
if mode == 'align':
	parser = argparse.ArgumentParser(description='flair-align parse options', \
				usage='python flair.py align -r <reads.fq>/<reads.fa> -g genome.fa [options]')
	parser.add_argument('align')
	required = parser.add_argument_group('required named arguments')
	required.add_argument('-r', '--reads', action='store', dest='a', \
		type=str, required=True, \
		help='FastA/FastQ files of raw reads')
	required.add_argument('-g', '--genome', action='store', dest='g', \
		type=str, required=True, help='FastA of reference genome')
	parser.add_argument('-m', '--minimap2', type=str, default='', \
		action='store', dest='m', help='Path to minimap2')
	parser.add_argument('-o', '--output', \
		action='store', dest='o', default='', help='Output file name')
	args = parser.parse_args()
	
	args.o = args.o if args.o else args.r[:-3]+'sam'
	if args.m[-8:] != 'minimap2':
		if args.m[-1] == '/':
			args.m += 'minimap2'
		else:
			args.m += '/minimap2'

	mm_out = subprocess.call([args.m, '-a', '-t', args.t, '--secondary=no', args.g, args.r], stdout=subprocess.PIPE)	
	minimap2_out = open(args.o, "w")
	minimap2_out.write(mm_out.stdout.decode('utf-8'))
	minimap2_out.close()
	subprocess.call(['python', path+'sam_to_psl.py', args.o, args.o[:-3]+'psl'])

elif mode == 'correct':
	parser = argparse.ArgumentParser(description='flair-correct parse options', \
				usage='python flair.py correct -a annotated.gp -g genome.fa -q query.psl [options]')
	parser.add_argument('correct')
	required = parser.add_argument_group('required named arguments')
	required.add_argument('-a', '--annotation', action='store', dest='a', \
		type=str, required=True, \
		help='Genepred format splice junctions from annotation')
	parser.add_argument('-s', '--shortread', action='store', dest='s', \
		type=str, default='', \
		help='Genepred format splice junctions from short-read sequencing')
	required.add_argument('-g', '--genome', action='store', dest='g', \
		type=str, required=True, help='FastA of reference genome')
	required.add_argument('-q', '--query', type=str, default='', required=True, \
		action='store', dest='q', help='Uncorrected PSL file')
	parser.add_argument('-w', '--window', \
		action='store', dest='w', default=10, \
		help='Window size for correcting splice sites (10)')
	parser.add_argument('-m', '--mergesize', \
		action='store', dest='m', default=30, help='Merge size for gaps within exons (30)')
	parser.add_argument('-o', '--output', \
		action='store', dest='o', default='', help='Output file name')
	args = parser.parse_args()
	args.s = args.s if args.s else args.a
	stranded_psl = args.q[:args.q.rfind('.')] + '_strand.psl'
	args.o = args.o if args.o else 'correction_output'

	sys.stderr.write('Inferring strand for reads in PSL\n')
	subprocess.call(['python', path+'infer_strand_for_psl.py', args.q, args.a, stranded_psl])

	sys.stderr.write('Correcting reads in {}\n'.format(stranded_psl))
	subprocess.call(['python', path+'correctSplice.py', '-a', args.a, '-j', args.j, \
	'-g', args.g, '-q', stranded_psl, '-w', args.w, '-m', args.m, '-o', args.o, '-mode', 'closest'])

	sys.stderr.write('Converting output gp to PSL\n')
	subprocess.call(['python', path+'genePredToPSL.py', stranded_psl, args.o+'/corrected.gp',\
		args.o+'/'+stranded_psl[:-3]+'_corrected.psl'])

elif mode == 'collapse':
	parser = argparse.ArgumentParser(description='flair-collapse parse options', \
				usage='python flair.py correct -r <reads.fq>/<reads.fa> -q query.psl -g genome.fa [options]')
	parser.add_argument('collapse')
	required = parser.add_argument_group('required named arguments')
	required.add_argument('-r', '--reads', action='store', dest='a', \
		type=str, required=True, \
		help='FastA/FastQ files of raw reads')
	required.add_argument('-q', '--query', type=str, default='', required=True, \
		action='store', dest='q', help='PSL file of aligned/corrected reads')
	required.add_argument('-g', '--genome', action='store', dest='g', \
		type=str, required=True, help='FastA of reference genome')
	parser.add_argument('-m', '--minimap2', type=str, default='', \
		action='store', dest='m', help='Path to minimap2')
	parser.add_argument('-t', '--threads', type=str, \
		action='store', dest='t', default=12, help='Minimap2 number of threads (12)')
	parser.add_argument('-w', '--window', default=20, \
		action='store', dest='w', help='Window size for comparing TSS/TES (20)')
	parser.add_argument('-s', '--support', default=3, \
		action='store', dest='s', help='Minimum number of supporting reads for an isoform (3)')
	args = parser.parse_args()

	if args.m[-8:] != 'minimap2':
		if args.m[-1] == '/':
			args.m += 'minimap2'
		else:
			args.m += '/minimap2'

	sys.stderr.write('Removing isoforms containing unsupported junctions\n')
	# subprocess.call(['python', path+'remove_novels.py', sys.argv[2], '20', '1', sys.argv[3]])

	sys.stderr.write('Collapsing isoforms\n')
	subprocess.call(['python', path+'collapse_isoforms_precise.py', args.q, args.w, '1', args.q[:-3]+'collapse1.psl'])
	subprocess.call(['python', path+'psl_to_sequence.py', args.q[:-3]+'collapse1.psl', args.g, args.q[:-3]+'collapse1.fa'])
	
	sys.stderr.write('Aligning reads to collapsed reference\n')
	mm_out = subprocess.call([args.m, '-a', '-t', args.t, '--secondary=no', args.q[:-3]+'collapse1.fa', args.r], stdout=subprocess.PIPE)	
	minimap2_out = open(args.q[:-3]+'collapse1.sam', "w")
	minimap2_out.write(mm_out.stdout.decode('utf-8'))
	minimap2_out.close()
	subprocess.call(['rm', args.q[:-3]+'collapse1.fa'])

	# sys.stderr.write('Keeping only primary alignments\n')
	# print('samtools', 'view', '-h', '-F', '256', '-S', sys.argv[3][:-3]+'sam')
	# st_out = subprocess.call(['samtools', 'view', '-h', '-F', '256', '-S', sys.argv[3][:-3]+'sam'], stdout=subprocess.PIPE)
	# samtools_out = open(sys.argv[3][:-3]+'p.sam', "w")
	# samtools_out.write(st_out.stdout.decode('utf-8'))
	# samtools_out.close()

	subprocess.call(['python', path+'count_sam_genes.py', args.q[:-3]+'collapse1.sam', args.q[:-3]+'collapse1.counts'])
	sys.stderr.write('Removing alignment file\n')
	subprocess.call(['rm', args.q[:-3]+'collapse1.sam'])
	sys.stderr.write('Filtering isoforms by read coverage\n')
	subprocess.call(['python', path+'match_counts.py', args.q[:-3]+'collapse1.counts', args.q[:-3]+'collapse1.psl', args.s, args.q[:-3]+'final_isoforms.psl'])
	subprocess.call(['rm', args.q[:-3]+'collapse1.counts'])
