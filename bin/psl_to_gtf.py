import sys, csv

try:
	psl = open(sys.argv[1])
except:
	sys.stderr.write('usage: script.py pslfile > outfile.gtf \n')
	sys.stderr.write('Entry name must contain underscore-delimited transcriptid and geneid like so:\
	 ENST00000318842.11_ENSG00000156313.12 or a4bab8a3-1d28_chr8:232000\n')
	sys.exit(1)

for line in psl:
	line = line.rstrip().split('\t')
	chrom, strand, score, name = line[13], line[8], line[0], line[9]
	tstarts = [int(n) for n in line[20].split(',')[:-1]]  # target starts
	bsizes = [int(n) for n in line[18].split(',')[:-1]]  # block sizes
	
	if '_' not in name:
		sys.stderr.write('Please first run bin/identify_annotated_gene.py or \
			bin/identify_gene_isoform.py prior for best results\n')
		sys.exit()

	if 'ENSG' in name:
		gene_id = name[name.find('ENSG'):]
	elif 'chr' in name:
		gene_id = name[name.find('chr'):]
	else:
		gene_id = name[name.find('_')+1:]
	if '-' in gene_id:
		gene_id = gene_id[:gene_id.find('-')]
	transcript_id = name[:name.rfind('_')]	

	for b in range(len(tstarts)):
		exon_assignment = transcript_id + '_' + str(b)
		endstring = 'gene_id \"{}\"; transcript_id \"{}\"; exon_assignment \"{}\";'\
				.format(gene_id, transcript_id, exon_assignment)
		print('\t'.join([chrom, 'FLAIR', 'exon', str(tstarts[b]), \
			str(tstarts[b]+bsizes[b]), '.', strand, str(score), endstring]))
