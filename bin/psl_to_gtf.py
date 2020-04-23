#!/usr/bin/env python3
import sys, csv

try:
	psl = open(sys.argv[1])
	isbed = sys.argv[1][-3:].lower() != 'psl'
	force = len(sys.argv) > 2 and 'force' in sys.argv[2]
except:
	sys.stderr.write('usage: script.py file.psl|file.bed > outfile.gtf \n')
	sys.stderr.write('Entry name must contain underscore-delimited transcriptid and geneid like so:\
	 ENST00000318842.11_ENSG00000156313.12 or a4bab8a3-1d28_chr8:232000\n')
	sys.exit(1)

def split_iso_gene(iso_gene):
	if '_chr' in iso_gene:
		iso = iso_gene[:iso_gene.rfind('_chr')]
		gene = iso_gene[iso_gene.rfind('_chr')+1:]
	elif '_XM' in iso_gene:
		iso = iso_gene[:iso_gene.rfind('_XM')]
		gene = iso_gene[iso_gene.rfind('_XM')+1:]
	elif '_XR' in iso_gene:
		iso = iso_gene[:iso_gene.rfind('_XR')]
		gene = iso_gene[iso_gene.rfind('_XR')+1:]
	elif '_NM' in iso_gene:
		iso = iso_gene[:iso_gene.rfind('_NM')]
		gene = iso_gene[iso_gene.rfind('_NM')+1:]
	elif '_NR' in iso_gene:
		iso = iso_gene[:iso_gene.rfind('_NR')]
		gene = iso_gene[iso_gene.rfind('_NR')+1:]
	elif '_R2_' in iso_gene:
		iso = iso_gene[:iso_gene.rfind('_R2_')]
		gene = iso_gene[iso_gene.rfind('_R2_')+1:]		
	else:
		iso = iso_gene[:iso_gene.rfind('_')]
		gene = iso_gene[iso_gene.rfind('_')+1:]
	return iso, gene

for line in psl:
	line = line.rstrip().split('\t')
	if isbed:
		start = int(line[1])
		chrom, strand, score, name, start = line[0], line[5], line[4], line[3], int(line[1])
		tstarts = [int(n) + start for n in line[11].split(',')[:-1]]
		bsizes = [int(n) for n in line[10].split(',')[:-1]]
	else:
		chrom, strand, score, name, start = line[13], line[8], line[0], line[9], int(line[15])
		tstarts = [int(n) for n in line[20].split(',')[:-1]]  # target starts
		bsizes = [int(n) for n in line[18].split(',')[:-1]]  # block sizes
	
	if '_' not in name and not force:
		sys.stderr.write('No GTF conversion was done. Please run bin/identify_gene_isoform.py first\n')
		sys.stderr.write('for best results, or run with --force\n')
		sys.exit(1)

	if ';' in name:
		name = name.replace(';', ':')

	transcript_id, gene_id = split_iso_gene(name)

	endstring = 'gene_id \"{}\"; transcript_id \"{}\";'\
				.format(gene_id, transcript_id)
	print('\t'.join([chrom, 'FLAIR', 'transcript', str(start+1), str(tstarts[-1]+bsizes[-1]), '.', strand, '.', \
		endstring]))
	# if strand == '-':  # to list exons in 5'->3'
	# 	for b in range(len(tstarts)):  # exon number
	# 		bi = len(tstarts) - 1 - b  # block index
	# 		endstring = 'gene_id \"{}\"; transcript_id \"{}\"; exon_number \"{}\";'\
	# 						.format(gene_id, transcript_id, b)
	# 		print('\t'.join([chrom, 'FLAIR', 'exon', str(tstarts[bi]+1), \
	# 			str(tstarts[bi]+bsizes[bi]), '.', strand, str(score), endstring]))			
	# else:
	for b in range(len(tstarts)):
		endstring = 'gene_id \"{}\"; transcript_id \"{}\"; exon_number \"{}\";'\
				.format(gene_id, transcript_id, b)
		print('\t'.join([chrom, 'FLAIR', 'exon', str(tstarts[b]+1), \
			str(tstarts[b]+bsizes[b]), '.', strand, str(score), endstring]))
