#!/usr/bin/env python3

import sys, csv, os
import scipy.stats as sps

try:
	counts_matrix = open(sys.argv[1])
	colname1 = sys.argv[2]
	colname2 = sys.argv[3]
	outfilename = sys.argv[4]
except:
	sys.stderr.write('usage: diff_iso_usage.py counts_matrix colname1 colname2 diff_isos.txt\n')
	sys.exit()

header = counts_matrix.readline().rstrip().split('\t')

if colname1 in header:
	col1 = header.index(colname1)
else:
	sys.stderr.write('Could not find {} in {}\n'.format(colname1, ' '.join(header)))
	sys.exit(1)

if colname2 in header:
	col2 = header.index(colname2)
else:
	sys.stderr.write('Could not find {} in {}\n'.format(colname2, ' '.join(header)))
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
	else:
		iso = iso_gene[:iso_gene.rfind('_')]
		gene = iso_gene[iso_gene.rfind('_')+1:]
	return iso, gene

counts = {}
for line in counts_matrix:
	line = line.rstrip().split('\t')
	iso_gene, count1, count2 = line[0], float(line[col1]), float(line[col2])
	if '_' not in iso_gene:
		sys.stderr.write('Please run bin/identify_annotated_gene.py first so that isoforms\
			can be grouped by their parent genes\n')
		sys.exit(1)
	iso, gene = split_iso_gene(iso_gene)
	if gene not in counts:
		counts[gene] = {}
	counts[gene][iso] = [count1, count2] 

with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep)
	geneordered = sorted(counts.keys())
	for gene in geneordered:
		generes = []
		for iso in counts[gene]:
			thesecounts = counts[gene][iso]
			othercounts = [0, 0]
			ctable = [thesecounts, othercounts]
			if thesecounts[0] == 0 or thesecounts[1] == 0:  # do not test this isoform
				continue
			for iso_ in counts[gene]:  # count up for all other isoforms of this gene
				if iso_ == iso:
					continue
				othercounts[0] += counts[gene][iso_][0]
				othercounts[1] += counts[gene][iso_][1]
			ctable[1] = othercounts
			if ctable[0][0] + ctable[1][0] == 0 or ctable[0][1] + ctable[1][1] == 0 or not sum(ctable[1]):
				continue
			generes += [[gene, iso, sps.fisher_exact(ctable)[1]] + \
						 ctable[0] + ctable[1]]
		if not generes:
			writer.writerow([gene, iso, 'NA'] + ctable[0] + ctable[1])
			continue

		for res in generes:
			writer.writerow(res)

		# generes = sorted(generes, key=lambda x: x[2])
		# writer.writerow(generes[0])  # biggest differential for this gene

