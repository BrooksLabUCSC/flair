#!/usr/bin/env python3
import sys, csv, os

try:
	isoforms = open(sys.argv[1])
	isbed = sys.argv[1][-3:].lower() != 'psl' 
	s = float(sys.argv[2])
	outfilename = sys.argv[3]
except:
	sys.stderr.write('usage: filter_isoforms_by_proportion_of_gene_expr.py isoforms support_percentage outfile\n')
	sys.exit(1)

if s >= 1:
	sys.stderr.write('Support percentage should be a decimal e.g. 0.1 for 10%\n')
	sys.exit(1)

def split_iso_gene(iso_gene):
    if '_chr' in iso_gene:
        splitchar = '_chr'
    elif '_XM' in iso_gene:
        splitchar = '_XM'
    elif '_XR' in iso_gene:
        splitchar = '_XR'
    elif '_NM' in iso_gene:
        splitchar = '_NM'
    elif '_NR' in iso_gene:
        splitchar = '_NR'
    elif '_R2_' in iso_gene:
        splitchar = '_R2_'
    elif '_NC_' in iso_gene:
        splitchar = '_NC_'
    else:
        splitchar = '_'
    iso = iso_gene[:iso_gene.rfind(splitchar)]
    gene = iso_gene[iso_gene.rfind(splitchar)+1:]
    return iso, gene

genes = {}
for line in isoforms:
	line = line.rstrip().split('\t')
	if isbed:
		name = line[3]
	else:
		name = line[9]
	iso, gene = split_iso_gene(name)
	if gene not in genes:
		genes[gene] = set()
	genes[gene].add(tuple(line))
	# print(tuple(line))

with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep)
	for gene in genes:
		gene_total = float(sum([float(iso[0][-1]) for iso in genes[gene]]))
		if gene_total == 0:
			continue
		for iso in genes[gene]:
			if float(iso[0][-1])/gene_total >= s:
				writer.writerow(iso)

