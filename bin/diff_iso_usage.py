import sys, csv
import scipy.stats as sps

try:
	counts_matrix = open(sys.argv[1])
	colname1 = sys.argv[2]
	colname2 = sys.argv[3]
	outfilename = sys.argv[4]
except:
	sys.stderr.write('usage: script.py counts_matrix colname1 colname2 diff_isos.txt\n')
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

counts = {}
for line in counts_matrix:
	line = line.rstrip().split('\t')
	iso_gene, count1, count2 = line[0], float(line[col1]), float(line[col2])
	if '_' not in iso_gene:
		sys.stderr.write('Please run bin/identify_annotated_gene.py first so that isoforms\
			can be grouped by their parent genes\n')
		sys.exit(1)
	iso = iso_gene[:iso_gene.rfind('_')]
	gene = iso_gene[iso_gene.rfind('_')+1:]
	if gene not in counts:
		counts[gene] = {}
	if count1 != 0 and counts2 != 0:
		counts[gene][iso] = [count1, count2] 

with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	geneordered = sorted(counts.keys())
	for gene in geneordered:
		generes = []
		for iso in counts[gene]:
			thesecounts = counts[gene][iso]
			othercounts = [0, 0]
			for iso_ in counts[gene]:  # count up for all other isoforms of this gene
				if iso_ == iso:
					continue
				othercounts[0] += counts[gene][iso_][0]
				othercounts[1] += counts[gene][iso_][1]
			ctable = [thesecounts, othercounts]
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

