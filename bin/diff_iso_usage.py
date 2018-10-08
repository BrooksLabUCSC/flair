import sys, csv
import scipy.stats as sps

try:
	psl = open(sys.argv[1])
	col = int(sys.argv[2])  # 25
	outfilename = sys.argv[3]
except:
	sys.stderr.write('usage: script.py isoforms.psl colnum diff_isos.txt\n')
	sys.exit()

counts = {}

for line in psl:
	line = line.rstrip().split('\t')
	if '_' not in line[9]:
		sys.stderr.write('Please run bin/infer_strand_for_psl.py first,\
			necessary for grouping isoforms into genes.\n')
		sys.exit()
	gene = line[9][line[9].rfind('_')+1:]
	if '-' in gene:
		gene = gene[:gene.find('-')]
	if gene.count('.') == 2:
		gene = gene[:gene.rfind('.')]
	if gene not in counts:
		counts[gene] = {}
	if line[col] == 'NA':
		line[col] = 0
	if line[col+1] == 'NA':
		line[col+1] = 0

	counts[gene][line[9]] = [float(line[col]), float(line[col+1])]


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
			continue
		generes = sorted(generes, key=lambda x: x[2])
		writer.writerow(generes[0])  # biggest differential for this gene

