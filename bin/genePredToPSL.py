import sys, csv

try:
	psl = open(sys.argv[1])
	gp = open(sys.argv[2])
	outfilename = sys.argv[3]
	chop = len(sys.argv) > 4
except:
	sys.stderr.write('usage: script.py psl_in corrected.gp pslout\n')
	sys.exit(1)

pslentries = {}
for line in psl:
	line = line.rstrip().split('\t')
	if line[9] in pslentries:
		pslentries[line[9]] += [line]  # multiple alignments for the same entry
	else:
		pslentries[line[9]] = [line]

with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	for line in gp:
		name, chrom, strand, txStart, txEnd, \
			cdsStart, cdsEnd, exonCount, exonStarts, \
			exonEnds = line.rstrip().split('\t')
		if name[:-2] in pslentries:
			name = name[:-2]
		pslentry = pslentries[name]
		pslentry = sorted(pslentry, key=lambda x: x[0])
		exonStarts = exonStarts.split(',')[:-1]
		exonEnds = exonEnds.split(',')[:-1]
		blockSizes = ','.join([str(int(e) - int(s)) for (e, s) in zip(exonEnds, exonStarts)])+','
		exonStartsRel = ','.join([str(int(s) - int(txStart)) for s in exonStarts])+','
		if '-' in blockSizes + exonStartsRel:
			print(line.rstrip())
			continue
		exonStartsAbs = ','.join(exonStarts)+','
		i = 0
		while i < len(pslentry) and pslentry[i][13] != chrom:
			i += 1
		if pslentry[i][13] != chrom:
			sys.exit()
		pslentry[i][15] = txStart
		pslentry[i][16] = txEnd
		row = pslentry[i][:17] + [exonCount, blockSizes, exonStartsRel, exonStartsAbs]
		writer.writerow(pslentry[i][:17] + [exonCount, blockSizes, exonStartsRel, exonStartsAbs])