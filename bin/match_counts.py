import sys, csv
counts = {}

try:
	countsfile = open(sys.argv[1])
	psl = open(sys.argv[2])
	min_reads = int(sys.argv[3])
	outfilename = sys.argv[4]  # sample_expression_pca_mtsf1.pdf
except:
	sys.stderr.write('usage: script.py countsfile psl_to_append_counts min_read_threshold outfilename\n')
	sys.exit()


for line in countsfile:
	line = line.rstrip().split('\t')
	counts[line[0]] = int(line[1])

with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	for line in psl:
		line = line.rstrip().split('\t')
		if line[9] in counts:
			count = counts[line[9]]
		else:
			count = 0
		if count >= min_reads:
			writer.writerow(line+[count])
