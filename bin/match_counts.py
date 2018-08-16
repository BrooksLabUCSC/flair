import sys, csv
counts = {}
for line in open(sys.argv[1]):
	line = line.rstrip().split('\t')
	counts[line[0]] = int(line[1])

with open(sys.argv[4], 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	for line in open(sys.argv[2]):
		line = line.rstrip().split('\t')
		if line[9] in counts:
			count = counts[line[9]]
		else:
			count = 0
		if count >= int(sys.argv[3]):
			writer.writerow(line+[count])
