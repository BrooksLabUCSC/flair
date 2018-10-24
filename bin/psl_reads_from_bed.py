import sys, csv
keepnames = set()
for line in open(sys.argv[1]):  # bed of promoter-supported reads
	line = line.rstrip().split('\t')
	keepnames.add(line[3])

with open(sys.argv[3], 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	for line in open(sys.argv[2]):  # psl
		line = line.rstrip().split('\t')
		if line[9] in keepnames:
			writer.writerow(line)
