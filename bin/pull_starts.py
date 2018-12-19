import sys,csv
try:
    psl = open(sys.argv[1])
    outfilename = sys.argv[2]
except:
    sys.stderr.write('script.py psl outfilename\n')
    sys.exit(1)

with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	for line in psl:
		line = line.rstrip().split('\t')
		if line[20].count(',') == 1:
			writer.writerow([line[13], int(line[15]), int(line[15]), line[9]])
			writer.writerow([line[13], int(line[16]), int(line[16]), line[9]])
		if line[8] == '+':
			start = int(line[15])
		elif line[8] == '-':
			start = int(line[16])
		else:
			start = int(line[15])
			writer.writerow([line[13], start, start, line[9]])
			start = int(line[16])		
		writer.writerow([line[13], start, start, line[9]])