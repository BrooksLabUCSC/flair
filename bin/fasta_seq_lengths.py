import sys, csv

try:
	fasta = open(sys.argv[1])
	outfilename = sys.argv[2]
except:
	sys.stderr.write('usage: script.py fasta outfilename\n')
	sys.exit(1)

with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	seqlen = 0
	for line in fasta:
		line = line.rstrip()
		if line.startswith('>'):
			if seqlen:
				writer.writerow([name, seqlen])
			name = line[1:]
#			print(name)

			seqlen = 0
			continue
		seqlen += len(line.rstrip())
	writer.writerow([name, seqlen])
