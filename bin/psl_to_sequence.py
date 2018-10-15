import sys, csv

try:
	psl = open(sys.argv[1])
	genome = open(sys.argv[2])
	outfilename = sys.argv[3]
	if outfilename[-2:] == 'fq' or outfilename[-5:] == 'fastq':
		fastq = True
	else:
		fastq = False
except:
	sys.stderr.write('usage: script.py psl genome.fa outfilename\n')
	sys.exit(1)

psldata = {}
for line in psl:
	line = line.rstrip().split('\t')
	if line[13] in psldata:
		psldata[line[13]] += [line]
	else:
		psldata[line[13]] = [line]

def get_sequence(entry, seq):
	blocksizes = [int(x) for x in entry[18].split(',')[:-1]]
	blockstarts = [int(x) for x in entry[20].split(',')[:-1]]
	pulledseq = ''
	for block in range(len(blockstarts)):
		pulledseq += seq[blockstarts[block]:blockstarts[block]+blocksizes[block]]
	return pulledseq

with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	seq, chrom = '', ''
	for line in genome:
		line = line.rstrip()
		if line.startswith('>'):
			if not chrom:
				chrom = line.split()[0][1:]
				continue
			if chrom in psldata:
				for entry in psldata[chrom]:
					if fastq:
						writer.writerow(['@' + entry[9]])
					else:
						writer.writerow(['>' + entry[9]])
					pulledseq = get_sequence(entry, seq)
					writer.writerow([pulledseq])
					if fastq:
						writer.writerow(['+'])
						writer.writerow(['@'*len(pulledseq)])
			chrom = line.split()[0][1:]
			seq = ''
		else:
			seq += line
	if chrom in psldata:
		for entry in psldata[chrom]:
			if fastq:
				writer.writerow(['@'+entry[9]])
			else:
				writer.writerow(['>'+entry[9]])
			pulledseq = get_sequence(entry,seq)
			writer.writerow([pulledseq])
			if fastq:
				writer.writerow(['+'])
				writer.writerow(['@'*len(pulledseq)])

