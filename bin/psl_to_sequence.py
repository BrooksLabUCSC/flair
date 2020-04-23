#!/usr/bin/env python3
import sys, csv, os

try:
	psl = open(sys.argv[1])
	isbed = sys.argv[1][-3:].lower() != 'psl'
	genome = open(sys.argv[2])
	outfilename = sys.argv[3]
	if outfilename[-2:].lower() == 'fq' or outfilename[-5:].lower() == 'fastq':
		fastq = True
	else:
		fastq = False
except:
	sys.stderr.write('usage: script.py psl|bed genome.fa outfilename\n')
	sys.exit(1)

psldata = {}
for line in psl:  # or bed
	line = line.rstrip().split('\t')
	chrom = line[0] if isbed else line[13]
	if chrom not in psldata:
		psldata[chrom] = []
	psldata[chrom] += [line]

def get_sequence(entry, seq):
	if isbed:
		start = int(entry[1])
		blockstarts = [int(n) + start for n in entry[11].split(',')[:-1]]
		blocksizes = [int(n) for n in entry[10].split(',')[:-1]]
		strand = entry[5]
	else:
		blocksizes = [int(n) for n in entry[18].split(',')[:-1]]
		blockstarts = [int(n) for n in entry[20].split(',')[:-1]]
		strand = entry[8]
	pulledseq = ''
	for block in range(len(blockstarts)):
		pulledseq += seq[blockstarts[block]:blockstarts[block]+blocksizes[block]]
	if strand == '-':
		pulledseq = revcomp(pulledseq)
	return pulledseq

def revcomp(seq):
	seq = seq.replace('A', 'X').replace('T', 'A').replace('X', 'T')
	seq = seq.replace('G', 'X').replace('C', 'G').replace('X', 'C')
	return seq[::-1]

with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep)
	seq, chrom = '', ''
	for line in genome:
		line = line.rstrip()
		if line.startswith('>'):
			if not chrom:
				chrom = line.split()[0][1:]
				continue
			if chrom in psldata:  # or bed
				for entry in psldata[chrom]:
					name = entry[3] if isbed else entry[9]
					if fastq:
						writer.writerow(['@' + name])
					else:
						writer.writerow(['>' + name])
					pulledseq = get_sequence(entry, seq)
					writer.writerow([pulledseq])
					if fastq:
						writer.writerow(['+'])
						writer.writerow(['@'*len(pulledseq)])
			chrom = line.split()[0][1:]
			seq = ''
		else:
			seq += line

	if chrom in psldata:  # last chromosome
		for entry in psldata[chrom]:
			name = entry[3] if isbed else entry[9]
			if fastq:
				writer.writerow(['@'+name])
			else:
				writer.writerow(['>'+name])
			pulledseq = get_sequence(entry,seq)
			writer.writerow([pulledseq])
			if fastq:
				writer.writerow(['+'])
				writer.writerow(['@'*len(pulledseq)])

