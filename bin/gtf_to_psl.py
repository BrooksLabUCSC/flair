#!/usr/bin/env python3
import sys, csv, os

try:
	gtf = open(sys.argv[1])
	outfilename = sys.argv[2]
	isbed = sys.argv[2][-3:].lower() != 'psl'
	if len(sys.argv) > 3:
		chrom_sizes_file = open(sys.argv[3])
	else:
		chrom_sizes_file = ''
except:
	sys.stderr.write('usage: script.py in.gtf out.psl|bed \n')
	sys.stderr.write('converts a gtf to a bed or psl, depending on the output filename extension\n')
	sys.stderr.write('a chrom sizes file is recommended for psl conversion:\n')
	sys.stderr.write('usage: script.py in.gtf out.psl chrom.sizes \n')

	sys.exit(1)

chrom_to_size = {}
for line in chrom_sizes_file:
	line = line.rstrip().split('\t')
	chrom_to_size[line[0]] = line[1]

missing_chroms = set()

with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep)

	prev_transcript = ''
	for line in gtf:  # extract all exons from the gtf, keep exons grouped by transcript
		if line.startswith('#'):
			continue
		line = line.rstrip().split('\t')
		chrom, ty, start, end, strand = line[0], line[2], int(line[3]) - 1 , int(line[4]), line[6]
		if ty != 'exon':
			continue

		this_transcript = line[8][line[8].find('transcript_id')+15:]
		this_transcript = this_transcript[:this_transcript.find('"')]

		# once all the exons for a transcript are read, write the psl/bed entry
		if this_transcript != prev_transcript:
			if prev_transcript:
				blockcount = len(blockstarts)
				if blockcount > 1 and blockstarts[0] > blockstarts[1]:  # need to reverse exons
					blocksizes = blocksizes[::-1]
					blockstarts = blockstarts[::-1]

				tstart, tend = blockstarts[0], blockstarts[-1] + blocksizes[-1]  # target (e.g. chrom)
				qsize = sum(blocksizes)  # query (e.g. transcript)
				qname = prev_transcript#+'_'+this_gene
				if not isbed:  # psl specific
					pos = 0
					qstarts = [pos]
					for b in blocksizes[:-1]:
						pos += b
						qstarts += [pos]
					qstarts = ','.join([str(b) for b in qstarts]) + ',' 

				blocksizes = ','.join([str(b) for b in blocksizes]) + ','

				if isbed:
					relblockstarts = [block - tstart for block in blockstarts]
					relblockstarts = ','.join([str(b) for b in relblockstarts]) + ','
					writer.writerow([prev_chrom, tstart, tend, qname, 1000, prev_strand, tstart, \
						tend, 0, blockcount, blocksizes, relblockstarts])
				else:
					blockstarts = ','.join([str(b) for b in blockstarts]) + ','
					if chrom_to_size and prev_chrom in chrom_to_size:
						writer.writerow([0, 0, 0, 0, 0, 0, 0, 0, prev_strand, qname, qsize, 0, qsize, \
							prev_chrom, chrom_to_size[prev_chrom], tstart, tend, blockcount, blocksizes, qstarts, blockstarts])
					else:
						if chrom_to_size and prev_chrom not in chrom_to_size:
							missing_chroms.add(prev_chrom)
						writer.writerow([0, 0, 0, 0, 0, 0, 0, 0, prev_strand, qname, qsize, 0, qsize, \
							prev_chrom, 0, tstart, tend, blockcount, blocksizes, qstarts, blockstarts])

			blockstarts, blocksizes = [], []
			prev_transcript = this_transcript
			prev_gene = line[8][line[8].find('gene_id')+9:]
			prev_gene = prev_gene[:prev_gene.find('"')]
			prev_chrom = chrom
			prev_strand = strand

		blockstarts += [start]
		blocksizes += [end-start]

	# last entry...
	# this_gene = line[8][line[8].find('gene_id')+9:]
	# this_gene = this_gene[:this_gene.find('"')]
	if blockcount > 1 and blockstarts[0] > blockstarts[1]:  # need to reverse exons
		blocksizes = blocksizes[::-1]
		blockstarts = blockstarts[::-1]
	qsize = sum(blocksizes)  # query (e.g. transcript)
	tstart, tend = blockstarts[0], blockstarts[-1] + blocksizes[-1]  # target (e.g. chrom)
	blocksizes = ','.join([str(b) for b in blocksizes]) + ','
	qname = this_transcript#+'_'+this_gene
	if isbed:
		relblockstarts = [block - tstart for block in blockstarts]		
		relblockstarts = ','.join([str(b) for b in relblockstarts]) + ','
		writer.writerow([chrom, tstart, tend, qname, 1000, strand, start, end, 0, \
			blockcount, blocksizes, relblockstarts])
	else:
		blockstarts = ','.join([str(b) for b in blockstarts]) + ','
		if chrom_to_size and prev_chrom in chrom_to_size:
			writer.writerow([0, 0, 0, 0, 0, 0, 0, 0, strand, qname, qsize, 0, qsize, \
				chrom, chrom_to_size[chrom], tstart, tend, blockcount, blocksizes, qstarts, blockstarts])
		else:
			if chrom_to_size and chrom not in chrom_to_size:
				missing_chroms.add(chrom)
			writer.writerow([0, 0, 0, 0, 0, 0, 0, 0, strand, qname, qsize, 0, qsize, \
				chrom, 0, tstart, tend, blockcount, blocksizes, qstarts, blockstarts])
if missing_chroms:
	sys.stderr.write('chromosomes found in gtf but not in chrom_sizes file: {}\n'.format(missing_chroms))

