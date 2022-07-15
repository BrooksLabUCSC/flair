#!/usr/bin/env python3
import sys, csv, os, argparse

parser = argparse.ArgumentParser(description='''converts a gtf to a bed or psl, depending on the output filename extension;
	gtf exons need to be grouped by transcript and sorted by coordinate w/in a transcript''',
	usage='gtf_to_psl.py in.gtf out.psl|bed [options]')
required = parser.add_argument_group('required named arguments')
required.add_argument('gtf', action='store',
	type=str, help='annotated gtf')
required.add_argument('psl', action='store',
	type=str, help='psl or bed file')
parser.add_argument('--include_gene', action='store_true', dest='include_gene',
	required=False, help='''Include gene name in the isoform name''')
parser.add_argument('--chrom_sizes', action='store', dest='chrom_sizes', default='',
	required=False, help='''chrom sizes file for psl conversion, recommended''')
args = parser.parse_args()

isbed = args.psl[-3:].lower() != 'psl'

chrom_to_size = {}
if args.chrom_sizes:
	for line in open(args.chrom_sizes):
		line = line.rstrip().split('\t')
		chrom_to_size[line[0]] = line[1]

missing_chroms = set()
iso_to_cds = {}
with open(args.psl, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep)

	prev_transcript = ''
	for line in open(args.gtf):  # extract all exons from the gtf, keep exons grouped by transcript
		if line.startswith('#'):
			continue
		line = line.rstrip().split('\t')
		chrom, ty, start, end, strand = line[0], line[2], int(line[3]) - 1 , int(line[4]), line[6]
		this_transcript = line[8][line[8].find('transcript_id')+15:]
		this_transcript = this_transcript[:this_transcript.find('"')]

		if ty == 'CDS':
			if this_transcript not in iso_to_cds:
				iso_to_cds[this_transcript] = [start, end]
			elif end > iso_to_cds[this_transcript][1]:
				iso_to_cds[this_transcript][1] = end
		if ty != 'exon':
			continue

		# once all the exons for a transcript are read, write the psl/bed entry
		if this_transcript != prev_transcript:
			if prev_transcript:
				blockcount = len(blockstarts)
				if blockcount > 1 and blockstarts[0] > blockstarts[1]:  # need to reverse exons
					blocksizes = blocksizes[::-1]
					blockstarts = blockstarts[::-1]

				tstart, tend = blockstarts[0], blockstarts[-1] + blocksizes[-1]  # target (e.g. chrom)
				qsize = sum(blocksizes)  # query (e.g. transcript)
				if args.include_gene:
					qname = prev_transcript+'_'+prev_gene
				else:
					qname = prev_transcript
				if not isbed: # psl specific
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
					if qname in iso_to_cds:
						cds_start, cds_end = iso_to_cds[qname]
					else:
						cds_start, cds_end = tstart, tend
					writer.writerow([prev_chrom, tstart, tend, qname, 1000, prev_strand, cds_start,
						cds_end, 0, blockcount, blocksizes, relblockstarts])
				else:
					blockstarts = ','.join([str(b) for b in blockstarts]) + ','
					if chrom_to_size and prev_chrom in chrom_to_size:
						writer.writerow([0, 0, 0, 0, 0, 0, 0, 0, prev_strand, qname, qsize, 0, qsize,
							prev_chrom, chrom_to_size[prev_chrom], tstart, tend, blockcount, blocksizes, qstarts, blockstarts])
					else:
						if chrom_to_size and prev_chrom not in chrom_to_size:
							missing_chroms.add(prev_chrom)
						writer.writerow([0, 0, 0, 0, 0, 0, 0, 0, prev_strand, qname, qsize, 0, qsize,
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
	this_gene = line[8][line[8].find('gene_id')+9:]
	this_gene = this_gene[:this_gene.find('"')]
	if blockcount > 1 and blockstarts[0] > blockstarts[1]:  # need to reverse exons
		blocksizes = blocksizes[::-1]
		blockstarts = blockstarts[::-1]
	qsize = sum(blocksizes)  # query (e.g. transcript)
	tstart, tend = blockstarts[0], blockstarts[-1] + blocksizes[-1]  # target (e.g. chrom)
	blocksizes = ','.join([str(b) for b in blocksizes]) + ','
	if args.include_gene:
		qname = this_transcript+'_'+this_gene
	else:
		qname = this_transcript
	if isbed:
		relblockstarts = [block - tstart for block in blockstarts]		
		relblockstarts = ','.join([str(b) for b in relblockstarts]) + ','
		if qname in iso_to_cds:
			cds_start, cds_end = iso_to_cds[qname]
		else:
			cds_start, cds_end = tstart, tend
		writer.writerow([chrom, tstart, tend, qname, 1000, strand, cds_start, cds_end, 0,
			blockcount, blocksizes, relblockstarts])
	else:
		blockstarts = ','.join([str(b) for b in blockstarts]) + ','
		if chrom_to_size and prev_chrom in chrom_to_size:
			writer.writerow([0, 0, 0, 0, 0, 0, 0, 0, strand, qname, qsize, 0, qsize,
				chrom, chrom_to_size[chrom], tstart, tend, blockcount, blocksizes, qstarts, blockstarts])
		else:
			if chrom_to_size and chrom not in chrom_to_size:
				missing_chroms.add(chrom)
			writer.writerow([0, 0, 0, 0, 0, 0, 0, 0, strand, qname, qsize, 0, qsize,
				chrom, 0, tstart, tend, blockcount, blocksizes, qstarts, blockstarts])
if missing_chroms:
	sys.stderr.write('chromosomes found in gtf but not in chrom_sizes file: {}\n'.format(missing_chroms))

