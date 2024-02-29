#!/usr/bin/env python3
import sys
import csv
import os
import argparse
import textwrap
from bed import Bed, BedBlock, intArraySplit

def main():
	parser = argparse.ArgumentParser(
	formatter_class=argparse.RawDescriptionHelpFormatter,
	description=textwrap.dedent('''\
Converts gtf to bed format.

Please note: GTF includes the stop codon  as part of the CDS, which
this method doesn't catch correctly. You cannot just move the end by 3 bases
because some genes have split stop codons.
We don't need this for Flair; if you want to have bed with included CDS, 
please download http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred

        '''))
	required = parser.add_argument_group('required named arguments')
	required.add_argument('gtf', type=str, help='annotated gtf')
	required.add_argument('bed', type=str, help='bed file')
	parser.add_argument('--include_gene', action='store_true', dest='include_gene', required=False, 
		help='''Include gene name in the isoform name''')
	args = parser.parse_args()

	gtf_to_bed(args.bed, args.gtf, args.include_gene)

class Cds(object):
	def __init__(self, start, end):
		self.start = start
		self.end = end
	def update(self, start, end):
		self.start = min(start, self.start)
		self.end = max(end, self.end)

class Gene(object):
	def __init__(self, gene, tx, chrom, start, end, strand):
		self.geneID = gene
		self.txID = tx
		self.chrom = chrom
		self.geneStart = start
		self.geneEnd = end
		self.strand = strand
		self.blockList = [BedBlock(start, end)]
	def add(self, start, end):
		self.geneStart = min(start, self.geneStart)
		self.geneEnd = max(end, self.geneEnd)
		self.blockList.append(BedBlock(start, end))
	def write(self, FH, include_gene, thickRegion=None):
		if include_gene:
			name = ('_'.join([self.txID, self.geneID]))
		else:
			name = self.txID
		blocks = sorted(self.blockList, key=lambda x: x.start)
		if thickRegion:
			bedObj = Bed(self.chrom, self.geneStart, self.geneEnd, name=name, strand=self.strand,
			blocks=blocks, numStdCols=12, thickStart=thickRegion.start, thickEnd=thickRegion.end)
		else:
			bedObj = Bed(self.chrom, self.geneStart, self.geneEnd, name=name, strand=self.strand, 
			blocks=blocks, numStdCols=12)
		bedObj.write(FH)


def ids_from_gtf(descrField):
	'''Parse the last field of a gtf to extract transcript and gene IDs'''

	pairs = descrField.split("; ")
	data_dict = {pair.split(" ", 1)[0].strip(): pair.split(" ", 1)[1].strip('"') for pair in pairs}
	return data_dict['gene_id'], data_dict['transcript_id']



def gtf_to_bed(outputfile, gtf, include_gene=False):
	iso_to_cds = {}
	geneObj = None
	with open(outputfile, 'wt') as outfile:
		writer = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep)

		prev_transcript, blockstarts, blocksizes, prev_gene, prev_chrom, prev_strand = [None, None, None, None, None, None]
		blockList = []
		# extract all exons from the gtf, keep exons grouped by transcript
		for line in open(gtf):  
			if line.startswith('#'):
				continue
			line = line.rstrip().split('\t')
			chrom, ty, start, end, strand = line[0], line[2], int(line[3]) - 1, int(line[4]), line[6]

			# please note: GTF includes the stop codon  as part of the CDS, which
			# this method doesn't catch correctly. You cannot just move the end by 3 bases
			# because some genes have split stop codons.
#			if ty == 'CDS':
#				if tx not in iso_to_cds:
#					iso_to_cds[tx] = Cds(start, end)
#				else:
#					iso_to_cds[tx].update(start, end)
			if ty != 'exon':
				continue
			gene, tx = ids_from_gtf(line[8])
			if geneObj is None:
				geneObj = Gene(gene, tx, chrom, start, end, strand)
			elif geneObj.txID == tx:
				geneObj.add(start, end)
			else:
				if geneObj.txID in iso_to_cds:
					geneObj.write(outfile, include_gene, thickRegion = iso_to_cds[geneObj.txID])
				else:
					geneObj.write(outfile, include_gene)
				geneObj = Gene(gene, tx, chrom, start, end, strand)

		# last entry...
		if geneObj.txID in iso_to_cds:
			geneObj.write(outfile, include_gene, thickRegion = iso_to_cds[geneObj.txID])
		else:
			geneObj.write(outfile, include_gene)


if __name__ == "__main__":
    main()
