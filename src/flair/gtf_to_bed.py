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

        '''))
	required = parser.add_argument_group('required named arguments')
	required.add_argument('gtf', type=str, help='annotated gtf')
	required.add_argument('bed', type=str, help='bed file')
	parser.add_argument('--include_gene', action='store_true', dest='include_gene', required=False, 
		help='''Include gene name in the isoform name''')
	args = parser.parse_args()

	gtf_to_bed(args.bed, args.gtf, args.include_gene)

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
	def write(self, FH, include_gene, cdsObject):
		if include_gene:
			name = ('_'.join([self.txID, self.geneID]))
		else:
			name = self.txID
		blocks = sorted(self.blockList, key=lambda x: x.start)
		if cdsObject is None:
			bedObj = Bed(self.chrom, self.geneStart, self.geneEnd, name=name, strand=self.strand, 
			blocks=blocks, numStdCols=12)
		else:
			cdsObject.finalize()
			bedObj = Bed(self.chrom, self.geneStart, self.geneEnd, name=name, strand=self.strand,
			blocks=blocks, numStdCols=12, thickStart=cdsObject.thickStart, thickEnd=cdsObject.thickEnd)
		bedObj.write(FH)


def ids_from_gtf(descrField):
	'''Parse the last field of a gtf to extract transcript and gene IDs'''
	pairs = descrField.split("; ")
	data_dict = {pair.split(" ", 1)[0].strip(): pair.split(" ", 1)[1].strip('"') for pair in pairs}
	return data_dict['gene_id'], data_dict['transcript_id']

class Cds(object):
	'''We keep track of CDS because GTF puts the stop codon outside the ORF and bed does not.
	Stop codons may be split, so we cannot just add 3'''
	def __init__(self, name, strand):
		self.name = name
		self.strand = strand
	def add(self, ty, start, end):
		if ty == 'stop_codon':
			# stop codons may be split, and gtf does exclude them from the 
			# ORF while bed does not
			self.cdsstop = end if self.strand == '+' else start
		elif ty == 'CDS':
			# for the main cds body the strand doesn't matter
			# but the gtf exons may be in reverse
			self.start = min(start, self.start) if hasattr(self, 'start') else start
			self.end = max(end, self.end) if hasattr(self, 'end') else end
	def finalize(self):
		'''Not all genes have annotated start and stop codons even if they have CDS'''
		if self.strand == '+':
			self.thickStart, self.thickEnd = self.start, getattr(self, 'cdsstop', self.end)
		elif self.strand == '-':
			self.thickStart, self.thickEnd = getattr(self, 'cdsstop', self.start), self.end

def gtf_to_bed(outputfile, gtf, include_gene=False):
	iso_to_cds = {}
	geneObj = None
	cdsObj = None
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

			if ty in ['gene', 'transcript', 'UTR', 'start_codon']:
				continue

			gene, tx = ids_from_gtf(line[8])
			if geneObj is None:
				geneObj = Gene(gene, tx, chrom, start, end, strand)
				if ty in ['CDS', 'stop_codon']: 
					if cdsObj is None:
						cdsObj = Cds(tx, strand)
					cdsObj.add(ty, start, end)
			elif geneObj.txID == tx:
				if ty in ['CDS', 'stop_codon']: 
					if cdsObj is None:
						cdsObj = Cds(tx, strand)
					cdsObj.add(ty, start, end)
				elif ty == 'exon':
					geneObj.add(start, end)
			# this line is a new gene, so print the previous
			else:
				geneObj.write(outfile, include_gene, cdsObj)
				geneObj = Gene(gene, tx, chrom, start, end, strand)
				cdsObj = None

		# last entry...
		geneObj.write(outfile, include_gene, cdsObj)


if __name__ == "__main__":
    main()
