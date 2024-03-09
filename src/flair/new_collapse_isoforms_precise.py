#!/usr/bin/env python3

import pyranges as pr
import pandas as pd
import sys
import argparse
from collections import Counter
from bed import Bed, BedBlock

def main():
	parser = argparse.ArgumentParser(description='collapse parse options',
				usage='python collapse_isoforms_precise.py -q <query.psl>/<query.bed> [options]')
	required = parser.add_argument_group('required named arguments')
	required.add_argument('-q', '--query', type=str, required=True, 
		help='BED12 or PSL file of aligned/corrected reads. PSL should end in .psl')
	parser.add_argument('-o', '--output', type=str, 
		help='specify output file, should agree with query file type')
	parser.add_argument('-w', '--window', default=200, type=int,
		 help='window size for grouping TSS/TES (200)')
	parser.add_argument('-s', '--support', default=0.25, type=float, 
		help='minimum proportion(s<1) or number of supporting reads(s>=1) for an isoform (0.25)')
	parser.add_argument('-f', '--gtf', type=str,
		 help='GTF annotation file for selecting annotated TSS/TES')
	parser.add_argument('-m', '--max_results', default=2, type=int, 
		help='maximum number of novel TSS or TES picked per isoform unless --no_redundant is specified (2)')
	parser.add_argument('-t', '--threads', default=2, type=int,
		 help='number of threads to use (2)')
	parser.add_argument('-n', '--no_redundant', default='best_only',
		help='For each unique splice junction chain, report options include: \
		none: multiple best TSSs/TESs chosen for each unique set of splice junctions, see M; \
		longest: TSS/TES chosen to maximize length; \
		best_only: single best TSS/TES used in conjunction chosen; \
		longest/best_only override max_results argument immediately before output \
		resulting in one isoform per unique set of splice junctions (default: best_only)')
	parser.add_argument('-c', '--clean', default=False, action='store_true',
		help='Specify this to not append read support to the end of each entry (default: not specified)')
	parser.add_argument('-i', '--isoformtss', default=False, action='store_true',
		help='when specified, TSS/TES for each isoform will be determined from supporting reads \
		for individual isoforms (default: not specified, determined at the gene level)')
	parser.add_argument('--nosplice', default='chrM',
		help='Comma separated list of chromosomes that should not have spliced isoforms (default: chrM)')
	
	args = parser.parse_args()
	
	if args.output:
		if args.output[-3:].lower() != args.query[-3:].lower():
			sys.stderr.write('Make sure input and output file extensions agree\n')
			sys.exit(1)
	else:  # default output name
		args.output = args.query[:-3]+'collapsed.bed'
	
	# This renaming of arguments is in preparation of turning this program into a function
	max_results=args.max_results
	threads=args.threads
	isoformtss=args.isoformtss
	gtfname=args.gtf
	clean=args.clean
	nosplice=args.nosplice
	collapse_isoforms_precise(queryfile=args.query, outfile=args.output, window=args.window, minsupport=args.support, 
		no_redundant=args.no_redundant)

def collapse_isoforms_precise(queryfile, outfile, window=200, minsupport=0.25, no_redundant='best_only'):
	'''Collapses input reads with identical intron patterns into one (best or longest) isoform'''
	# read file into pyranges object
	gr = pr.read_bed(queryfile)
	outFH = open(outfile, 'w')
	
	# this is used PER ISOFORM, so 'none' allows multiple isoforms with different starts and ends
	# readsupport only works when it's an integer (in original)
	
	# # single exon genes
	# singleExonGenes = gr[gr.BlockCount == 1]
	# if len(singleExonGenes) > 0:
	# # group overlapping genes on the same strand
	# 	clustered_ranges = singleExonGenes.cluster(strand=True, by=['Chromosome', 'Strand'])
	# 	df = clustered_ranges.df
	# 	# merge overlapping genes with starts and ends within the window and print
	# 	df.groupby(['Cluster']).apply(split_clusters, window)
	
	# for multiexon overlaps we want to know if the intron pattern agrees
	# This can easily be done by bed.py's getGaps
	multiExonGenes = gr[gr.BlockCount > 1]
	if len(multiExonGenes) == 0:
		sys.exit()
	clustered_ranges = multiExonGenes.cluster(strand=True, by=['Chromosome', 'Strand'])
	df = clustered_ranges.df
	intronsAndObject_df = df.apply(RowToBedObjs, axis=1)
	df = pd.concat([df, intronsAndObject_df], axis=1)
	for groupid, mdf in df.groupby('Cluster'):
		# gene is the combination of all isoforms
		gene = geneObj(mdf)
		# Group isoforms with identical introns and print if they pass minsupport criterium
		mdf.groupby(['Introns']).apply(printBestIsoform, outfilehandle=outFH, window=window, minsupport=minsupport, method=no_redundant, gene=gene)
		

	


def RowToBedObjs(row):
	'''Read a bed df row into a bed object and its introns in genome coordinates'''
	bedObj = Bed.parse(row[:12])
	introns = bedObj.getGaps()
	firstIntronStart = introns[0].start
	lastIntronEnd = introns[-1].end
	introns = bedObj.istring
	return pd.Series([introns, firstIntronStart, lastIntronEnd, bedObj], 
		index=['Introns', 'firstIntronStart', 'lastIntronEnd', 'BedObject'])

def update_bedObj(bedObj, newstart, newend):
	'''Extend chrom start and end and blocks'''
	bedBlocks = bedObj.blocks
	bedBlocks[0] = BedBlock(start=newstart, end=bedBlocks[0].end)
	bedBlocks[-1] = BedBlock(start=bedBlocks[-1].start, end=newend)
	bedObj.blocks = bedBlocks
	bedObj.chromStart = newstart
	bedObj.chromEnd = newend
	return bedObj

class geneObj(object):
	'''Information on the whole gene. Please note that this assumes that all overlapping
	multiexon isoforms that overlap belong to the same gene.'''
	def __init__(self, df):
		# get all isoforms that start at the gene start
		self.geneStart = df['Start'].min()
		rows_with_geneStart = df[df['Start'] == self.geneStart]
		# and get the intron boundary of their first exon
		self.firstExonEnds = rows_with_geneStart['firstIntronStart'].unique()
		# do the same for gene ends
		self.geneEnd = df['End'].max()
		rows_with_geneEnd = df[df['End'] == self.geneEnd]
		self.lastExonStarts = rows_with_geneEnd['lastIntronEnd'].unique()
		self.support = df.shape[0]

def printIfEnoughSupport(bedObj, outfilehandle, geneReadCount, minsupport=0.25):
	'''Prints bedObject if it passes minsupport requirements. Number of reads supporting
	this isoform is in the score field'''
	if minsupport < 1:
		if bedObj.score/geneReadCount < minsupport:
			return False
	# minimum required number of reads
	elif bedObj.score < minsupport:
		return False
	bedObj.write(outfilehandle, noCDS=True)

def printBestIsoform(df, outfilehandle, window=100, minsupport=0.25, method='best_only', gene=False):
	'''Input is a df with isoforms. Determines most occurring starts and ends, and prints that version of the isoform'''
	#method='best_only'  # I don't know how to do this without making a hybrid read.
	#method='longest'   # should this be the single longest read, or a combination?
	bedObj = df['BedObject'].iloc[0]
	# if there's only one read supporting this isoform, just print it
	if df.shape[0] == 1:
		bedObj.score = 1 # 1 read
		printIfEnoughSupport(bedObj, outfilehandle=outfilehandle, geneReadCount=gene.support, 
			minsupport=minsupport)
		return True
	# if there are multiple reads, determine the best start and end
	if method == 'best_only':   # must supported (and after that longest) isoform
		groupStart = int(df['Start'].mode().min())
		groupEnd = int(df['End'].mode().max())
	elif  method == 'longest':
		groupStart = int(df['Start'].min())
		groupEnd = int(df['End'].max())
	else:
		return False   # TODO: maybe split as before (for 'none')
	bedObj = update_bedObj(bedObj, groupStart, groupEnd)
	bedObj.score = df.shape[0]
	printIfEnoughSupport(bedObj, outfilehandle=outfilehandle, geneReadCount=gene.support, minsupport=minsupport)


if __name__ == "__main__":
    main()



