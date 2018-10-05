#!/usr/bin/env python3


########################################################################
# File: ssCorrect.py
#  executable: ssCorrect.py
# Purpose: 
#
#          
# Author: Cameron M. Soulette
# History:      cms 05/01/2018 Created
#
########################################################################


########################################################################
# Hot Imports & Global Variable
########################################################################


import os, sys
import numpy as np
from multiprocessing import Pool
from intervaltree import Interval, IntervalTree
import random
from tqdm import *

########################################################################
# CommandLine
########################################################################

class CommandLine(object) :
	'''
	Handle the command line, usage and help requests.
	CommandLine uses argparse, now standard in 2.7 and beyond. 
	it implements a standard command line argument parser with various argument options,
	and a standard usage and help,
	attributes:
	myCommandLine.args is a dictionary which includes each of the available command line arguments as
	myCommandLine.args['option'] 
	
	methods:
	
	'''
	
	def __init__(self, inOpts=None) :
		'''
		CommandLine constructor.
		Implements a parser to interpret the command line argv string using argparse.
		'''
		import argparse
		self.parser = argparse.ArgumentParser(description = ' ssCorrect.py - a tool to leverage annotation and short read data to correct misaligned splice junctions in short read data.',
											 epilog = 'Please feel free to forward any questions/concerns to /dev/null', 
											 add_help = True, #default is True 
											 prefix_chars = '-', 
											 usage = '%(prog)s -b reads.bed -g annotations.gtf -j other_junctions.bed')
		# Add args
		self.parser.add_argument('-b', "--bedFile", action = 'store', required=True, help='Input reads in bed12 format.')
		self.parser.add_argument('-g', "--gtf", action = 'store', required=False, help='Gencode annotation file.')
		self.parser.add_argument('-j', "--junctionsBed", action = 'store', required=False, help='Junction bed file.')
		self.parser.add_argument('-w', '--wiggleWindow', action = 'store', required=False, default = 15, help='Splice site correction window flank size.')
		self.parser.add_argument('--report_junctions', action = 'store', required=False, default= "norm", choices=['strict', 'norm', 'all'],
																						help = """Choose which types of nanopore splice sites to report:
																								"strict" - SS must be in gtf or junctionsBed, and have a single closest hit.
																								"norm" - SS must be in gtf or junctionsBed (tie for closest hit is random).
																								"all" - SS can be in gtf, junctionsBed, or unique to nanopore data.""")
		self.parser.add_argument('-p', "--threads", action = 'store', required=False, default = 2, help='Num threads.')
		self.parser.add_argument("--quiet", action = 'store_false', required=False, default = True, help='Do not display progress')
		if inOpts is None :
			self.args = vars(self.parser.parse_args())
		else :
			self.args = vars(self.parser.parse_args(inOpts))


########################################################################
# BED File
########################################################################

class BED12(object):
	'''
	Handles BED format file input and output.

	# BED12 has some built in functions
	# for formally defining elemnts of a bed12
	# and to preform coordinate conversions

	Attribute names are stable, but the value refreshes while iterating through
	the bed file.

	attributes:
	chrom, start, end, strand = reference aligment descriptors
	read = query/feature id
	score = integer 
	c1, c2 = (formally) describe where the open reading frame starts and stop
	exons, size, starts = bed12 alignment blocks descibing where the aligmnents matches in the reference.

	methods:
	getLine - gets line from bed file, and defines values for attributes
	bed12to(Juncs|Exons) - converts bed12 aligmnet blocks to reference coordinate positions

	getLine must be called before bed12to(Juncs|Exons) can be called since it relies on attributes defined in getLine.

	'''
	def __init__(self, fname=None):
		self.fname = fname
		
		if not os.path.isfile(fname):
			print("%s does not exist. Exiting.", file=sys.stderr)
			sys.exit(1)

	def getLine(self):

		with open(self.fname,'r') as entries:
			for entry in entries:
				cols = entry.rstrip().split()
				self.chrom, self.start, self.end, self.read = cols[0], int(cols[1]), int(cols[2]), cols[3]
				self.score, self.strand, self.c1, self.c2 = int(cols[4]), cols[5], int(cols[6]), int(cols[7])
				self.color, self.exons, self.sizes, self.starts = cols[8], int(cols[9]), [int(x) for x in cols[10].split(",")[:-1]], [int(x) for x in cols[11].split(",")[:-1]] 
				yield cols

	def bed12toJuncs(self):
		'''
		Take bed12 entry and convert block/sizes to junction coordinates.
		'''
		junctions = list()
		for num, st in enumerate(self.starts,0):

			if num+1 >= len(self.starts):
				break
			ss1 = self.start + st + self.sizes[num]
			ss2 = self.start + self.starts[num+1]
			junctions.append((ss1,ss2))

		return junctions

	def bed12toExons(self):
		'''
		Take bed12 entry and convert block/sizes to exon coordinates.
		'''
		exons = list()
		for num, st in enumerate(self.starts,0):
			c1 = self.start + st
			c2 = c1 + self.sizes[num]
			exons.append((c1,c2))
		return exons

	 

########################################################################
# Functions
########################################################################

def juncsToBed12(start, end, coords):
	'''
	Take alignment start, end, and junction coords and convert to block/size bed12 format.
	start = integer
	end = integer
	coords = list formatted like so [(j1_left,j1_right),(j2_left,j2_right)]
	'''
	
	sizes, starts = [],[]
	# initial start is 0
	if len(coords) > 0:
		for num,junc in enumerate(coords,0):
			ss1, ss2 = junc
			if num == 0:
				st = 0
				size = abs(start-ss1)
			else:
				st = coords[num-1][1] - start
				size =  ss1 - (st + start)
			starts.append(st)
			sizes.append(size)
		st = coords[-1][1] - start
		size =  end - (st + start)
		starts.append(st)
		sizes.append(size)
		return len(starts), sizes, starts
	else:
		return 1, [end-start], [0] 


def buildOtherDB(ssDB, spliceSites, bedJuncs, wiggle):

	ss = ssDB
	lineNum = 0
	with open(bedJuncs,'r') as bedLines:
		for i in bedLines:
			lineNum += 1


	with open(bedJuncs,'r') as bedLines:
		for line in tqdm(bedLines, total=lineNum, desc="Adding juncs from juncs.bed to interval tree") if verbose else bedLines:
			cols = line.rstrip().split()
			if len(cols)<4:
				print("Other junctions.bed misformatted. Expecting at least 4 cols, got %s" % len(cols), file=sys.stderr)
		
			chrom, c1, c2, strand = cols[0], int(cols[1]), int(cols[2]), cols[-1]
			
			

			if chrom not in ss:
				ss[chrom] = IntervalTree()

			if (chrom,c1) not in spliceSites:
				ss[chrom][c1-wiggle:c1+wiggle] = ('other',c1, strand)
			if (chrom,c2) not in spliceSites:
				ss[chrom][c2-wiggle:c2+wiggle] = ('other',c2, strand)
	return ss

def buildGTFDB(file, wiggle):

	ss = dict()
	exons = dict()
	junctionSet = set()

	if verbose: print("reading gtf %s" % (file), file=sys.stderr) 

	with open(file,'r') as lines:
		for l in lines:
			if l[0] == "#":
				continue

			cols = l.split()

			if "exon" == cols[2]:
				
				chrom, c1, c2, strand =  cols[0], int(cols[3]), int(cols[4]), cols[6]
				
				if c2-c1 < 50:
					continue

				txn = cols[11]
				key = (chrom, txn, strand)
				if key not in exons:
					exons[key] = list()
				exons[key].append((c1,c2))
	
	txnList = list(exons.keys())			
	for exonInfo in tqdm(txnList, total=len(txnList), desc="Building junction interval tree from GTF") if verbose else txnList:
		chrom, txn, strand = exonInfo
		coords = exons[exonInfo]
		if strand == "-":
			coords = coords[::-1]

		if chrom not in ss:
			ss[chrom] = IntervalTree()

		for pos,x in enumerate(coords,0):
			if pos+1<len(coords):
				junctionSet.add((chrom,x[1]))
				
				junctionSet.add((chrom,coords[pos+1][0]-1))
				
				# Need to add +1 to x[1], otherwise i catch the end/begining of the exon.
				ss[chrom][x[1]-wiggle:x[1]+wiggle] = ('gtf',x[1],strand)
				ss[chrom][coords[pos+1][0]-wiggle-1:coords[pos+1][0]+wiggle-1] = ('gtf',coords[pos+1][0]-1,strand)

	return ss, junctionSet

def resolveHits(spliceSite, hits):
	hits = list(hits)


	if len(hits)<1:
		return (np.nan, spliceSite, np.nan, np.nan)

	elif len(hits)<2:
		return(hits[0].data[1],spliceSite,spliceSite-hits[0].data[1], hits[0].data[0])
	else:
		distances = [(abs(spliceSite-hit.data[1]),hit.data[1],hit.data[0]) for hit in hits]
		sortedDist = sorted(distances,key=lambda x: x[0])
		top, second = sortedDist[0], sortedDist[1]

		if top[0] == second[0]:
			best = random.choice([top,second])
			distance, ssCord, refType = best
		else:
			distance, ssCord, refType = top

		return(ssCord,spliceSite,spliceSite-ssCord,refType)	
	

def main():
	'''
	maine
	'''

	# Command Line Stuff...
	myCommandLine = CommandLine()
	bed = myCommandLine.args['bedFile']
	gtf = myCommandLine.args['gtf']
	otherJuncs = myCommandLine.args['junctionsBed']
	wiggle = myCommandLine.args['wiggleWindow']
	threads = myCommandLine.args['threads']
	report = myCommandLine.args['report_junctions']
	

	# There are a few functions that evaluate what verbose is defined as.
	# Instead of passing it around, just global it.
	global verbose
	verbose = myCommandLine.args['quiet']

	# Build interval tree with user splice site input data.
	ssDB, spliceSites = buildGTFDB(gtf, wiggle)
	ssDB = buildOtherDB(ssDB, spliceSites, otherJuncs, wiggle)
	

	# Check how large bed file to be correctd is.
	entries = 0
	with open(bed) as f:
		for l in f:
			entries += 1
			continue


	# BED12 has some built in functions
	# for formally defining elemnts of a bed12
	# and to preform coordinate conversions
	data = BED12(bed)

	statsOut = open("ssCorrectionInfo.tsv",'w') if not os.path.isfile("ssCorrectionInfo.tsv") else open("ssCorrectionInfo%s.tsv" % random.randint(0,9999),'w')

	for line in tqdm(data.getLine(), total=entries, desc="Correcting junctions") if verbose else data.getLine():
		junctionCoords = data.bed12toJuncs()
		if len(junctionCoords)<1:
			continue

		ch, st, end = data.chrom, data.start, data.end
		leftHits = [resolveHits(x[0],ssDB[ch][x[0]]) for x in junctionCoords]
		rightHits = [resolveHits(x[1], ssDB[ch][x[1]]) for x in junctionCoords]

	

		correctedJuncs = [(x[0][0],x[1][0]) for x in zip(leftHits,rightHits)
							if x[0][0] != np.nan and x[1][0] != np.nan]
		if len(correctedJuncs)>0:
			exons, sizes, starts = juncsToBed12(data.start,data.end,correctedJuncs)
			print(data.chrom, data.start, data.end, data.read[:7], 
				data.score, data.strand, data.c1, data.c2, data.color,
				exons, ",".join(map(str,sizes))+",", ",".join(map(str,starts))+",", sep="\t")

		for i in zip(leftHits,rightHits):
			left,right = i
			if data.strand == "+":
				print(data.read[:7], ch, "\t".join(map(str,left)), "5'", data.strand, sep="\t", file=statsOut)
				print(data.read[:7], ch, "\t".join(map(str,right)), "3'", data.strand, sep="\t", file=statsOut)
			else:
				print(data.read[:7], ch, "\t".join(map(str,left)), "3'",data.strand,  sep="\t", file=statsOut)
				print(data.read[:7], ch, "\t".join(map(str,right)), "5'",data.strand, sep="\t", file=statsOut)

	statsOut.close()

if __name__ == "__main__":
	main()
