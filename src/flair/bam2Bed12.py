#!/usr/bin/env python3

import sys
from flair.samJuncs import SAM

########################################################################
# CommandLine
########################################################################


class CommandLine(object):
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

	def __init__(self, inOpts=None):
		'''
		CommandLine constructor.
		Implements a parser to interpret the command line argv string using argparse.
		'''
		import argparse
		self.parser = argparse.ArgumentParser(description='A tool to convert minimap2 BAM to Bed12.',
											add_help=True, #default is True
											prefix_chars='-',
											usage='%(prog)s -i sorted.aligned.bam ')
		# Add args
		self.parser.add_argument('-i', "--input_bam", action='store', required=True, help='Input bam file.')
		self.parser.add_argument('--keep_supplementary', action='store_true', required=False, default=False, help='Keep supplementary alignments')

		if inOpts is None:
			self.args = vars(self.parser.parse_args())
		else:
			self.args = vars(self.parser.parse_args(inOpts))

########################################################################
# Functions
########################################################################


def juncsToBed12(start, end, coords):
	'''
	junctToBed12 takes in alignment start position, end position, and genomic junction coordinates
	and converts them to start, end, and length blocks for bed12.
	'''

	sizes, starts = [],[]

	#coords with 0 length are reads without introns
	if len(coords) > 0:
		for num,junc in enumerate(coords,0):
			# a junc is 2 Splice Sites
			ss1, ss2 = junc

			# initial start is 0
			if num == 0:
				st = 0
				size = abs(start-ss1)
			else:
				st = coords[num-1][1] - start
				size = ss1 - (st + start)
			starts.append(st)
			sizes.append(size)

		# Here is the computation for the BED end coordinate
		st = coords[-1][1] - start
		size = end - (st + start)
		starts.append(st)
		sizes.append(size)

		return len(starts), sizes, starts
	else:
		return 1, [end-start], [0]


def bam2Bed12(alignmentFile, bedoutput='/dev/stdout', keep_supplementary=False):

	#Color codes for positive and negative stranded read transcripts
        #To show in UCSC browser add itemRgb=On to bed header
	positiveTxn = "27,158,119" # green
	negativeTxn = "217,95,2"   # orange
	unknownTxn = "99,99,99"

	# SAM Object allows for execution of many SAM-related functions.
	sObj = SAM(alignmentFile, keep_supplementary=keep_supplementary)

	with open(bedoutput, 'w') as outf:
		for num, readData in enumerate(sObj.readJuncs(),0):
			read, chrom, startPos, junctions, endPos, flags, juncDirection, score = readData
			blocks, sizes, starts = juncsToBed12(startPos, endPos, junctions)
			flags = str(flags)

			# set correct color
                        # Bed format needs a strand. The readJuncs function infers strand based on
                        # alignment orientation and splice direction (ts tag in SAM) if possible, or
                        # is set to None if not
			rgbcolor = unknownTxn
			if juncDirection == "+":
				rgbcolor = positiveTxn
			elif juncDirection == "-":
				rgbcolor = negativeTxn
			else:
                                # if the juncDirection is None, use the SAM flag but don't change the color
                                # the SAM object only retains reads with flag 0 or 16
				juncDirection = "+" if flags == "0" else "-"

			sizes = ",".join(str(x) for x in sizes) + ','
			starts = ",".join(str(x) for x in starts) + ','

			outitems = [chrom, startPos, endPos, read, score, juncDirection, startPos,
			    	endPos, rgbcolor, blocks, sizes, starts]
			outstring = ('\t').join([str(x) for x in outitems])
			outf.write(f'{outstring}\n')
	outf.close()

if __name__ == "__main__":
	'''When bam2Bed12 is called as a standalone program it prints to STDOUT'''
	myCommandLine = CommandLine()
	bam2Bed12(alignmentFile = myCommandLine.args['input_bam'],
	   keep_supplementary = myCommandLine.args['keep_supplementary'])
