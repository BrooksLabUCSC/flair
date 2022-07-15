#!/usr/bin/env python3

from __future__ import print_function


########################################################################
# File: bam2Bed12.py
#  executable: bam2Bed12.py
# Purpose: Conver minimap2 aligned bam file to bed12
#
#          
# Author: Cameron M. Soulette
# History:      cms 03/22/2018 Created
#
########################################################################


########################################################################
# Hot Imports & Global Variable
########################################################################


import os, sys
from samJuncs import SAM
import pysam

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
        self.parser = argparse.ArgumentParser(description = 'A tool to convert minimap2 BAM to Bed12.',
                                            add_help = True, #default is True 
                                            prefix_chars = '-', 
                                            usage = '%(prog)s -i sorted.aligned.bam ')
        # Add args
        self.parser.add_argument('-i', "--input_bam", action = 'store', required=True, help='Input bam file.')
        self.parser.add_argument('--keep_supplementary', action = 'store_true', required=False, default=False,  help='Keep supplementary alignments')

        if inOpts is None :
            self.args = vars(self.parser.parse_args())
        else :
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
                size =  ss1 - (st + start)
            starts.append(st)
            sizes.append(size)

        # Here is the computation for the BED end coordinate
        st = coords[-1][1] - start
        size =  end - (st + start)
        starts.append(st)
        sizes.append(size)


        return len(starts), sizes, starts
    else:
        return 1, [end-start], [0] 


########################################################################
# Main
# Here is the main program
# 
########################################################################
def main():
    '''
    stuff...
    '''

    myCommandLine = CommandLine()
    
    alignmentFile = myCommandLine.args['input_bam']

    #Color codes for positive and negative stranded read transcripts
    positiveTxn = "27,158,119"
    negativeTxn = "217,95,2"
    unknownTxn = "99,99,99"

    # SAM Object allows for execution of many SAM-related functions.
    sObj = SAM(alignmentFile, keep_supplementary = myCommandLine.args['keep_supplementary'])


    for num, readData in enumerate(sObj.readJuncs(),0):
        read, chrom, startPos, junctions, endPos, flags, tags, score = readData
        blocks, sizes, starts = juncsToBed12(startPos, endPos, junctions)
        flags = str(flags)

        if tags == "+":

            print(chrom, startPos, endPos, read + ";" + flags , score, tags, startPos, endPos, positiveTxn, blocks, 
                ",".join(str(x) for x in sizes) + ",", ",".join(str(x) for x in starts) + ",", sep="\t")
        elif tags == "-":
            print(chrom, startPos, endPos, read + ";" + flags , score, tags, startPos, endPos, negativeTxn, blocks, 
                ",".join(str(x) for x in sizes) + ",", ",".join(str(x) for x in starts) + ",", sep="\t")                    

        else:
            tags = "+" if flags == "0" else "-"
            print(chrom, startPos, endPos, read + ";" + flags , score, tags, startPos, endPos, unknownTxn, blocks, 
            ",".join(str(x) for x in sizes) + ",", ",".join(str(x) for x in starts) + ",", sep="\t")                    


if __name__ == "__main__":
    main();        
