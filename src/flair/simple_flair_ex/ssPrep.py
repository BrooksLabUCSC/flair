#!/usr/bin/env python3

import os
import sys
from ncls import NCLS

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

    def __init__(self):
        '''
        CommandLine constructor.
        Implements a parser to interpret the command line argv string using argparse.
        '''
        import argparse
        self.parser = argparse.ArgumentParser(description=' ssPrep.py - a tool to leverage annotation and short read data to correct misaligned splice junctions in short read data.',
                                             add_help=True, #default is True
                                             prefix_chars='-',
                                             usage='%(prog)s -i reads.bed -j known_junctions.bed -o out_file.bed --working_dir dir')
        # Add args
        self.parser.add_argument('-i', "--input_bed", required=True, help='Input reads in bed12 format.')
        self.parser.add_argument('-j', "--juncs", required=True, help='KnownJunction.bed.')
        self.parser.add_argument('-w', '--wiggleWindow', type=int, required=False, default=15, help='Splice site correction window flank size (default 15).')
        self.parser.add_argument('-o', "--output_fname", required=True, help='Output file name.')
        self.parser.add_argument('-f', "--genome_fasta", required=True, help='Genome Fasta.')
        self.parser.add_argument("--workingDir", required=True, help='Working directory.')
        self.parser.add_argument('--correctStrand', action='store_true', required=False, default=False, help='Try to resolve read strand by using annotated splice site strand.')
        self.parser.add_argument('--check_file', required=False, default=False, help='Write file for print_check')

        #self.parser.add_argument('--keepZero', action='store_true', required=False, default=False, help='Keep alignments with no spliced junctions (single exon txns).')
# TODO: is juncs always a real bed file?

        if len(sys.argv) == 1:
                self.parser.print_help()
                sys.exit(1)
        self.args = vars(self.parser.parse_args())


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
                self.chrom, self.start, self.end, self.name = cols[0], int(cols[1]), int(cols[2]), cols[3]
                self.score, self.strand, self.c1, self.c2 = int(cols[4]), cols[5], int(cols[6]), int(cols[7])
                self.color, self.exons = cols[8], int(cols[9])
                self.sizes = [int(x) for x in cols[10].split(",")[:-1]] if cols[10][-1] == "," else [int(x) for x in (cols[10]+",").split(",")[:-1]]
                self.starts = [int(x) for x in cols[11].split(",")[:-1]] if cols[11][-1] == "," else [int(x) for x in (cols[11]+",").split(",")[:-1]]
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
# Junc
########################################################################


class junctObj(object):
    '''
    Handles Junc data.
    '''
    def __init__(self, coord=None, strand=None, ssType=None):
        self.coord   = coord
        self.strand  = strand
        self.ssType  = ssType # donor or acceptor

        # Descriptive attributes.
        self.support = set()
        self.ssCorr  = None


########################################################################
# Functions
########################################################################


def juncsToBed12(start, end, coords):
    '''
    Take alignment start, end, and junction coords and convert to block/size bed12 format.
    start = integer
    end = integer
    coords = list formatted like so [(j1_left,j1_right),(j2_left,j2_right)]
    returns num_exons, sizes, starts
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
                size = ss1 - (st + start)
            starts.append(st)
            sizes.append(size)
        st = coords[-1][1] - start
        size = end - (st + start)
        starts.append(st)
        sizes.append(size)
        return len(starts), sizes, starts
    else:
        return 1, [end-start], [0]


def ssCorrect(c,strand,ssType,intTree,junctionBoundaryDict, errFile):
    '''
    correct un-annotated splice sites.
    '''

    if intTree == False:
        return junctionBoundaryDict

    hits = [h for h in intTree.find_overlap(c,c)]

    if len(hits) < 1:
        ss = junctObj(c,strand,None)
        junctionBoundaryDict[c] = ss
        ss.ssCorr = ss
        return junctionBoundaryDict
    else:

        distances = [abs(c-x[-1]) for x in hits]
        minVal    = min(distances)
        count     = distances.count(minVal)

        if count > 1: ###multiple splice sites with equal distance to read ss, how to pick?

            sortedvals = []
            for x in range(len(hits)):
                if distances[x] == minVal:
                    annot = junctionBoundaryDict[hits[x][-1]]
                    sortedvals.append(('both' in annot.support or 'gtf' in annot.support, annot.strand == strand, hits[x][-1]))
            sortedvals.sort()

            # with open(errFile, 'a+') as fo:
            #     print(c, hits, sortedvals, file=fo)
            cCorr = sortedvals[-1][-1]
        else:
            cCorr = hits[distances.index(minVal)][-1]
        ss = junctObj(c,strand,ssType)
        ss.ssCorr = junctionBoundaryDict[cCorr]
        junctionBoundaryDict[c] = ss
        return junctionBoundaryDict


        #correctReads(readsBed, intervalTree, junctionBoundaryDict, chrom, resolveStrand, workingDir, chrom, errFile)
def correctReads(readsBed, intTree, junctionBoundaryDict, filePrefix, correctStrand, wDir, currentChr, errFile):
    ''' Builds read and splice site objects '''

    if errFile:
        with open(errFile,'a+') as fo:
            print("** Creating temporary correction files for chromosome %s: %s & %s" % (currentChr, os.path.join(wDir, "%s_inconsistent.bed" % filePrefix), os.path.join(wDir,"%s_corrected.bed" % filePrefix)), file=fo)

    inconsistent = open(os.path.join(wDir, "%s_inconsistent.bed" % filePrefix),'w')
    corrected = open(os.path.join(wDir,"%s_corrected.bed" % filePrefix),'w')

    bedObj = BED12(readsBed)
    for line in bedObj.getLine():
        juncs  = bedObj.bed12toJuncs()
        strand = bedObj.strand
        c1Type,c2Type = ("donor","acceptor") if strand == "+" else ("acceptor","donor")

        newJuncs  = list()
        ssTypes   = list()
        ssStrands = set()
        novelSS   = False

        for x in juncs:
            c1, c2 = x[0], x[1]
            if c1 not in junctionBoundaryDict:
                junctionBoundaryDict = ssCorrect(c1,strand,c1Type,intTree,junctionBoundaryDict, errFile)
            if c2 not in junctionBoundaryDict:
                junctionBoundaryDict = ssCorrect(c2,strand,c2Type,intTree,junctionBoundaryDict, errFile)

            #c1Obj, c2Obj = junctionBoundaryDict[c1], junctionBoundaryDict[c2] unused

            c1Corr = junctionBoundaryDict[c1].ssCorr.coord
            c2Corr = junctionBoundaryDict[c2].ssCorr.coord

            ssTypes = [junctionBoundaryDict[c1].ssCorr.ssType, junctionBoundaryDict[c2].ssCorr.ssType]

            # with open(errFile, 'a+') as fo:
            #     print(c1, c1Corr, ssTypes[0], c2, c2Corr, ssTypes[1], None in ssTypes, file=fo)

            ssStrands.add(junctionBoundaryDict[c1].ssCorr.strand)
            ssStrands.add(junctionBoundaryDict[c2].ssCorr.strand)

            if None in ssTypes: #or ssTypes[0] == ssTypes[1]:
                # Either two donors or two acceptors or both none.
                novelSS = True

            newJuncs.append((c1Corr,c2Corr))

        blocks, sizes, starts = juncsToBed12(bedObj.start,bedObj.end,newJuncs)

        if correctStrand:
            # if len(ssStrands) > 1:
            #     novelSS = True
            if len(ssStrands) == 1:
                strand = list(ssStrands)[0]
            elif len(ssStrands) == 0:
                strand = strand

        # 0 length exons, remove them.
        minSize = min(sizes)
        if minSize == 0: novelSS = True

        # with open(errFile, 'a+') as fo:
        #     # print(correctStrand, len(ssStrands) > 1, minSize == 0, file=fo)
        #     print(minSize == 0, file=fo)

        if novelSS:
            print(bedObj.chrom, bedObj.start, bedObj.end, bedObj.name,
                    bedObj.score, strand, bedObj.c1, bedObj.c2, bedObj.color,
                    blocks, ",".join(map(str,sizes))+",", ",".join(map(str,starts))+",", sep="\t", file=inconsistent)
        else:

            print(bedObj.chrom, bedObj.start, bedObj.end, bedObj.name,
                    bedObj.score, strand, bedObj.c1, bedObj.c2, bedObj.color,
                    blocks, ",".join(map(str,sizes))+",", ",".join(map(str,starts))+",", sep="\t", file=corrected)
    corrected.close()
    inconsistent.close()

    inc = os.path.isfile(os.path.join(wDir, "%s_inconsistent.bed" % filePrefix))
    cor = os.path.isfile(os.path.join(wDir, "%s_corrected.bed" % filePrefix))

    if errFile:
        with open(errFile,'a+') as fo:
            print("** Checking inc/corr files for chromsome %s: %s %s" % (currentChr,inc,cor), file=fo)


def buildIntervalTree(knownJuncsFile, wiggle, currentChr, errFile):
    ''' Builds read and splice site objects '''

    if errFile:
        with open(errFile,'a+') as fo:
            print("** Initializing int tree for chromosome %s" % (currentChr), file=fo)

    intronranges = [] # holds list of windowStart, windowEnd, junction
    junctionBoundaryDict = dict() # holds junction objects

    # the below adds support for splice boundaries, found in the 4th column of the known junctions file
    # extract junctions from the file
    with open(knownJuncsFile, 'r') as lines:
        for line in lines:
            cols     = line.rstrip().split()
            intronBoundary1, intronBoundary2   = int(cols[1]), int(cols[2])
            strand   = cols[-1]

            annoType = cols[3] # 'sr' in example TODO what else?
            intronBoundary1Type,intronBoundary2Type = ("donor","acceptor") if strand == "+" else ("acceptor","donor")

            # add intronBoundary1 first
            if intronBoundary1 not in junctionBoundaryDict:
                ss = junctObj(intronBoundary1,strand,intronBoundary1Type)
                ss.support.add(annoType)
		# TODO: turtles all the way down?! Why.
                ss.ssCorr = ss

                # SS window
                junctionWindowStart, junctionWindowEnd = max(intronBoundary1-wiggle,1), intronBoundary1+wiggle

                # Add to tree and object to junctionBoundaryDict
                junctionBoundaryDict[intronBoundary1] = ss
                intronranges.append([junctionWindowStart,junctionWindowEnd,intronBoundary1])

            else:
                junctionBoundaryDict[intronBoundary1].support.add(annoType)

            # now add intronBoundary2
            if intronBoundary2 not in junctionBoundaryDict:
                ss = junctObj(intronBoundary2,strand,intronBoundary2Type)
                ss.support.add(annoType)
                ss.ssCorr = ss

                # SS window
                junctionWindowStart, junctionWindowEnd = max(intronBoundary2-wiggle,1), intronBoundary2+wiggle

                # Add to tree and object to junctionBoundaryDict
                junctionBoundaryDict[intronBoundary2] = ss
                intronranges.append([junctionWindowStart,junctionWindowEnd,intronBoundary2])
            else:
                junctionBoundaryDict[intronBoundary2].support.add(annoType)
    # ncls does not allow empty creation, so make one nonsense entry
    intervalTree = NCLS([-1], [-1], [-1])
    if len(intronranges) > 0:
	# x is a new reference to the same dict
        x = intronranges
        intervalTree = NCLS([val[0]-1 for val in x], [val[1]+1 for val in x], [val[2] for val in x])
    if errFile:
        with open(errFile,'a+') as fo:
            print("** Tree Initialized. %s data points added for chromosome %s." % (len(list(junctionBoundaryDict.keys())),currentChr), file=fo)

    return intervalTree, junctionBoundaryDict


def main():
    '''
    main
    '''

    # Command Line Stuff...
    myCommandLine = CommandLine()
    readsBed        = myCommandLine.args['input_bed']
    knownJuncsFile  = myCommandLine.args['juncs']
    genomeFa            = myCommandLine.args['genome_fasta']
    wiggle        = myCommandLine.args['wiggleWindow']
    out           = myCommandLine.args['output_fname']
    resolveStrand = myCommandLine.args['correctStrand']
    workingDir    = myCommandLine.args['workingDir']
    errFile    = myCommandLine.args['check_file']

    # WARNING: out is a chromosome name and interpreted as such.
    ssPrep([readsBed, knownJuncsFile, genomeFa, wiggle, out, resolveStrand, workingDir, errFile])


def ssPrep(x):
    '''one argument so we can run p.map from the main program'''
    readsBed, knownJuncsFile, genomeFa, wiggle, chrom, resolveStrand, workingDir, errFile = x

    if errFile:
        with open(errFile,'a+') as fo:
            print("** Correcting %s with a wiggle of %s against %s. Checking splice sites with genome %s." % (readsBed, wiggle, knownJuncsFile, genomeFa), file=fo)

    # Build interval tree of known juncs
    intervalTree, junctionBoundaryDict = buildIntervalTree(knownJuncsFile, wiggle, chrom, errFile)

    if errFile:
        with open(errFile,'a+') as fo:
            print("** SS Correction DB for  %s against %s Built. Moving to correction. Writing files to " % (knownJuncsFile, readsBed), file=fo)
    # Build read objects.
    try:
        correctReads(readsBed, intervalTree, junctionBoundaryDict, chrom, resolveStrand, workingDir, chrom, errFile)
    except Exception as ex:
        raise Exception("** correctReads FAILED for %s" % (readsBed)) from ex



if __name__ == "__main__":
    main()
