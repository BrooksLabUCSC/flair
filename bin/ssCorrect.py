from __future__ import print_function

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
import multiprocessing
from multiprocessing import Pool
from contextlib import closing
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
                                             usage = '%(prog)s -i reads.bed -g annotations.gtf -j other_junctions.bed -o out_file.bed')
        # Add args
        self.parser.add_argument('-i', "--input_bed", action = 'store', required=True, help='Input reads in bed12 format.')
        self.parser.add_argument('-g', "--gtf", action = 'store', required=False, help='Gencode annotation file.')
        self.parser.add_argument('-j', "--junctionsBed", default=None, action = 'store', required=False, help='Short-read supported junctions in bed6 format (Optiona) [BED entries must be UNIQUE and have strand information].')
        self.parser.add_argument('-w', '--wiggleWindow', action = 'store', type=int, required=False, default = 15, help='Splice site correction window flank size.')
        self.parser.add_argument('-o', "--output_fname", action = 'store', required=True, help='Output file name.')
        
        # (Under Development!!!)
        #self.parser.add_argument('--report_junctions', action = 'store', required=False, default= "norm", choices=['strict', 'norm', 'all'],
        #                                                                                help = """Choose which types of nanopore splice sites to report:
        #                                                                                        "strict" - SS must be in gtf or junctionsBed, and have a single closest hit.
        #                                                                                        "norm" - SS must be in gtf or junctionsBed (tie for closest hit is random).
        #                                                                                        "all" - SS can be in gtf, junctionsBed, or unique to nanopore data.""")
        self.parser.add_argument('-p', "--threads", action = 'store', required=False, type=int, default = 2, help='Number of threads.')
        self.parser.add_argument('--keepZero', action = 'store_true', required=False, default = False, help='Keep alignments with no spliced junctions (single exon txns).')
        self.parser.add_argument("--quiet", action = 'store_false', required=False, default = True, help='Do not display progress')
        self.parser.add_argument("--cleanup", action = 'store_false', required=False, default = True, help='Remove teomprary files with correction info.')
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
                self.chrom, self.start, self.end, self.name = cols[0], int(cols[1]), int(cols[2]), cols[3]
                self.score, self.strand, self.c1, self.c2 = int(cols[4]), cols[5], int(cols[6]), int(cols[7])
                self.color, self.exons = cols[8], int(cols[9]),  
                self.sizes =  [int(x) for x in cols[10].split(",")[:-1]] if cols[10][-1] == "," else [int(x) for x in (cols[10]+",").split(",")[:-1]]
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

    def determinType(self):
        pass


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
    novelFlag = False
    sizes, starts = [],[]
    # initial start is 0
    if len(coords) > 0:
        for num,junc in enumerate(coords,0):
            ss1, ss2, novel1, novel2 = junc
            if novel1 or novel2: novelFlag = True
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
        return len(starts), sizes, starts, novelFlag
    else:
        return 1, [end-start], [0], novelFlag 


def buildOtherDB(acceptors, donors, spliceSites, bedJuncs, wiggle, proc):

    
    lineNum = 0

    with open(bedJuncs,'r') as bedLines:
        lineNum +=1
        cols = next(bedLines).split()

        for i in bedLines:
            lineNum += 1

    # guess what kind of bedFile
    if cols[-1] == "+" or cols[-1] == "-":
        # normal bed
        reverseSS = "-"
        strandCol = -1
        starOffset = 0

    elif len(cols) == 12:
        # bed12
        bedType   = "bed12"
        print("ERROR: Bed12 not currently supported for other_juncs.bed. Please convert to bed6. Exiting.", file=sys.stderr)

    elif cols[3] == "0" or cols[3] == "1" or cols[3] == "2":
        # star junc.tab
        reverseSS = "2"
        strandCol = 3
        starOffset = 1

    else:
        print("ERROR: Cannot find strand info for %s. Is this bed6 or STAR_juncs.tab file? Exit." % bedJuncs, file=sys.stderr)

    proc += 1 
    with open(bedJuncs,'r') as bedLines:
        for line in tqdm(bedLines, total=lineNum, desc="Adding splice sites from %s  to interval tree" % bedJuncs, dynamic_ncols=True, position = proc) if verbose else bedLines:
            cols = line.rstrip().split()
            chrom, c1, c2, strand = cols[0], int(cols[1])-starOffset, int(cols[2]), cols[strandCol]
            if strand == reverseSS: c1,c2 = c2,c1

            if chrom not in donors:
                donors[chrom] = IntervalTree()
                acceptors[chrom] = IntervalTree()

            if (chrom,c1) not in spliceSites:
                donors[chrom][c1-wiggle:c1+wiggle] = ('other',c1, strand)
            if (chrom,c2) not in spliceSites:
                acceptors[chrom][c2-wiggle:c2+wiggle] = ('other',c2, strand)
    return donors, acceptors

def buildGTFDB(file, wiggle, proc):

    exons       = dict()
    junctionSet = set()

    if verbose: print("reading gtf %s ..." % (file), file=sys.stderr) 

    with open(file,'r') as lines:
        for l in lines:
            if l[0] == "#": # skip header lines
                continue

            cols = l.split("\t")

            if "exon" == cols[2]:
                
                # -1 for 1 to 0 based conversion
                chrom, c1, c2, strand =  cols[0], int(cols[3])-1, int(cols[4]), cols[6]

                #txn info is in the SECOND position of the shoutout column
                txn = cols[-1].split(";")[1].split()[-1].replace('"','')
                key = (chrom, txn, strand)
                if key not in exons:
                    exons[key] = list()
                exons[key].append(c1)
                exons[key].append(c2)

    txnList = list(exons.keys())
    proc += 1            
    for exonInfo in tqdm(txnList, total=len(txnList), desc="Building junction interval tree from GTF", dynamic_ncols=True, position=proc) if verbose else txnList:
        chrom, txn, strand = exonInfo

        coords = exons[exonInfo]
        
        # assume lowest and highest as TSS and TES, and remove them
        coords.sort()
        coords = coords[1:-1]

        if len(coords)<2: continue

        if chrom not in donors: #make the interval tree
            donors[chrom]    = IntervalTree()
            acceptors[chrom] = IntervalTree()
        
        for pos in range(0,len(coords)-1,2):
            c1 = coords[pos]
            c2 = coords[pos+1]
            
            if strand == "-": c1,c2 = c2,c1 
            
            if (chrom,c1) not in junctionSet:
                donors[chrom][c1-wiggle:c1+wiggle] = ('gtf',c1,strand)
                junctionSet.add((chrom,c1))
            if (chrom,c2) not in junctionSet:
                acceptors[chrom][c2-wiggle:c2+wiggle] = ('gtf',c2,strand)
                junctionSet.add((chrom,c2))
   
    return donors, acceptors, junctionSet

def resolveHits(spliceSite, hits, read):
    hits = list(hits)


    if len(hits)<1:
        # novel
        return (spliceSite, spliceSite, np.nan, np.nan, True)

    elif len(hits) == 1:
        # single closest hit
        return (hits[0].data[1],spliceSite,spliceSite-hits[0].data[1], hits[0].data[0], False)
    else:
        # multi hit
        distances = [(abs(spliceSite-hit.data[1]),hit.data[1],hit.data[0]) for hit in hits]
        sortedDist = sorted(distances,key=lambda x: x[0])
        top, second = sortedDist[0], sortedDist[1]

        if top[0] == second[0]:
            best = random.choice([top,second])
            distance, ssCord, refType = best
        else:
            distance, ssCord, refType = top

        return (ssCord,spliceSite,spliceSite-ssCord,refType, False)
    


def ssCorrect(chrom, bedFile, fileSize, procNum):
    
    data         = BED12(bedFile)
    statsOut     = open("%s_ssCorrectionInfo.tsv" % chrom,'w')
    tempOut      = open("%s_corrected.bed" % chrom ,'w')
    tempNovelOut = open("%s_uncorrected.bed" % chrom ,'w')

    for line in tqdm(data.getLine(), total=fileSize, desc="Working on %s" % chrom, position=procNum,dynamic_ncols=True) if verbose else data.getLine():
        junctionCoords = data.bed12toJuncs()

        ch, st, end, blocks, strand = data.chrom, data.start, data.end, data.exons, data.strand
        readID = data.name[:18]


        if ch in donors and ch in acceptors:
            if strand == "+":
                hits = [ (resolveHits(x[0],donors[ch][x[0]],readID), resolveHits(x[1], acceptors[ch][x[1]],readID))
                                 for x in junctionCoords]
            else:
                hits = [ (resolveHits(x[0],acceptors[ch][x[0]],readID), resolveHits(x[1], donors[ch][x[1]],readID))
                                 for x in junctionCoords]

            correctedJuncs = [(x[0][0],x[1][0],x[0][-1],x[1][-1]) for x in hits]
        else:
            # un annotated reference chrom/contig
            if len(junctionCoords) < 1:
                if keepZero:
                     print(data.chrom, data.start, data.end, readID, 
                        data.score, data.strand, data.c1, data.c2, data.color,
                        data.exons, "%s," % data.sizes[0], "%s," % data.starts[0], sep="\t",file=tempOut)
                else:
                     print(data.chrom, data.start, data.end, readID, 
                        data.score, data.strand, data.c1, data.c2, data.color,
                        data.exons, "%s," % data.sizes[0], "%s," % data.starts[0], sep="\t",file=tempNovelOut)
                    
            else:
                print(data.chrom, data.start, data.end, readID, 
                    data.score, data.strand, data.c1, data.c2, data.color,
                    data.exons, "%s," % data.sizes[0], "%s," % data.starts[0], sep="\t",file=tempNovelOut)
            continue

        if len(junctionCoords) < 1:
            if keepZero:
                print(data.chrom, data.start, data.end, readID,
                data.score, data.strand, data.c1, data.c2, data.color,
                data.exons, "%s," % data.sizes[0], "%s," % data.starts[0], sep="\t",file=tempOut)
            continue

        else:

            exons, sizes, starts, novelJuncs = juncsToBed12(data.start,data.end,correctedJuncs)
            if novelJuncs:
                print(data.chrom, data.start, data.end, readID, 
                    data.score, data.strand, data.c1, data.c2, data.color,
                    exons, ",".join(map(str,sizes))+",", ",".join(map(str,starts))+",", sep="\t",file=tempNovelOut)

            else:
                print(data.chrom, data.start, data.end, readID, 
                    data.score, data.strand, data.c1, data.c2, data.color,
                    exons, ",".join(map(str,sizes))+",", ",".join(map(str,starts))+",", sep="\t",file=tempOut)

            for i in hits:
                left,right = i

                if data.strand == "+":
                    print(readID, ch, "\t".join(map(str,left)), "5'", data.strand, sep="\t", file=statsOut)
                    print(readID, ch, "\t".join(map(str,right)), "3'", data.strand,  sep="\t", file=statsOut)
                else:
                    print(readID, ch, "\t".join(map(str,left)), "3'",data.strand,  sep="\t", file=statsOut)
                    print(readID, ch, "\t".join(map(str,right)), "5'",data.strand, sep="\t", file=statsOut)

    os.remove("%s_reads.temp.bed" % chrom)
    statsOut.close()
    tempOut.close()
    tempNovelOut.close()
    return ("%s_ssCorrectionInfo.tsv" % chrom, "%s_corrected.bed" % chrom, "%s_uncorrected.bed" % chrom)


def runCMD(x):
    chrom, bedFile, fileSize, procNum = x
    
    return ssCorrect(chrom, bedFile, fileSize,procNum)

def main():
    '''
    maine
    '''

    # Command Line Stuff...
    myCommandLine = CommandLine()
    bed           = myCommandLine.args['input_bed']
    gtf           = myCommandLine.args['gtf']
    otherJuncs    = myCommandLine.args['junctionsBed']
    wiggle        = myCommandLine.args['wiggleWindow']
    threads       = myCommandLine.args['threads']
    out           = myCommandLine.args['output_fname']
    cleanup       = myCommandLine.args['cleanup']

    # keep track of process number
    proc = -1

    global keepZero
    keepZero      = myCommandLine.args['keepZero']
    
    #report        = myCommandLine.args['report_junctions']

    # There are a few functions that evaluate what verbose is defined as.
    # Instead of passing it around, just global it.
    global verbose
    verbose = myCommandLine.args['quiet']

    # Build interval tree with user splice site input data.
    global donors, acceptors
    donors, acceptors = dict(), dict()
    donors, acceptors, spliceSites = buildGTFDB(gtf, wiggle, proc)

    if otherJuncs != None: donors, acceptors = buildOtherDB(donors, acceptors, spliceSites, otherJuncs, wiggle, proc)



    # Split bed by CHROMOSOME!!!!!
    # and Check how large bed file to be correctd is.
    fileDict = dict()
    fileSizeDict = dict()
    proc += 1
    with open(bed) as f:
        for l in tqdm(f, desc="Splitting isoforms by chrom for multiprocessing", position=proc,dynamic_ncols=True) if verbose else f:
            chrom = l.split()[0]
            if chrom not in fileDict:
                fileDict[chrom] = open("%s_reads.temp.bed" % chrom,'w')
                fileSizeDict[chrom] = 0
            print(l.rstrip(), file=fileDict[chrom])
            fileSizeDict[chrom] += 1
            

    cmdList = list()
    files = sorted(list(fileDict.keys()))
    for key in files:
        proc += 1
        cmdList.append((key, "%s_reads.temp.bed" % key, fileSizeDict[key], proc))
        fileDict[key].close()   

    
    
    # BED12 has some built in functions
    # for formally defining elemnts of a bed12
    # and to preform coordinate conversions

    finalFiles = list()
    proc += 1
    with closing(Pool(threads)) as p:
        for i in tqdm(p.imap(runCMD, cmdList), total=len(cmdList), desc="Correcting junctions per chromosome", position=proc,dynamic_ncols=True) if verbose else p.imap(runCMD, cmdList):
            finalFiles.append(i)
    proc += 1

    allUncorrected = open("uncorrected_%s" % out, 'w')
    with open(out,'w') as outF:
        for i in tqdm(finalFiles, total=len(finalFiles), desc="Writing corrected juncs to %s" % out, dynamic_ncols=True, position = proc) if verbose else finalFiles:
            stats,corrected, uncorrected = i
            with open(corrected) as lines:
                [print(x.rstrip(), file=outF) for x in lines]
            with open(uncorrected) as lines:
                [print(x.rstrip(), file=allUncorrected) for x in lines]
            os.remove(uncorrected)
            os.remove(corrected)
            if cleanup:
                os.remove(stats)
            


if __name__ == "__main__":
    main()
