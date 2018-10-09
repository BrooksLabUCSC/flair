#!/usr/bin/env python3


########################################################################
# File: annotatCollapsed.py
#  executable: annotateCollapsed
# Purpose: 
#
#          
# Author: Cameron M. Soulette
# History:      cms 10/04/2018 Created
#
########################################################################


########################################################################
# Hot Imports & Global Variable
########################################################################


import os, sys
from tqdm import *
from ssCorrect import BED12
import random

random.seed(15)

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
                                             usage = '%(prog)s [options] -i isoforms.bed -g annotations.gtf')
        # Add args
        self.parser.add_argument('-i', "--input_isoforms", action = 'store', required=True, help='Input collapsed isoforms in bed12 format.')
        self.parser.add_argument('-g', "--gtf", action = 'store', required=False, help='Gencode annotation file.')
        self.parser.add_argument("--quiet", action = 'store_false', required=False, default = True, help='Do not display progress')
        if inOpts is None :
            self.args = vars(self.parser.parse_args())
        else :
            self.args = vars(self.parser.parse_args(inOpts))

# ########################################################################
# # Gene
# ########################################################################

class Gene(object):
    def __init__(self, name=None, strand=None, chrom=None):
        self.gid = name
        self.chrom = chrom
        self.strand = strand
        self.transcripts = dict()
        self.exons = dict()
        

########################################################################
# Transcript
########################################################################

class Transcript(object):
    def __init__(self, parent=None, name=None, strand=None):
        self.tid = name
        self.exons = list()
        self.parent=parent

########################################################################
# Exon
########################################################################

class Exon(object):
    def __init__(self, parent=None, start=None, end=None):
        self.c1 = start
        self.c2 = end
        self.parent = parent


########################################################################
# Funktions
########################################################################
def buildGTFDB(gtf):

    genes = dict()
    if not os.path.isfile(gtf):
        print("%s GTF Not found. exiting..." % gtf, file=sys.stderr)

    # Verbosity stuff
    print("Reading annotation file - %s" % gtf, file=sys.stderr) if verbose else True
        
    # build gene objex
    with open(gtf) as lines:
        for line in lines:

            cols = line.rstrip().split("\t")
            if len(cols)<3 or cols[2] != "exon" :
                continue

            chrom, c1, c2, strand = cols[0], int(cols[3]), int(cols[4]), cols[6]

            data = cols[-1].split(";")
            gene = data[4].split()[-1].replace('"','')
            #txn = data[7].split()[-1].replace('"','') #hugo tid
            txn = data[1].split()[-1].replace('"','') #ensembl tid
            exon = (chrom, c1, c2)
            
            if gene not in genes:
                genes[gene] = Gene(gene, strand, chrom)

            geneObj = genes[gene]

            if txn not in geneObj.transcripts:
                geneObj.transcripts[txn] = Transcript(geneObj, txn, geneObj)

            if exon not in geneObj.exons:
                geneObj.exons[exon] = Exon(geneObj, c1, c2)

            exonObj = geneObj.exons[exon]
            txnObj = geneObj.transcripts[txn]

            txnObj.exons.append(exonObj)

    juncRef = {gObj.chrom:dict() for g,gObj in genes.items()}
    geneRef = {gObj.chrom:dict() for g,gObj in genes.items()}
    for gene in tqdm(genes, total=len(list(genes.keys())), desc="Building gene/transcript reference.") if verbose else genes:
        geneObj = genes[gene]
        for txn,tObj in geneObj.transcripts.items():
            juncs = [(tObj.exons[pos].c2,tObj.exons[pos+1].c1-1) for pos,e in enumerate(tObj.exons,0) if pos+1<len(tObj.exons)]
            leftJuncs = [x[0] for x in juncs]
            rightJuncs = [x[1] for x in juncs]
            juncs = tuple(sorted(leftJuncs+rightJuncs))

            geneRef[geneObj.chrom][juncs] = (geneObj,tObj)
            for junc in juncs:
                juncRef[geneObj.chrom][junc] = (geneObj,tObj)
            

    return geneRef, juncRef

def main():
    '''
    not the state.
    '''
    myCommandLine = CommandLine()
    bed     = myCommandLine.args['input_isoforms']
    gtf     = myCommandLine.args['gtf']
    
    global verbose
    verbose = myCommandLine.args['quiet']

    geneRef, juncRef = buildGTFDB(gtf)

    bedObj = BED12(bed)
    for bedCols in bedObj.getLine():
        juncs = bedObj.bed12toJuncs()
        chrom = bedCols[0]
        leftJuncs = [x[0] for x in juncs]
        rightJuncs = [x[1] for x in juncs]
        juncs = tuple(sorted(leftJuncs+rightJuncs))


        if juncs in geneRef[chrom]:
            # do stuff
            gid, tid = geneRef[chrom][juncs]
            gid, tid = gid.gid, tid.tid

        else:
            newSet = set([None])
            for i in juncs: 
                try:
                    newSet.add(juncRef[chrom][i][0].gid)
                except:
                    pass
                    #No gene

            newSet.remove(None)

            if len(newSet) < 1:
                # new "gene"
                gid = bedCols[3]
                tid = bedCols[3]
                
            else:
                gid = ",".join([x for x in list(newSet)])
                # Weird naming bug ...get better scheme. This code will fix for now.
                gids = set(gid.split(","))
                gid = ",".join([x for x in list(gids)])
                tid = bedCols[3]

            # now add to generef

            newGObj, newTObj = Gene(gid), Transcript(tid) 
            geneRef[juncs] = (newGObj, newTObj)
            for x in juncs: juncRef[x] = (newGObj, newTObj)

        bedName = "%s_%s" % (gid, tid)
        bedCols[3] = bedName
        print("\t".join(map(str, bedCols)))

if __name__ == "__main__":
    main()
