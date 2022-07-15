#!/usr/bin/env python3
from __future__ import print_function

########################################################################
# File: runFish.py
#  executable: runFish.py
# Purpose: 
#
#          
# Author: Cameron M. Soulette
# History:      cms 03/25/2019 Created
#
########################################################################


########################################################################
# Hot Imports & Global Variable
########################################################################

import os, sys, csv
from itertools import combinations
import scipy.stats as sps
import statsmodels.stats.multitest as sm
import numpy as np
import re

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
        self.parser = argparse.ArgumentParser(description = ' fish_suppaEvents.py ',
                                             epilog = '', 
                                             add_help = True, #default is True 
                                             prefix_chars = '-', 
                                             usage = '%(prog)s ')
        # Add args
        
        self.parser.add_argument("--count_matrix"    , action = 'store', required=True, 
                                    help='Count Matrix')
        self.parser.add_argument("--ioe_table"    , action = 'store', required=True, 
                                    help='ioe table from suppa2')
        self.parser.add_argument("--out_prefix"     , action = 'store', required=True,
                                    help='Output filename prefix.')

        if inOpts is None :
            self.args = vars(self.parser.parse_args())
        else :
            self.args = vars(self.parser.parse_args(inOpts))

# funktions

def getQuant(x):
    #takes in count matrix formatted tsv file, returns
    #dictionary where gene/iso id is key and count is value
    # also returns list of samples for fishers test
    gQ, iQ = dict(),dict()
    with open(x) as lines:
        samples = {pos: sample for pos, sample in enumerate(next(lines).rstrip().split()[1:],0)}
        for i in lines:
            cols = i.rstrip().split()
            vals = np.asarray(cols[1:], dtype=float)
            gid = re.search("(ENSG[^\.]+|chr[^\-]+)",cols[0]).group(1)
            tid = cols[0].split("_")[0]

            if gid not in gQ:
                gQ[gid] = np.zeros(len(samples))
            gQ[gid] = gQ[gid] + vals
            iQ[tid] = vals

    return gQ, iQ, samples


def fishAS(ioe, gQ, iQ, s, outPrefix):
    '''
    stuff
    '''

    pvalues = dict()
    events = list()
    sampleCombos = [x for x in combinations(list(s.keys()),2)] 
    numSamples = len(list(s.keys()))
    print(sampleCombos)
    with open(ioe) as lines:
        next(lines)
        for line in lines:
            c = line.rstrip().split()
            gene, event, inclusionTxn, allTxn  = c[1], c[2], set(c[3].split(",")), set(c[4].split(","))
            exclusionTxn = allTxn.difference(inclusionTxn)

            if len(exclusionTxn)<1 or len(inclusionTxn)<1:
                # weird bug for now.
                # has to to with isoform TES/TSS
                continue

            
            inclusionVal = np.asarray([iQ.get(x,np.zeros(numSamples)) for x in inclusionTxn]).sum(axis=0)
            exclusionVal = np.asarray([iQ.get(x,np.zeros(numSamples)) for x in exclusionTxn]).sum(axis=0)

            sumCount = inclusionVal + exclusionVal
            psi = inclusionVal / sumCount

            for combo in sampleCombos:
                i,j = combo
                key = "%s_v_%s" % (s[i],s[j])
                if key not in pvalues:
                    pvalues[key] = list()

                exclVal1, exclVal2 = exclusionVal[i], exclusionVal[j]
                inclVal1, inclVal2 = inclusionVal[i], inclusionVal[j]
                
                if exclVal1 + inclVal1 < 35 or exclVal2 + inclVal2 < 35:
                    continue
                if exclVal1<25 and exclVal2<25:
                    continue
                if inclVal1<25 and inclVal2<25:
                    continue

                table = [[inclVal1, inclVal2],[exclVal1, exclVal2]]
                pval  = sps.fisher_exact(table)[1]
                pvalues[key].append((event,inclVal1, exclVal1+inclVal1, inclVal2, exclVal2+inclVal2, pval))
    return pvalues
# main

def main():
    '''
    maine
    '''

    # Command Line Stuff...
    myCommandLine = CommandLine()

    matrix     = myCommandLine.args['count_matrix']
    ioeTable   = myCommandLine.args['ioe_table']
    outPrefix    = myCommandLine.args['out_prefix']
    
    geneQuants,isoQuants, samples = getQuant(matrix)
    
    if os.path.isfile(ioeTable):
        pvalues = fishAS(ioeTable,geneQuants,isoQuants,samples, outPrefix)
    else:
        isoUsage(geneQuants,isoQuants,samples) 

    for sample, plist in pvalues.items():
        corrected = sm.multipletests([x[-1] for x in plist], method='hs')[1]
        with open("%s_%s_fishersP.tsv" % (outPrefix,sample), 'w') as out1:
            print("event","samp1_inclusion", "samp1_total", "samp2_inclusion", "sampe2_total", "pvalue", "bh.adj.pval", sep="\t", file=out1)
            [print("\t".join(map(str,x)),corrected[num],sep="\t", file=out1) for num,x in enumerate(plist)]


if __name__ == "__main__":
    main()
