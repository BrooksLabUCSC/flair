#!/usr/bin/env python3
from __future__ import print_function

########################################################################
# File: runDS.py
#  executable: runDS.py
# Purpose: 
#
#          
# Author: Cameron M. Soulette + Alison Tang
# History: 01/31/2020 Created
#
########################################################################


########################################################################
# Hot Imports & Global Variable
########################################################################

import warnings
from rpy2.rinterface import RRuntimeWarning
warnings.filterwarnings("ignore", category=RRuntimeWarning)

import os, sys
import pandas as pd
import numpy as np

from rpy2 import robjects

import rpy2.robjects as ro 
from rpy2.robjects import r, pandas2ri, Formula
from rpy2.robjects.conversion import localconverter


from rpy2.robjects.lib import grid
R = robjects.r


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
        self.parser = argparse.ArgumentParser(description = 'runDS.py - a rpy2 convenience tool to run DRIMseq.',
                                             add_help = True, #default is True 
                                             prefix_chars = '-', 
                                             usage = '%(prog)s ')
        # Add args
        self.parser.add_argument("--matrix"    , action='store', required=True,
                                    help='Input DRIM-Seq formatted count files.')
        self.parser.add_argument("--outDir"    , action='store', required=False, 
                                    help='Write to specified output directory.', default='')
        self.parser.add_argument("--prefix"    , action='store', required=True, 
                                    help='Specify file prefix.')
        self.parser.add_argument('--drim1'     , action='store', type=int, required=False, default=6, \
            help='''The minimum number of samples that have coverage over an AS event inclusion/exclusion;
                                events with too few samples are filtered out and not tested (6)''')
        self.parser.add_argument('--drim2'     , action='store', type=int, required=False, default=3, \
            help='''The minimum number of samples expressing the inclusion of an AS event; events with
                                too few samples are filtered out and not tested (3)''')
        self.parser.add_argument('--drim3'     , action='store', type=int, required=False, default=15, \
            help='''The minimum number of reads covering an AS event inclusion/exclusion, events with
                                too few samples are filtered out and not tested (15)''')
        self.parser.add_argument('--drim4'     , action='store', type=int, required=False, default=5, \
            help='''The minimum number of reads covering an AS event inclusion, events with too few
                                samples are filtered out and not tested (5)''')
        self.parser.add_argument("--threads"   , action='store', type=int, default=4, required=False, 
                                    help='Number of threads for running DRIM-Seq. BBPARAM')
        self.parser.add_argument('--batch', action='store_true', dest='batch', required=False, default=False, \
                                    help='''If specified, batch correction will be performed''')
        self.parser.add_argument('--conditionA', action='store', dest='conditionA', required=False, default='', \
            help='''Specify one condition corresponding to samples in the counts_matrix to be compared against
                                condition2; by default, the first two unique conditions are used''')
        self.parser.add_argument('--conditionB', action='store', dest='conditionB', required=False, default='', \
            help='''Specify one condition corresponding to samples in the counts_matrix to be compared against
                                 condition1''')

        if inOpts is None :
            self.args = vars(self.parser.parse_args())
        else :
            self.args = vars(self.parser.parse_args(inOpts))


# main
def main():
    '''
    maine
    '''

    # Command Line Stuff...
    myCommandLine = CommandLine()

    outdir     = myCommandLine.args['outDir']
    matrix     = myCommandLine.args['matrix']
    prefix     = myCommandLine.args['prefix']
    threads    = myCommandLine.args['threads']
    drim1      = myCommandLine.args['drim1']
    drim2      = myCommandLine.args['drim2']
    drim3      = myCommandLine.args['drim3']
    drim4      = myCommandLine.args['drim4']
    usebatch   = myCommandLine.args['batch']
    conditionA = myCommandLine.args['conditionA']
    conditionB = myCommandLine.args['conditionB']

    samples = open(matrix).readline().rstrip().split('\t')[2:-1]
    groups  = [x.split("_")[1] for x in samples]
    batches = [x.split("_")[-1] for x in samples]

    if not conditionA:  # determine the first two unique conditions
        conditionA = groups[0]
        conditionB = [g for g in groups if g != conditionA][0]

    # import
    from rpy2.robjects.packages import importr
    import rpy2.robjects.lib.ggplot2 as ggplot2
    methods   = importr('methods')
    drim      = importr('DRIMSeq')

    # make formula (sample) dataframe
    formula = []
    for s,g,b in zip(samples, groups, batches):
        if conditionA not in g and conditionB not in g:
            continue
        formula += [(s, g, b)]
    formulaDF = pd.DataFrame(data=formula, columns = ['sample_id', 'condition', 'batch'])

    # get quant table
    quantDF  = pd.read_csv(matrix, header=0, sep='\t', index_col=False)

    # subset the counts_matrix for only the samples with condA/B
    quantDF  = quantDF[['feature_id', 'coordinate'] + list(formulaDF['sample_id']) + ['isoform_ids']]
    # rename coordinate column name to gene_id
    quantDF.columns = ['feature_id', 'gene_id'] + list(quantDF.columns[2:])

    # add a pseudocount of 1 to each event in each sample
    for col in list(formulaDF['sample_id']):
        quantDF[[col]] = quantDF[[col]] + 1

    # Convert pandas to R data frame.
    import rpy2
    rpy2_version = rpy2.__version__
    rpy2_version = float(rpy2_version[:rpy2_version.rfind('.')])
    if rpy2_version >= 3.4:
        with localconverter(ro.default_converter + pandas2ri.converter):
            samples = ro.conversion.py2rpy(formulaDF)
            counts = ro.conversion.py2rpy(quantDF)
    else:
        pandas2ri.activate()
        samples = pandas2ri.py2ri(formulaDF)
        counts = pandas2ri.py2ri(quantDF)

    # DRIMSEQ part.
    if "batch" in list(formulaDF): R.assign('batch', samples.rx2('batch'))
    R.assign('condition', samples.rx2('condition'))
    R.assign('counts', counts)
    R.assign('samples', samples)
    R.assign('numThread', threads)
    R.assign('cooef', "condition%s" % conditionA)
    R.assign('drim1', drim1)
    R.assign('drim2', drim2)
    R.assign('drim3', drim3)
    R.assign('drim4', drim4)

    R('data <- dmDSdata(counts = counts, samples = samples)')
    R('filtered <- dmFilter(data, min_samps_gene_expr = drim1, min_samps_feature_expr = drim2, min_gene_expr = drim3, min_feature_expr = drim4)')
    if usebatch:
        R('design_full <- model.matrix(~ condition + batch, data = samples(filtered))')
    else: 
        R('design_full <- model.matrix(~ condition, data = samples(filtered))')

    R('set.seed(123)')

    R('d <- dmPrecision(filtered, design = design_full, BPPARAM=BiocParallel::MulticoreParam(numThread))')
    R('d <- dmFit(d, design = design_full, verbose = 1, BPPARAM=BiocParallel::MulticoreParam(numThread))')

    R('contrast <- grep("condition",colnames(design_full),value=TRUE)')
    R('d <- dmTest(d, coef = contrast, verbose = 1, BPPARAM=BiocParallel::MulticoreParam(numThread))')


    res = R('merge(proportions(d),results(d,level="feature"), by=c("feature_id","gene_id"))')
    resOut = "%s.%s_v_%s_drimseq_results.tsv"  % (prefix, conditionA, conditionB)
    if rpy2_version >= 3.4:
        with localconverter(ro.default_converter + pandas2ri.converter):
            res = ro.conversion.rpy2py(res)

    res.to_csv(resOut, sep='\t')
    # pltFName = '%s_%s_v_%s_pval_histogram.pdf' % (prefix,conditionA,conditionB)
    # R.assign('fname', pltFName)
    # R('pdf(file=fname)')
    # R('plotPValues(res)')  # histogram of p-value distribution
    # R('dev.off()')



if __name__ == "__main__":
    main()
