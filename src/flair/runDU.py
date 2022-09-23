#!/usr/bin/env python3

########################################################################
# File: runDU.py
#  executable: runDU.py
# Purpose:
#
#
# Author: Cameron M. Soulette
# History:      cms 12/05/2018 Created
#
########################################################################


########################################################################
# Hot Imports & Global Variable
########################################################################


import os
import sys
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import pandas as pd
import numpy as np
import argparse

import rpy2
from rpy2 import robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
#    import rpy2.robjects.lib.ggplot2 as ggplot2
pandas2ri.activate()
R = robjects.r

import warnings
from rpy2.rinterface import RRuntimeWarning
warnings.filterwarnings("ignore", category=RRuntimeWarning)

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
        self.parser = argparse.ArgumentParser(description=' runDU.py - a rpy2 convenience tool to run DRIMseq.',
                                             add_help=True, #default is True
                                             prefix_chars='-',
                                             usage='%(prog)s ')
        # Add args
        self.parser.add_argument("--group1", action='store', required=True,
                                    help='Sample group 1.')
        self.parser.add_argument("--group2", action='store', required=True,
                                    help='Sample group 2.')
        self.parser.add_argument("--batch", action='store', required=False, default=None,
                                    help='Secondary sample attribute (used in design matrix).')
        self.parser.add_argument("--matrix", action='store', required=True,
                                    help='Input count files.')
        self.parser.add_argument("--outDir", action='store', required=True,
                                    help='Write to specified output directory.')
        self.parser.add_argument("--prefix", action='store', required=True,
                                    help='Specify file prefix.')
        self.parser.add_argument("--formula", action='store', required=True,
                                    help='Formula design matrix.')
        self.parser.add_argument("--threads", type=int, action='store',default=4, required=False,
                                    help='Number of threads for running DRIM-Seq. BBPARAM')

        if inOpts is None:
            self.args = vars(self.parser.parse_args())
        else:
            self.args = vars(self.parser.parse_args(inOpts))


# main
def main():
    '''
    '''

    # Command Line Stuff...
    myCommandLine = CommandLine()

    outdir     = myCommandLine.args['outDir']
    group1     = myCommandLine.args['group1']
    group2     = myCommandLine.args['group2']
#    batch      = myCommandLine.args['batch'] unused, check args
    matrix     = myCommandLine.args['matrix']
    prefix     = myCommandLine.args['prefix']
    formula    = myCommandLine.args['formula']
    threads    = myCommandLine.args['threads']

    print("running DRIMSEQ %s" % prefix, file=sys.stderr)

    rundrimseq(outdir, group1, group2, matrix, prefix, formula, threads)


def rundrimseq(outdir, group1, group2, matrix, prefix, formula, threads):
    '''Run DRIMSeq via rpy R emulator'''

    print(f'input file: {matrix}', file=sys.stderr)

    # create output working directory if it doesn't exist
    data_folder = os.path.join(os.getcwd(), outdir)
    workdir = os.path.join(data_folder, 'workdir')
    if not os.path.exists(workdir):
        os.makedirs(workdir)
    resOut = os.path.join(workdir, "%s_%s_v_%s_results.tsv"  % (prefix,group1,group2))
    # final file goes in the main output dir
    cleanOut = os.path.join(data_folder, "%s_%s_v_%s.tsv"  % (prefix,group1,group2))

    # clean up rpy2/R's stderr
    def f(x):
        print(x.rstrip(), file=sys.stderr)
    rpy2.rinterface_lib.callbacks.consolewrite_warnerror = f

    # import
    importr('methods')
    importr('DRIMSeq')

    # get quant table and formula table
    quantDF  = pd.read_csv(matrix, header=0, sep='\t', index_col=0)
    df       = pandas2ri.py2rpy(quantDF)

    formulaDF = pd.read_csv(formula,header=0, sep="\t")

    pydf      = pandas2ri.py2rpy(formulaDF)

    # Convert pandas to R data frame.
    samples = pydf
    counts  = df

    # DRIMSEQ part.
    # Formula
    if "batch" in list(formulaDF): R.assign('batch', samples.rx2('batch'))
    R.assign('condition', samples.rx2('condition'))
    R.assign('counts', counts)
    R.assign('samples',samples)
    R.assign('numThread', threads)
    R.assign("cooef", "condition%s" % group2)

    R('data <- dmDSdata(counts = counts, samples = samples)')
    R('filtered <- dmFilter(data, min_samps_gene_expr = 6, min_samps_feature_expr = 3, min_gene_expr = 15, min_feature_expr = 5)')
    if "batch" in list(formulaDF):
        R('design_full <- model.matrix(~ condition + batch, data = samples(filtered))')
    else:
        R('design_full <- model.matrix(~ condition, data = samples(filtered))')
    R('set.seed(123)')

    R('d <- dmPrecision(filtered, design = design_full, BPPARAM=BiocParallel::MulticoreParam(numThread))')
    R('d <- dmFit(d, design = design_full, verbose = 1, BPPARAM=BiocParallel::MulticoreParam(numThread))')

    #print(rpy2.__version__)
    #print(np.__version__)
    R('contrast = colnames(design_full)[2]')

    R('d <- dmTest(d, coef = contrast, verbose = 1, BPPARAM=BiocParallel::MulticoreParam(numThread))')
    #res = R('merge(proportions(d),results(d,level="feature"), by=c("feature_id","gene_id"))')
    #res.to_csv(resOut, sep='\t')
    R('res <- merge(proportions(d),results(d,level="feature"), by=c("feature_id","gene_id"))')
    # raw output
    res = R('res')
    res.to_csv(resOut, sep='\t')

    # order by adjusted p value
    R('res <- res[order(res[,"adj_pvalue"]),]')
    # keep only significant
    R('res <- subset(res, adj_pvalue <= 0.05)')
    # we don't know how many samples there are
    R('fina <- ncol(res)-4')
    # round sample and log ratio value, skip df and pvalue, round(=signif) the e-values of adj_pvalue
    R('outdf <- data.frame(res[1:2], round(res[3:fina],3), lr = round(res[["lr"]],2), adj_pvalue = signif(res[["adj_pvalue"]],3))')
    R.assign('outf', cleanOut)
    R('write.table(outdf, row.names=FALSE, quote=FALSE, file=outf, sep="\t")')

if __name__ == "__main__":
    main()
