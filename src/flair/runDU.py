#!/usr/bin/env python3
from __future__ import print_function

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


import os, sys
import pandas as pd
import numpy as np

from rpy2 import robjects
from rpy2.robjects import r,pandas2ri, Formula
from rpy2.robjects.lib import grid
pandas2ri.activate()
R = robjects.r

import warnings
from rpy2.rinterface import RRuntimeWarning
warnings.filterwarnings("ignore", category=RRuntimeWarning)

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
        self.parser = argparse.ArgumentParser(description = ' runDU.py - a rpy2 convenience tool to run DRIMseq.',
                                             add_help = True, #default is True
                                             prefix_chars = '-',
                                             usage = '%(prog)s ')
        # Add args
        self.parser.add_argument("--group1"    , action = 'store', required=True,
                                    help='Sample group 1.')
        self.parser.add_argument("--group2"    , action = 'store', required=True,
                                    help='Sample group 2.')
        self.parser.add_argument("--batch"     , action = 'store', required=False, default=None,
                                    help='Secondary sample attribute (used in design matrix).')
        self.parser.add_argument("--matrix"     , action = 'store', required=True,
                                    help='Input count files.')
        self.parser.add_argument("--outDir"    , action = 'store', required=True,
                                    help='Write to specified output directory.')
        self.parser.add_argument("--prefix"    , action = 'store', required=True,
                                    help='Specify file prefix.')
        self.parser.add_argument("--formula"    , action = 'store', required=True,
                                    help='Formula design matrix.')
        self.parser.add_argument("--threads"    , type=int, action = 'store',default=4, required=False,
                                    help='Number of threads for running DRIM-Seq. BBPARAM')


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
    group1     = myCommandLine.args['group1']
    group2     = myCommandLine.args['group2']
    batch      = myCommandLine.args['batch']
    matrix     = myCommandLine.args['matrix']
    prefix     = myCommandLine.args['prefix']
    formula    = myCommandLine.args['formula']
    threads    = myCommandLine.args['threads']

    print("running DRIMSEQ %s" % prefix, file=sys.stderr)

    # import
    from rpy2.robjects.packages import importr
    import rpy2.robjects.lib.ggplot2 as ggplot2
    methods   = importr('methods')
    drim      = importr('DRIMSeq')

    # get quant table and formula table
    quantDF  = pd.read_csv(matrix, header=0, sep='\t', index_col=0)
    df       = pandas2ri.py2ri(quantDF)

    formulaDF = pd.read_csv(formula,header=0, sep="\t")

    pydf      = pandas2ri.py2ri(formulaDF)

    # Convert pandas to R data frame.
    samples = pydf
    counts  = df

    # DRIMSEQ part.
    # Forumla
    print('here')
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

    # R('f = colnames(design_full)')
    # R("save.image(file='/private/groups/brookslab/atang/flair/testing/misc/colette/.RData')")
    # f = robjects.r['f']
    # print('here2')
    # g = robjects.r['condition']
    # print(f)
    # print(g)
    import rpy2
    print(rpy2.__version__)
    print(np.__version__)
    # R('contrast <- grep("condition",colnames(design_full),value=TRUE)')
    R('contrast = colnames(design_full)[2]')

    R('d <- dmTest(d, coef = contrast, verbose = 1, BPPARAM=BiocParallel::MulticoreParam(numThread))')
    res = R('merge(proportions(d),results(d,level="feature"), by=c("feature_id","gene_id"))')

    data_folder = os.path.join(os.getcwd(), outdir)
    resOut = os.path.join(data_folder, "%s_%s_v_%s_drimseq2_results.tsv"  % (prefix,group1,group2))


    res.to_csv(resOut, sep='\t')
    sys.exit(1)

    R('library(stageR)')
    R('pScreen <- results(d)$pvalue')
    R('names(pScreen) <- results(d)$gene_id')
    ## Assign transcript-level pvalues to the confirmation stage
    R('pConfirmation <- matrix(results(d, level = "feature")$pvalue, ncol = 1)')
    R('rownames(pConfirmation) <- results(d, level = "feature")$feature_id')
    ## Create the gene-transcript mapping
    R('tx2gene <- results(d, level = "feature")[, c("feature_id", "gene_id")]')
    ## Create the stageRTx object and perform the stage-wise analysis
    R('stageRObj <- stageRTx(pScreen = pScreen, pConfirmation = pConfirmation, pScreenAdjusted = FALSE, tx2gene = tx2gene)')
    R('stageRObj <- stageWiseAdjustment(object = stageRObj, method = "dtu", alpha = 0.05)')
    R('getSignificantGenes(stageRObj)')
    R('getSignificantTx(stageRObj)')
    R('padj <- getAdjustedPValues(stageRObj, order = TRUE, onlySignificantGenes = FALSE)')
    R('head(padj)')
    #     # import plotting
    # plotMA    = robjects.r['plotData']
    # plotPrec  = robjects.r['plotPrecision']
    # plotQQ    = robjects.r['qq']

    # # arrange
    # pltFName = './%s/%s_%s_vs_%s_%s_%s_cutoff_plots.pdf' % (outdir,prefix,group1,group2,str(batch),sFilter)
    # R.assign('fname',pltFName)
    # R('pdf(file=fname)')
    # R('plotData(filtered)')
    # R('ggp <- plotPrecision(d)')
    # R('ggp + geom_point(size = 4, alpha=0.3)')
    # R('plotPValues(d)')
    # R('dev.off()')
    # grdevices.pdf(file="./%s/%s_%s_vs_%s_%s_%s_cutoff_plots.pdf" % (outdir,prefix,group1,group2,str(batch),sFilter))
    # qqman(res['pvalue'])
    # grdevices.dev_off()


if __name__ == "__main__":
    main()
