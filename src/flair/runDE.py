#!/usr/bin/env python3
from __future__ import print_function

########################################################################
# File: runDE.py
#  executable: runDE.py
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
        self.parser = argparse.ArgumentParser(description = ' runDE.py - a rpy2 convenience tool to run DESeq2.',
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

    print("running DESEQ2 %s" % prefix, file=sys.stderr)

    # make the quant DF
    quantDF  = pd.read_csv(matrix, header=0, sep='\t', index_col=0)
    df = pandas2ri.py2ri(quantDF)

    # import formula
    formulaDF     = pd.read_csv(formula,header=0, sep="\t",index_col=0)
    sampleTable = pandas2ri.py2ri(formulaDF)


    if "batch" in list(formulaDF):
        design = Formula("~ condition + condition")
    else:
        design = Formula("~ condition")
   

    # import DESeq2
    from rpy2.robjects.packages import importr
    import rpy2.robjects.lib.ggplot2 as ggplot2
    methods   = importr('methods')
    deseq     = importr('DESeq2')
    grdevices = importr('grDevices')
    qqman     = importr('qqman')



    ### RUN DESEQ2 ###
    R.assign('df', df)
    R.assign('sampleTable', sampleTable)
    R.assign('design',design)
    R('dds <- DESeqDataSetFromMatrix(countData = df, colData = sampleTable, design = design)')
    R('dds <- DESeq(dds)')
    # R("save.image(file='/private/groups/brookslab/atang/flair/testing/misc/colette/.RData')")
    R('f = resultsNames(dds)')
    g = robjects.r['f']
    print(g)
    R('name <- grep("condition", resultsNames(dds), value=TRUE)')

    ###
    ###
    # Get Results and shrinkage values
    res    = R('results(dds, name=name)')
    resLFC = R('lfcShrink(dds, coef=name)')
    vsd    = R('vst(dds,blind=FALSE)')
    resdf  = robjects.r['as.data.frame'](res) 
    reslfc = robjects.r['as.data.frame'](resLFC)
    dds    = R('dds')

    
    ### Plotting section ###
    # plot MA and PC stats for the user
    plotMA    = robjects.r['plotMA']
    plotDisp  = robjects.r['plotDispEsts']
    plotPCA   = robjects.r['plotPCA']
    plotQQ    = robjects.r['qq']
    
    # get pca data
    if "batch" in list(formulaDF):
        pcaData    = plotPCA(vsd, intgroup=robjects.StrVector(("condition", "batch")), returnData=robjects.r['T'])
        percentVar = robjects.r['attr'](pcaData, "percentVar")
    else:
        print(vsd)
        pcaData    = plotPCA(vsd, intgroup="condition", returnData=robjects.r['T'])
        percentVar = robjects.r['attr'](pcaData, "percentVar")
    # arrange 


    data_folder = os.path.join(os.getcwd(), outdir)
    qcOut = os.path.join(data_folder, "%s_QCplots_%s_v_%s.pdf"  % (prefix,group1,group2))
    
    grdevices.pdf(file=qcOut)

    x = "PC1: %s" % int(percentVar[0]*100) + "%% variance"
    y = "PC2: %s" % int(percentVar[1]*100) + "%% variance"

    if "batch" in list(formulaDF):
        pp = ggplot2.ggplot(pcaData) + \
                ggplot2.aes_string(x="PC1", y="PC2", color="condition", shape="batch") + \
                ggplot2.geom_point(size=3) + \
                robjects.r['xlab'](x) + \
                robjects.r['ylab'](y) + \
                ggplot2.theme_classic() + \
                ggplot2.coord_fixed()

    else:
        pp = ggplot2.ggplot(pcaData) + \
                ggplot2.aes_string(x="PC1", y="PC2", color="condition") + \
                ggplot2.geom_point(size=3) + \
                robjects.r['xlab'](x) + \
                robjects.r['ylab'](y) + \
                ggplot2.theme_classic() + \
                ggplot2.coord_fixed()
    pp.plot()
    plotMA(res, ylim=robjects.IntVector((-3,3)), main="MA-plot results")
    plotMA(resLFC, ylim=robjects.IntVector((-3,3)), main="MA-plot LFCSrhinkage")    
    plotQQ(reslfc.rx2('pvalue'), main="LFCSrhinkage pvalue QQ")
    hh = ggplot2.ggplot(resdf) + \
            ggplot2.aes_string(x="pvalue") + \
            ggplot2.geom_histogram() + \
            ggplot2.theme_classic() + \
            ggplot2.ggtitle("pvalue distribution")
    hh.plot()
    plotDisp(dds, main="Dispersion Estimates")
    grdevices.dev_off()


    data_folder = os.path.join(os.getcwd(), outdir)
    lfcOut = os.path.join(data_folder, "%s_%s_v_%s_deseq2_results_shrinkage.tsv"  % (prefix,group1,group2))
    resOut = os.path.join(data_folder, "%s_%s_v_%s_deseq2_results.tsv"  % (prefix,group1,group2))
   
    robjects.r['write.table'](reslfc, file=lfcOut, quote=False, sep="\t")
    robjects.r['write.table'](resdf, file=resOut, quote=False, sep="\t")


if __name__ == "__main__":
    main()
