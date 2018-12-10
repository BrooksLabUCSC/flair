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
                                             epilog = 'Please feel free to forward any questions/concerns to /dev/null', 
                                             add_help = True, #default is True 
                                             prefix_chars = '-', 
                                             usage = '%(prog)s --workingdir dir_name --outdir out_dir --group1 conditionA --group2 conditionB --batch batch --files conditionA_batch1.txnCounts,conditionB_batch1.txnCounts')
        # Add args
        self.parser.add_argument("--workingdir", action = 'store', required=True,
                                    help='Input dir containing files.')
        self.parser.add_argument("--outdir"    , action = 'store', required=True, 
                                    help='Output dir to write files.')
        self.parser.add_argument("--filter"    , action = 'store', required=False, default = 20, type=int,
                                    help='Output file name prefix.')
        self.parser.add_argument("--group1"    , action = 'store', required=True, 
                                    help='Sample group 1.')
        self.parser.add_argument("--group2"    , action = 'store', required=True, 
                                    help='Sample group 2.')
        self.parser.add_argument("--batch"     , action = 'store', required=False, default=None,
                                    help='Secondary sample attribute (used in design matrix).')
        self.parser.add_argument("--files"     , action = 'store', required=True,
                                    help='Input count files.')
        self.parser.add_argument("--out_prefix", action = 'store', required=True,
                                    help='Output file name prefix.')
       


        if inOpts is None :
            self.args = vars(self.parser.parse_args())
        else :
            self.args = vars(self.parser.parse_args(inOpts))




def makeDir(out):
    try:
        os.mkdir("./%s" % out)
    except:
        #exists
        pass

def checkSamples(samples):
    samples = samples.split(",")
    if len(samples)<3:
        print("Error. At least 3 samples for each condition required. Exit.",file=sys.stderr)
        sys.exit(1)
    return samples


def filesToDF(f, thresh):
    
    data = dict()
    for num,i in enumerate(f,0):
        with open(i) as l:
            for line in l:
                name, count = line.rstrip().split()
                count = int(count) 
                if name not in data:
                    data[name] = np.zeros(len(f))
                data[name][num] += count

    indices         = np.asarray(list(data.keys()))
    values          = np.asarray([data[x] for x in indices], dtype=int)
    filtered        = values[np.min(values, axis=1) > thresh]
    filteredIndices = indices[np.min(values, axis=1) > thresh]
    print(len(filtered),len(values),len(filteredIndices))
    df = pd.DataFrame(filtered,columns=f)
    df['ids'] = filteredIndices
    df = df.set_index('ids')
    return df

# main



def main():
    '''
    maine
    '''

    # Command Line Stuff...
    myCommandLine = CommandLine()

    workingdir = myCommandLine.args['workingdir']
    outdir     = myCommandLine.args['outdir']
    group1     = myCommandLine.args['group1']
    group2     = myCommandLine.args['group2']
    batch      = myCommandLine.args['batch']  
    files      = myCommandLine.args['files']
    prefix     = myCommandLine.args['out_prefix']
    sFilter    = myCommandLine.args['filter']

    makeDir(outdir)

    files = checkSamples(files)

    df = filesToDF(files, sFilter)
      
    # DO DESEQ2
    from rpy2 import robjects
    from rpy2.robjects import r,pandas2ri, Formula
    from rpy2.robjects.lib import grid
    pandas2ri.activate()

    # Compile data for data frame
    data = list()
    for f in files:
        if group1 in f:
            if batch in f:
                data.append((f,group1,'1'))
            else:
                data.append((f,group1,'2'))
        else:
            if batch in f:
                data.append((f,group2,'1'))
            else:
                data.append((f,group2,'2'))

    # Make the Data Frame
    pydf = pd.DataFrame(data)
    pydf.columns = ['sampleName','condition','batch']
    pydf = pydf.set_index('sampleName')
    # Convert pandas to R data frame.
    sampleTable = pandas2ri.py2ri(pydf)

    # DESEQ2 part.
    # Forumla
    design = Formula("~ batch + condition")

    # import DESeq2
    from rpy2.robjects.packages import importr
    import rpy2.robjects.lib.ggplot2 as ggplot2
    methods   = importr('methods')
    deseq     = importr('DESeq2')
    grdevices = importr('grDevices')
    qqman     = importr('qqman')

    # dds = deseq.DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
    #                                         directory = workingdir,
    #                                         design= design)

    dds = deseq.DESeqDataSetFromMatrix(countData = df,
                                        colData = sampleTable,
                                        design = design)
    dds = deseq.DESeq(dds)
    
    # get results; orient the results for groupA vs B
    res = deseq.results(dds, contrast=robjects.StrVector(("condition",group2,group1)))
    # results with shrinkage
    resLFC = deseq.lfcShrink(dds, coef="condition_%s_vs_%s" % (group2,group1), type="apeglm")
    resdf  = robjects.r['as.data.frame'](res)
    reslfc  = robjects.r['as.data.frame'](resLFC)

    # plot MA and PC stats for the user
    plotMA    = robjects.r['plotMA']
    plotDisp  = robjects.r['plotDispEsts']
    plotPCA   = robjects.r['plotPCA']
    plotQQ    = robjects.r['qq']
    
    vsd       = robjects.r['vst'](dds, blind=robjects.r['F'])
    # get pca data
    pcaData    = plotPCA(vsd, intgroup=robjects.StrVector(("condition", "batch")), returnData=robjects.r['T'])
    percentVar = robjects.r['attr'](pcaData, "percentVar")

    # arrange 
    grdevices.pdf(file="./%s/%s_%s_vs_%s_%s_%s_cutoff_plots.pdf" % (outdir,prefix,group1,group2,str(batch),sFilter))


    x = "PC1: %s" % int(percentVar[0]*100) + "%% variance"
    y = "PC2: %s" % int(percentVar[1]*100) + "%% variance"
    
    pp = ggplot2.ggplot(pcaData) + \
            ggplot2.aes_string(x="PC1", y="PC2", color="condition", shape="batch") + \
            ggplot2.geom_point(size=3) + \
            robjects.r['xlab'](x) + \
            robjects.r['ylab'](y) + \
            ggplot2.theme_classic() + \
            ggplot2.coord_fixed()
    pp.plot()

    plotMA(res, ylim=robjects.IntVector((-3,3)), main="MA-plot results")
    #plotMA(res, main="MA-plot results")
    plotMA(resLFC, ylim=robjects.IntVector((-3,3)), main="MA-plot LFCSrrhinkage")
    #plotMA(resLFC, main="MA-plot LFCSrrhinkage")
    plotQQ(resdf.rx2('pvalue'), main="pvalue QQ")
    plotQQ(reslfc.rx2('pvalue'), main="LFCSrhinkage pvalue QQ")
    hh = ggplot2.ggplot(resdf) + \
            ggplot2.aes_string(x="pvalue") + \
            ggplot2.geom_histogram() + \
            ggplot2.theme_classic() 
    hh.plot()
    plotDisp(dds, main="Dispersion Estimates")
    grdevices.dev_off()


    reslsf = pandas2ri.ri2py(reslfc)
    res    = pandas2ri.ri2py(resdf)

    reslsf.to_csv("./%s/%s_%s_vs_%s_%s_deseq2_results_LFC.tsv" % (outdir,prefix,group1,group2,str(batch)), 
                    sep='\t')
    reslsf.to_csv("./%s/%s_%s_vs_%s_%s_deseq2_results.tsv" % (outdir,prefix,group1,group2,str(batch)), 
                    sep='\t')

if __name__ == "__main__":
    main()