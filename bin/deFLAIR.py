from __future__ import print_function

########################################################################
# File: diaFLAIR.py
#  executable: diaFLAIR.py
# Purpose: wrapper for Differential Isoform Analyses
#
#          
# Author: Cameron M. Soulette
# History:      cms 01/17/2019 Created
#
########################################################################


########################################################################
# Hot Imports & Global Variable
########################################################################


import os, sys
import pandas as pd
import numpy as np
import subprocess 

scriptPath = os.path.realpath(__file__)
path = "/".join(scriptPath.split("/")[:-1])
runDE = path + "/" + "runDE.py"
runDU = path + "/" + "runDU.py"
runAS = path + "/" + "runAS.py"
runAP = path + "/" + "runAP.py"
predictProd = path + "/" + "predictProductivity"
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
        self.parser = argparse.ArgumentParser(description = ' deFLAIR.py - a rpy2 convenience tool to run DESeq2.',
                                             epilog = 'Please feel free to forward any questions/concerns to /dev/null', 
                                             add_help = True, #default is True 
                                             prefix_chars = '-', 
                                             usage = '%(prog)s --manifest manifest.txt --workingdir dir_name --outdir out_dir --filter N')
        # Add args
        self.parser.add_argument("--outDir"    , action = 'store', required=True, 
                                    help='Write to specified output directory.')
        self.parser.add_argument("--filter"    , action = 'store', required=False, default = 10, type=int,
                                    help='Isoforms with less than specified read count for either Condition A or B are filtered (Default: 10 reads)')
        self.parser.add_argument("--matrix"    , action = 'store', required=True, 
                                    help='Count matrix from FLAIR quantify.')


        if inOpts is None :
            self.args = vars(self.parser.parse_args())
        else :
            self.args = vars(self.parser.parse_args(inOpts))


########################################################################
# Isoform
########################################################################


class Isoform(object) :
    '''
    Object to handle isoform related data.
    
    attributes:
        
    methods:
    
    ''' 

    def __init__(self, tid=None, gid=None):
        self.name = tid
        self.parent = gid

        self.exp  = None
        self.usage   = None

        self.deseq2AdjP  = float()
        self.drimseqAdjP = float()

        self.deseq2FC = float()
        self.deltaIPU = float()
        

    def computeUsage(self):

        self.usage = ["%.2f" % np.divide(iso,gene) for iso,gene in zip(self.exp,self.parent.exp)]
        #self.usage[np.isinf(self.usage)] = np.nan


########################################################################
# Gene
########################################################################


class Gene(object) :
    '''
    Object to handle gene related data.
    
    attributes:
        
    methods:
    
    ''' 

    def __init__(self, gid=None):
        self.transcripts = list()
        self.name = gid

        self.exp  = None

        self.deseq2AdjP  = float()
      
        self.deseq2FC = float()

########################################################################
# Funktions
########################################################################

def separateTables(f, thresh, samples, groups):

    genes, isoforms = dict(), dict()


    with open(f) as lines:
        next(lines)
        for line in lines:
            data = line.rstrip().split()
            name = data[0]
            counts = np.asarray(data[1:], dtype=float)
            iso, gene   = name, name.split("_")[-1]
            if "-" in gene:
                gene = gene.split("-")[0]
            m = iso.count("_")
            if m > 1:
                iso = iso.replace("_","",1)

            if gene not in genes:
                genes[gene] = Gene(gene)
                genes[gene].exp = np.zeros(len(counts))

            geneObj = genes[gene]            
            geneObj.exp += counts

            if iso not in isoforms:
                isoforms[iso] = Isoform(iso)
                isoforms[iso].exp = counts

            else:
                sys.stderr.write("Duplicate isoform entry %s. Remove isoform duplicate IDs. Exit.\n")
                sys.exit(1)

            isoformObj = isoforms[iso]
            isoformObj.parent = geneObj



    # get group indices for filterin tables
    groups = np.asarray(groups)
    g1Ind  = np.where(groups == groups[0])[0]
    g2Ind  = np.where(groups == groups[-1])[0]


    #make gene table first
    geneIDs = np.asarray(list(genes.keys()))
    vals    = np.asarray([genes[x].exp for x in geneIDs])
    filteredRows = (np.min(vals[:,g1Ind],axis=1) > thresh) |  (np.min(vals[:,g2Ind],axis=1) > thresh)
    filteredGeneVals = vals[filteredRows]
    filteredGeneIDs  = geneIDs[filteredRows]

    # now do isoforms
    isoformIDs = np.asarray(list(isoforms.keys()))
    vals = np.asarray([isoforms[x].exp for x in isoformIDs])
    filteredRows = (np.min(vals[:,g1Ind],axis=1) > thresh) |  (np.min(vals[:,g2Ind],axis=1) > thresh)
    filteredIsoVals = vals[filteredRows]
    filteredIsoIDs  = isoformIDs[filteredRows]

    geneDF  = pd.DataFrame(filteredGeneVals,columns=samples, index=filteredGeneIDs)
    isoDF = pd.DataFrame(filteredIsoVals,columns=samples, index=filteredIsoIDs)
    

    geneDF.to_csv("filtered_gene_counts_ds2.tsv", sep="\t")
    isoDF.to_csv("filtered_iso_counts_ds2.tsv", sep="\t")

    # also make table for drimm-seq
    isoformIDs = np.asarray([[y.parent.name,x] for x,y in isoforms.items()])
    vals = np.asarray([isoforms[x[-1]].exp for x in isoformIDs])
    #isoformIDs = isoformIDs.reshape(len(isoformIDs),1)
    allIso = np.hstack((isoformIDs,vals))
    isoDF  = pd.DataFrame(allIso, columns=['gene_id','feature_id']+samples)
    isoDF.to_csv("filtered_iso_counts_drim.tsv", sep="\t")

    return genes, isoforms

def makeDir(out):
    try:
        os.mkdir("./%s" % out)
    except:
        #exists
        pass


def main():

    '''
    maine
    '''

    # Command Line Stuff...
    myCommandLine = CommandLine()

    outDir     = myCommandLine.args['outDir']
    quantTable = myCommandLine.args['matrix']
    sFilter    = myCommandLine.args['filter']

    # Get sample data info
    with open(quantTable) as l:
        header = next(l).split()[1:]

    samples = header
    groups  = [x.split("_")[1] for x in header]
    batches = [x.split("_")[-1] for x in header]

    # Create output directory.
    makeDir(outDir)

    # Convert count tables to dataframe and update isoform objects.
    genes, isoforms = separateTables(quantTable, sFilter, samples, groups)

    # Compute percent isoform usgae
    [o.computeUsage() for i,o in isoforms.items()]

    header        = ['sample_id','condition','batch']
    formulaMatrix = [[x,y,z] for x,y,z in zip(samples,groups,batches)]
    formulaDF     = pd.DataFrame(formulaMatrix,columns=header)
    formulaDF     = formulaDF.set_index('sample_id')

    formulaMatrixFile = "formula_matrix.tsv"
    isoMatrixFile     = "filtered_iso_counts_ds2.tsv"
    geneMatrixFile    = "filtered_gene_counts_ds2.tsv"
    drimMatrixFile    = "filtered_iso_counts_drim.tsv"
    
    formulaDF.to_csv( formulaMatrixFile, sep='\t')

    subprocess.call([sys.executable, runDE, "--group1", groups[0], "--group2", groups[-1], 
                        "--batch", batches[0], "--matrix", geneMatrixFile, "--outDir", outDir,
                        "--prefix", "dge", "--formula", formulaMatrixFile], stderr=open("%s/dge_stderr.txt" % outDir,"w"))
    subprocess.call([sys.executable, runDE, "--group1", groups[0], "--group2", groups[-1], 
                        "--batch", batches[0], "--matrix", isoMatrixFile, "--outDir", outDir,
                        "--prefix", "die", "--formula", formulaMatrixFile], stderr=open("%s/dge_stderr.txt" % outDir,"w+"))
    subprocess.call([sys.executable, runDU, "--group1", groups[0], "--group2", groups[-1], 
                        "--batch", batches[0], "--matrix", drimMatrixFile, "--outDir", outDir,
                        "--prefix", "die", "--formula", formulaMatrixFile], stderr=open("%s/dge_stderr.txt" % outDir,"w+"))


if __name__ == "__main__":
    main()
