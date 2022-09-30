#!/usr/bin/env python3

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


import os
import sys
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import pandas as pd
import numpy as np
import subprocess
import codecs
import errno
from collections import Counter

scriptPath = os.path.realpath(__file__)
path = "/".join(scriptPath.split("/")[:-1])
runDE = path + "/" + "runDE.py"
runDU = path + "/" + "runDU.py"

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
        import argparse
        self.parser = argparse.ArgumentParser(description=' deFLAIR.py - a rpy2 convenience tool to run DESeq2.',
                                             add_help=True, #default is True
                                             prefix_chars='-',
                                             usage='%(prog)s --manifest manifest.txt --workingdir dir_name --outdir out_dir --filter N')
        # Add args
        self.parser.add_argument("--outDir", action='store', required=True,
                                    help='Write to specified output directory.')
        self.parser.add_argument("--filter", action='store', required=False, default=10, type=int,
                                    help='Isoforms with less than specified read count for either Condition A or B are filtered (Default: 10 reads)')
        self.parser.add_argument("--matrix", action='store', required=True,
                                    help='Count matrix from FLAIR quantify.')
        self.parser.add_argument("--threads", type=int, action='store', required=False, default=4,
                                    help='Number of threads for running DRIM-Seq.')
        self.parser.add_argument('-of', '--out_dir_force', action='store_true', dest='of',
            required=False, help='''Specify this argument to force overwriting of
            an existing output directory for tables and plots.''')

        if inOpts is None:
            self.args = vars(self.parser.parse_args())
        else:
            self.args = vars(self.parser.parse_args(inOpts))


########################################################################
# Isoform
########################################################################


class Isoform(object):
    '''
    Object to handle isoform related data.

    attributes:

    methods:

    '''

    def __init__(self, tid=None, gid=None):
        self.name = tid
        self.parent = gid

        self.exp  = None
#        self.usage   = None

#        self.deseq2AdjP  = float()
#        self.drimseqAdjP = float()

#        self.deseq2FC = float()
#        self.deltaIPU = float()


########################################################################
# Gene
########################################################################


class Gene(object):
    '''
    Object to handle gene related data.

    attributes:

    methods:

    '''

    def __init__(self, gid=None):
#        self.transcripts = list()
        self.name = gid

        self.exp  = None

#        self.deseq2AdjP  = float()

#        self.deseq2FC = float()

########################################################################
# Functions
########################################################################


def separateTables(f, thresh, samples, groups, outDir):

    genes, isoforms = dict(), dict()
    duplicateID = 1

    with codecs.open(f, "r", encoding='utf-8', errors='ignore') as lines:
        cols = next(lines).split("\t")

        if len(cols) < 7:
            print("** Error. Found %s columns in counts matrix, expected >6. Exiting." % len(cols),file=sys.stderr)
            sys.exit(1)

        for num, line in enumerate(lines):

            data = line.rstrip().split("\t")
            if len(data) < 7:
                print("** Error. Found %s columns in matrix on line %s, expected >6. Exiting." % (len(data), num),file=sys.stderr)
                sys.exit(1)

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
                duplicateID += 1
                iso = iso + "-" + str(duplicateID)
                isoforms[iso] = Isoform(iso)
                isoforms[iso].exp = counts

            isoformObj = isoforms[iso]
            isoformObj.parent = geneObj

    # get group indices for filterin tables
    groups = np.asarray(groups)
    g1Ind  = np.where(groups == groups[0])[0]
    g2Ind  = np.where(groups == groups[-1])[0]

    #make gene table first
    geneIDs = np.asarray(list(genes.keys()))
    vals    = np.asarray([genes[x].exp for x in geneIDs])
    # genes must be expressed in all samples of at least one group
    filteredRows = (np.min(vals[:,g1Ind],axis=1) > thresh) | (np.min(vals[:,g2Ind],axis=1) > thresh)
    filteredGeneVals = vals[filteredRows]
    filteredGeneIDs  = geneIDs[filteredRows]

    # now do isoforms
    isoformIDs = np.asarray(list(isoforms.keys()))
    vals = np.asarray([isoforms[x].exp for x in isoformIDs])
    filteredRows = (np.min(vals[:,g1Ind],axis=1) > thresh) | (np.min(vals[:,g2Ind],axis=1) > thresh)
    filteredIsoVals = vals[filteredRows]
    filteredIsoIDs  = isoformIDs[filteredRows]

    geneDF  = pd.DataFrame(filteredGeneVals,columns=samples, index=filteredGeneIDs)
    isoDF = pd.DataFrame(filteredIsoVals,columns=samples, index=filteredIsoIDs)

    geneDF.to_csv(outDir + "/filtered_gene_counts_ds2.tsv", sep="\t")
    isoDF.to_csv(outDir + "/filtered_iso_counts_ds2.tsv", sep="\t")

    # also make table for drimm-seq
    isoformIDs = np.asarray([[y.parent.name,x] for x,y in isoforms.items()])
    vals = np.asarray([isoforms[x[-1]].exp for x in isoformIDs])
    #isoformIDs = isoformIDs.reshape(len(isoformIDs),1)
    allIso = np.hstack((isoformIDs,vals))
    isoDF  = pd.DataFrame(allIso, columns=['gene_id','feature_id']+samples)
    isoDF.to_csv(outDir + "/filtered_iso_counts_drim.tsv", sep="\t")
    return genes, isoforms


def main():

    '''
    '''

    myCommandLine = CommandLine()

    outDir     = myCommandLine.args['outDir']
    quantTable = myCommandLine.args['matrix']
    sFilter    = myCommandLine.args['filter']
    threads    = myCommandLine.args['threads']
    force_dir  = myCommandLine.args['of']

    # Get sample data info
    with open(quantTable) as l:
        header = next(l).split()[1:]

    samples = ["%s_%s" % (h,num) for num,h in enumerate(header)]
    try:
        groups  = [x.split("_")[1] for x in header]
        batches = [x.split("_")[-1] for x in header]
        combos  = set([(groups.index(x),batches.index(y)) for x,y in zip(groups,batches)])
    except IndexError:
        print("** Error. diffExp requires column headers to contain sample, group, and batch, separated by '_'", file=sys.stderr)
        sys.exit(1)

    groupCounts = Counter(groups)
    if len(list(groupCounts.keys())) != 2:
        print("** Error. diffExp requires exactly 2 condition groups. Maybe group name formatting is incorrect. Exiting.", file=sys.stderr)
        sys.exit(1)
    elif min(list(groupCounts.values())) < 3:
        print("** Error. diffExp requires >2 samples per condition group. Use diff_iso_usage.py for analyses with <3 replicates.", file=sys.stderr)
        sys.exit(1)
    elif set(groups).intersection(set(batches)):
        print("** Error. Sample group/condition names and batch descriptor must be distinct. Try renaming batch descriptor in count matrix.", file=sys.stderr)
        sys.exit(1)
    elif sum([1 if x.isdigit() else 0 for x in groups]) > 0 or sum([1 if x.isdigit() else 0 for x in batches]) > 0:
        print("** Error. Sample group/condition or batch names are required to be strings not integers. Please change formatting.", file=sys.stderr)
        sys.exit(1)

    # Create output directory including a working directory for intermediate files.
    workdir = os.path.join(outDir, 'workdir')
    
    if force_dir:
        if not os.path.exists(workdir):
            os.makedirs(workdir)
        pass
    elif not os.path.exists(outDir):
        try:
            os.makedirs(workdir, 0o700)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise
    else:
        print("** Error. Name '%s' already exists. Choose another name for out_dir" % outDir, file=sys.stderr)
        sys.exit(1)
    # Convert count tables to dataframe and update isoform objects.
    genes, isoforms = separateTables(quantTable, sFilter, samples, groups, workdir)

    # checks linear combination
    if len(combos) == 2:
        header        = ['sample_id','condition']
        formulaMatrix = [[x,y] for x,y in zip(samples,groups)]
        formulaDF     = pd.DataFrame(formulaMatrix,columns=header)
        formulaDF     = formulaDF.set_index('sample_id')

    elif len(set(batches)) > 1:
        header        = ['sample_id','condition','batch']
        formulaMatrix = [[x,y,z] for x,y,z in zip(samples,groups,batches)]
        formulaDF     = pd.DataFrame(formulaMatrix,columns=header)
        formulaDF     = formulaDF.set_index('sample_id')

    else:
        header        = ['sample_id','condition']
        formulaMatrix = [[x,y] for x,y in zip(samples,groups)]
        formulaDF     = pd.DataFrame(formulaMatrix,columns=header)
        formulaDF     = formulaDF.set_index('sample_id')

    formulaMatrixFile = workdir + "/formula_matrix.tsv"
    isoMatrixFile     = workdir + "/filtered_iso_counts_ds2.tsv"
    geneMatrixFile    = workdir + "/filtered_gene_counts_ds2.tsv"
    drimMatrixFile    = workdir + "/filtered_iso_counts_drim.tsv"

    formulaDF.to_csv(formulaMatrixFile, sep='\t')

    with open("%s/dge_stderr.txt" % workdir,"w") as out1:

        try:
            subprocess.check_call([sys.executable, runDE, "--group1", groups[0], "--group2", groups[-1],
                            "--batch", batches[0], "--matrix", geneMatrixFile, "--outDir", outDir,
                            "--prefix", "genes_deseq2", "--formula", formulaMatrixFile], stderr=out1)
        except subprocess.CalledProcessError:
            print(f'WARNING, running DESeq2 on genes failed, please check {workdir}/dge_stderr.txt for details')
        sys.stdout.flush()
    
        try:
            subprocess.check_call([sys.executable, runDE, "--group1", groups[0], "--group2", groups[-1],
                            "--batch", batches[0], "--matrix", isoMatrixFile, "--outDir", outDir,
                            "--prefix", "isoforms_deseq2", "--formula", formulaMatrixFile], stderr=out1)
        except subprocess.CalledProcessError:
            print(f'WARNING, running DESeq2 on isoforms failed, please check {workdir}/dge_stderr.txt for details')
        sys.stdout.flush()
    
        try:
            subprocess.check_call([sys.executable, runDU, "--threads", str(threads), "--group1", groups[0], "--group2", groups[-1],
                             "--batch", batches[0], "--matrix", drimMatrixFile, "--outDir", outDir,
                             "--prefix", "isoforms_drimseq", "--formula", formulaMatrixFile], stderr=out1)
        except subprocess.CalledProcessError:
            print(f'WARNING, running DRIMSeq failed, please check {workdir}/dge_stderr.txt for details')
        sys.stdout.flush()


if __name__ == "__main__":
    myCommandLine = CommandLine()

    outDir     = myCommandLine.args['outDir']
    quantTable = myCommandLine.args['matrix']
    sFilter    = myCommandLine.args['filter']
    threads    = myCommandLine.args['threads']
    force_dir  = myCommandLine.args['of']
    main()
