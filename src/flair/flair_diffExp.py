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
import argparse
import codecs
import errno
from collections import Counter
from statistics import median,stdev
from scipy.stats import ttest_ind, ranksums
from statsmodels.stats.multitest import multipletests
from statistics import median, mean
import pipettor

from flair import FlairError, check_diffexp_dependencies, set_unix_path

check_diffexp_dependencies()

os.environ['OPENBLAS_NUM_THREADS'] = '1'
import numpy as np
import pandas as pd



scriptPath = os.path.realpath(__file__)
path = "/".join(scriptPath.split("/")[:-1])
runDE = path + "/" + "runDE.py"
runDU = path + "/" + "runDU.py"


########################################################################
# Isoform
########################################################################


class Isoform(object):
    '''
    Object to handle isoform related data.

    '''

    def __init__(self, tid=None, gid=None):
        self.name = tid
        self.parent = gid

        self.exp  = None


########################################################################
# Gene
########################################################################


class Gene(object):
    '''
    Object to handle gene related data.

    '''

    def __init__(self, gid=None):
#        self.transcripts = list()
        self.name = gid

        self.exp  = None


########################################################################
# Functions
########################################################################


def get_gene_to_counts(filename):
    genetototcounts = {}
    for line in open(filename):
        line = line.rstrip().split('\t')
        if line[0] != 'ids':
            gene = line[0].split('_')[-1]
            counts = [int(x) for x in line[1:]]
            if gene not in genetototcounts:
                genetototcounts[gene] = [0 for x in range(len(counts))]
            genetototcounts[gene] = [genetototcounts[gene][x] + counts[x] for x in range(len(counts))]
    return genetototcounts


def do_mtc_ttest(filename, genetototcounts):
    allids, allpval, alldeltas = [], [], []
    for line in open(filename):
        line = line.rstrip().split('\t')
        if line[0] != 'ids':
            id = line[0]
            gene = id.split('_')[-1]
            wtcounts = [int(x) for x in line[1:4]]
            varcounts = [int(x) for x in line[4:]]
            deltaval = median(varcounts)-median(wtcounts)
            wttot, vartot = mean(genetototcounts[gene][:3]), mean(genetototcounts[gene][3:])
            deltausage = (mean(varcounts)/vartot if vartot > 0 else 0)-(mean(wtcounts)/wttot if wttot > 0 else 0)
            if abs(deltaval) > 3:
                pval = ttest_ind(wtcounts, varcounts).pvalue
                # pval = ranksums(wtcounts, varcounts).pvalue ###doesn't work, too stringent
                allids.append(id)
                allpval.append(pval)
                alldeltas.append(deltausage)
    corrpval = list(multipletests(allpval)[1])
    return allids, alldeltas, corrpval


def get_sig_from_norm_by_gene(outname, filename):
    """
    This function runs t-tests with multiple testing correction on a file of isoforms counts normalized by gene
    This method essentially does differential isoform usage testing, but accounts for differences in gene expression
    This is better for detecting novel transcripts than DRIM-seq
    """

    genetototcounts = get_gene_to_counts(filename)
    allids, alldeltas, corrpval = do_mtc_ttest(filename, genetototcounts)

    out = open(outname, 'w')
    for i in range(len(allids)):
        gene = allids[i]
        if corrpval[i] < 0.05:
            out.write('\t'.join([gene, str(round(alldeltas[i], 3)), str(corrpval[i])]) + '\n')



def separateTables(f, thresh, samples, groups, outDir):

    genes, isoforms = dict(), dict()
    duplicateID = 1

    with codecs.open(f, "r", encoding='utf-8', errors='ignore') as lines:
        cols = next(lines).split("\t")

        if len(cols) < 7:
            raise ValueError("** Error. Found %s columns in counts matrix, expected >6. Exiting." % len(cols))

        for num, line in enumerate(lines):

            data = line.rstrip().split("\t")
            if len(data) < 7:
                raise ValueError("** Error. Found %s columns in counts matrix, expected >6. Exiting." % len(cols))

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


def calc_gene_norm_sig(workdir, quantTable):
    """
    Make a file of counts normalized by gene.
    This is not a standard normalization method, only used for downstream stats.
    """
    genetosampletotot = {}
    lines = []
    out = open(workdir + '/counts.normbygene.tsv', 'w')
    for line in open(quantTable):
        line = line.rstrip().split('\t')
        if line[0] == 'ids':
            out.write('\t'.join(line) + '\n')
        else:
            oggenes = line[0].split('_')[-1]
            genes = oggenes.split('--')
            isoname = '_'.join(line[0].split('_')[:-1])
            l = len(line)
            for gene in genes:
                if '--' in oggenes:  ##is fusion
                    line[0] = isoname + '_' + oggenes + '_' + gene

                counts = [float(x) for x in line[1:]]
                lines.append([line[0]] + counts)
                if gene not in genetosampletotot: genetosampletotot[gene] = [0 for x in range(len(counts))]
                genetosampletotot[gene] = [genetosampletotot[gene][x] + counts[x] for x in range(len(counts))]
    for l in lines:
        gene = l[0].split('_')[-1]
        thisgenetot = genetosampletotot[gene]
        geneavg = sum(thisgenetot) / len(thisgenetot)
        thesecounts = l[1:]
        thesecounts = [(thesecounts[x] / thisgenetot[x]) * geneavg if thisgenetot[x] > 0 else 0 for x in
                       range(len(thisgenetot))]
        thesecounts = [str(round(x)) for x in thesecounts]
        out.write('\t'.join([l[0]] + thesecounts) + '\n')
    out.close()


def calculate_sig(args):

    '''
    '''

    outDir     = args.o
    quantTable = args.q
    sFilter    = args.e
    threads    = args.t
    force_dir  = args.of

    # Get sample data info
    with open(quantTable) as l:
        header = next(l).split()[1:]

    samples = ["%s_%s" % (h,num) for num,h in enumerate(header)]
    try:
        groups  = [x.split("_")[1] for x in header]
        batches = [x.split("_")[-1] for x in header]
        combos  = set([(groups.index(x),batches.index(y)) for x,y in zip(groups,batches)])
    except IndexError:
        raise Exception("** Error. diffExp requires column headers to contain sample, group, and batch, separated by '_'")
    except Exception as ex:
        raise Exception("** deFlair FAILED for %s" % (quantTable)) from ex

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


    calc_gene_norm_sig(workdir, quantTable)
    get_sig_from_norm_by_gene(outDir + '/isoforms_sig_exp_change_norm_by_gene.tsv', workdir + '/counts.normbygene.tsv')

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
            pipettor.run([sys.executable, runDE, "--group1", groups[0], "--group2", groups[-1],
                          "--batch", batches[0], "--matrix", geneMatrixFile, "--outDir", outDir,
                          "--prefix", "genes_deseq2", "--formula", formulaMatrixFile], stderr=out1)
        except pipettor.ProcessException as exc:
            raise FlairError(f'running DESeq2 on genes failed, please check {workdir}/dge_stderr.txt for details') from exc

        try:
            pipettor.run([sys.executable, runDE, "--group1", groups[0], "--group2", groups[-1],
                          "--batch", batches[0], "--matrix", isoMatrixFile, "--outDir", outDir,
                          "--prefix", "isoforms_deseq2", "--formula", formulaMatrixFile], stderr=out1)
        except pipettor.ProcessException as exc:
            raise FlairError(f'running DESeq2 on isoforms failed, please check {workdir}/dge_stderr.txt for details') from exc
        sys.stdout.flush()

        try:
            pipettor.run([sys.executable, runDU, "--threads", str(threads), "--group1", groups[0], "--group2", groups[-1],
                          "--batch", batches[0], "--matrix", drimMatrixFile, "--outDir", outDir,
                          "--prefix", "isoforms_drimseq", "--formula", formulaMatrixFile], stderr=out1)
        except pipettor.ProcessException as exc:
            raise FlairError(f'running DRIMSeq failed, please check {workdir}/dge_stderr.txt for details') from exc
        sys.stdout.flush()


def diffExp(counts_matrix=''):
    set_unix_path()
    parser = argparse.ArgumentParser(description='flair-diffExp parse options',
            usage='flair diffExp -q counts_matrix.tsv --out_dir out_dir [options]')
#       parser.add_argument('diffExp')
    required = parser.add_argument_group('required named arguments')
    if not counts_matrix:
        required.add_argument('-q', '--counts_matrix', action='store', dest='q',
                type=str, required=True, help='Tab-delimited isoform count matrix from flair quantify module.')
    required.add_argument('-o', '--out_dir', action='store', dest='o',
            type=str, required=True, help='Output directory for tables and plots.')
    parser.add_argument('-t', '--threads', action='store', dest='t',
            type=int, required=False, default=4, help='Number of threads for parallel DRIMSeq.')
    parser.add_argument('-e', '--exp_thresh', action='store', dest='e', type=int, required=False,
            default=10, help='''Read count expression threshold. Isoforms in which
            both conditions contain fewer than E reads are filtered out (Default E=10)''')
    parser.add_argument('-of', '--out_dir_force', action='store_true', dest='of',
            required=False, help='''Specify this argument to force overwriting of files in
            an existing output directory''')
    args = parser.parse_args()
    if counts_matrix:
        args.q = counts_matrix
        args.o = args.o+'.diffExp'

    if not os.path.exists(args.q):
        sys.stderr.write('Counts matrix file path does not exist\n')
        return 1

    calculate_sig(args)



if __name__ == "__main__":
    exit(diffExp())
