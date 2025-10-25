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


import os
import os.path as osp
import sys
import argparse
import errno
import csv
from collections import Counter
from scipy.stats import ttest_ind
from statistics import median, mean
import pipettor
import numpy as np

from flair import FlairError, FlairInputDataError, set_unix_path

os.environ['OPENBLAS_NUM_THREADS'] = '1'
import numpy as np


pkgdir = osp.dirname(osp.realpath(__file__))
diffExp_deseq2 = osp.join(pkgdir, "diffExp_deseq2.R")
diffExp_drimseq = osp.join(pkgdir, "diffExp_drimseq.R")


##
# Isoform and gene classes
##

class Isoform:
    '''
    Object to handle isoform related data.
    '''

    def __init__(self, tid, parent, counts):
        self.name = tid
        self.parent = parent
        self.exp = counts


class Gene(object):
    '''
    Object to handle gene related data.

    '''

    def __init__(self, gid, counts):
        self.name = gid
        self.exp = counts


########################################################################
# Functions
########################################################################


def multipletests(pvals, alpha=0.05, method='hs'):
    """adapted from statsmodels.stats.multitest
    does holm-sidak correction"""
    pvals = np.asarray(pvals)
    alphaf = alpha  # Notation ?

    sortind = np.argsort(pvals)
    pvals = np.take(pvals, sortind)

    ntests = len(pvals)
    alphacSidak = 1 - np.power((1. - alphaf), 1. / ntests)
    alphacBonf = alphaf / float(ntests)

    alphacSidak_all = 1 - np.power((1. - alphaf),
                                   1. / np.arange(ntests, 0, -1))
    notreject = pvals > alphacSidak_all
    del alphacSidak_all

    nr_index = np.nonzero(notreject)[0]
    if nr_index.size == 0:
        # nonreject is empty, all rejected
        notrejectmin = len(pvals)
    else:
        notrejectmin = np.min(nr_index)
    notreject[notrejectmin:] = True
    reject = ~notreject
    del notreject

    pvals_corrected_raw = 1 - np.power((1. - pvals),
                                       np.arange(ntests, 0, -1))
    pvals_corrected = np.maximum.accumulate(pvals_corrected_raw)
    del pvals_corrected_raw

    if not pvals_corrected is None:  # not necessary anymore
        pvals_corrected[pvals_corrected > 1] = 1
    pvals_corrected_ = np.empty_like(pvals_corrected)
    pvals_corrected_[sortind] = pvals_corrected
    del pvals_corrected
    reject_ = np.empty_like(reject)
    reject_[sortind] = reject
    return reject_, pvals_corrected_, alphacSidak, alphacBonf


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
            deltaval = median(varcounts) - median(wtcounts)
            wttot, vartot = mean(genetototcounts[gene][:3]), mean(genetototcounts[gene][3:])
            deltausage = (mean(varcounts) / vartot if vartot > 0 else 0) - (mean(wtcounts) / wttot if wttot > 0 else 0)
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

def quant_row_check(linenum, row):
    if len(row) < 7:
        raise FlairInputDataError(f"line {linenum}: found {len(row)} columns in counts matrix, expected >6")

def quant_table_reader(quant_table_tsv):
    """Generator for rows of (name, counts) from counts file"""
    try:
        with open(quant_table_tsv, "r", encoding='utf-8', errors='ignore') as fh:
            csvreader = csv.reader(fh, delimiter='\t')
            cols = next(csvreader)
            quant_row_check(1, cols)
            for linenum, row in enumerate(csvreader, start=2):
                quant_row_check(linenum, row)
                yield row[0], np.asarray(row[1:], dtype=float)
    except Exception as exc:
        raise FlairInputDataError(f"error parsing counts table: {quant_table_tsv}") from exc

def write_name_values_tsv(samples, names, values, out_tsv):
    "write gene or isoform matrix tsv"
    with open(out_tsv, 'w') as fh:
        writer = csv.writer(fh, delimiter='\t', dialect='unix', quoting=csv.QUOTE_NONE)
        writer.writerow([''] + samples)
        for name, value in zip(names, values):
            writer.writerow([name] + list(value))

def write_tsv(columns, rows, out_tsv):
    """write a TSV.  """
    with open(out_tsv, 'w') as fh:
        writer = csv.writer(fh, delimiter='\t', dialect='unix', quoting=csv.QUOTE_NONE)
        writer.writerow(columns)
        for row in rows:
            writer.writerow(row)

def separate_tables(quant_table_tsv, thresh, samples, groups, outDir):
    genes, isoforms = dict(), dict()
    duplicateID = 1

    for name, counts in quant_table_reader(quant_table_tsv):
        # FIXME: gene name parsing needs to be moved to a common module
        iso, gene = name, name.split("_")[-1]
        if "-" in gene:
            gene = gene.split("-")[0]
        m = iso.count("_")
        if m > 1:
            iso = iso.replace("_", "", 1)

        if gene not in genes:
            genes[gene] = Gene(gene, np.zeros(len(counts)))

        geneObj = genes[gene]
        geneObj.exp += counts

        if iso not in isoforms:
            isoforms[iso] = Isoform(iso, geneObj, counts)
        else:
            duplicateID += 1
            iso = iso + "-" + str(duplicateID)
            isoforms[iso] = Isoform(iso, geneObj, counts)

    # get group indices for filterin tables
    groups = np.asarray(groups)
    g1Ind = np.where(groups == groups[0])[0]
    g2Ind = np.where(groups == groups[-1])[0]

    # make gene table first
    geneIDs = np.asarray(list(genes.keys()))
    vals = np.asarray([genes[x].exp for x in geneIDs])
    if len(geneIDs) == 0:
        raise FlairInputDataError(f"no genes parsed from {quant_table_tsv}")

    # genes must be expressed in all samples of at least one group
    filteredRows = (np.min(vals[:, g1Ind], axis=1) > thresh) | (np.min(vals[:, g2Ind], axis=1) > thresh)
    filteredGeneVals = vals[filteredRows]
    filteredGeneIDs = geneIDs[filteredRows]
    write_name_values_tsv(samples, filteredGeneIDs, filteredGeneVals,
                          outDir + "/filtered_gene_counts_ds2.tsv")

    # now do isoforms
    isoformIDs = np.asarray(list(isoforms.keys()))
    vals = np.asarray([isoforms[x].exp for x in isoformIDs])
    filteredRows = (np.min(vals[:, g1Ind], axis=1) > thresh) | (np.min(vals[:, g2Ind], axis=1) > thresh)
    filteredIsoVals = vals[filteredRows]
    filteredIsoIDs = isoformIDs[filteredRows]

    write_name_values_tsv(samples, filteredIsoIDs, filteredIsoVals,
                          outDir + "/filtered_iso_counts_ds2.tsv")

    # also make table for drimm-seq.  It must have a unique row undex
    # added to prevent 'DataFrame contains duplicated elements in the index'
    isoformIDs = np.asarray([[y.parent.name, x] for x, y in isoforms.items()])
    vals = np.asarray([isoforms[x[-1]].exp for x in isoformIDs])
    indices = np.arange(isoformIDs.shape[0]).reshape(-1, 1)
    allIso = np.hstack((indices, isoformIDs, vals))
    write_tsv(['irow', 'gene_id', 'feature_id'] + samples, allIso,
              outDir + "/filtered_iso_counts_drim.tsv")
    return genes, isoforms


def calc_gene_norm_sig(workdir, quant_table_tsv):
    """
    Make a file of counts normalized by gene.
    This is not a standard normalization method, only used for downstream stats.
    """
    genetosampletotot = {}
    lines = []
    # FIXME: just load into memory rather than reading three timnes!!
    out = open(workdir + '/counts.normbygene.tsv', 'w')
    for line in open(quant_table_tsv):
        line = line.rstrip().split('\t')
        if line[0] == 'ids':
            out.write('\t'.join(line) + '\n')
        else:
            oggenes = line[0].split('_')[-1]
            genes = oggenes.split('--')
            isoname = '_'.join(line[0].split('_')[:-1])
            l = len(line)
            for gene in genes:
                if '--' in oggenes:  # is fusion ?
                    line[0] = isoname + '_' + oggenes + '_' + gene

                counts = [float(x) for x in line[1:]]
                lines.append([line[0]] + counts)
                if gene not in genetosampletotot:
                    genetosampletotot[gene] = [0 for x in range(len(counts))]
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

def run_deseq2(prefix, workdir, groups, batches, matrixFile, outDir, formulaMatrixFile):
    stderr = f"{workdir}/{prefix}.txt"
    try:
        with open(stderr, "w") as stderr_fh:
            pipettor.run(["Rscript", diffExp_deseq2, "--group1", groups[0], "--group2", groups[-1],
                          "--batch", batches[0], "--matrix", matrixFile, "--outDir", outDir,
                          "--prefix", prefix, "--formula", formulaMatrixFile], stderr=stderr_fh)
    except pipettor.ProcessException as exc:
        raise FlairError(f'running {prefix} failed, please check {stderr} for details') from exc

def run_dirmseq(prefix, workdir, threads, groups, batches, matrixFile, outDir, formulaMatrixFile):
    stderr = f"{workdir}/{prefix}.txt"
    try:
        with open(stderr, "w") as stderr_fh:
            pipettor.run(["Rscript", diffExp_drimseq, "--threads", threads, "--group1", groups[0], "--group2", groups[-1],
                          "--batch", batches[0], "--matrix", matrixFile, "--outDir", outDir,
                          "--prefix", prefix, "--formula", formulaMatrixFile], stderr=stderr_fh)
    except pipettor.ProcessException as exc:
        raise FlairError(f'running {prefix} failed, please check {stderr} for details') from exc


def calculate_sig(args):
    outDir = args.out_dir
    quant_table_tsv = args.counts_matrix
    sFilter = args.exp_thresh
    threads = args.threads
    force_dir = args.out_dir_force

    # FIXME convert to just loading table upfront
    # Get sample data info
    with open(quant_table_tsv) as l:
        header = next(l).split()[1:]

    samples = ["%s_%s" % (h, num) for num, h in enumerate(header)]
    try:
        groups = [x.split("_")[1] for x in header]
        batches = [x.split("_")[-1] for x in header]
        combos = set([(groups.index(x), batches.index(y)) for x, y in zip(groups, batches)])
    except IndexError:
        raise FlairInputDataError("diffExp requires column headers to contain sample, group, and batch, separated by '_'")

    groupCounts = Counter(groups)
    if len(list(groupCounts.keys())) != 2:
        raise FlairInputDataError("** Error. diffExp requires exactly 2 condition groups. Maybe group name formatting is incorrect")
    elif min(list(groupCounts.values())) < 3:
        raise FlairInputDataError("** Error. diffExp requires >2 samples per condition group. Use diff_iso_usage.py for analyses with <3 replicates.")
    elif set(groups).intersection(set(batches)):
        raise FlairInputDataError("** Error. Sample group/condition names and batch descriptor must be distinct. Try renaming batch descriptor in count matrix.")
    elif sum([1 if x.isdigit() else 0 for x in groups]) > 0 or sum([1 if x.isdigit() else 0 for x in batches]) > 0:
        raise FlairInputDataError("** Error. Sample group/condition or batch names are required to be strings not integers. Please change formatting.")

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
        raise FlairInputDataError(f"** Error. Name {outDir} already exists. Choose another name for out_dir")

    calc_gene_norm_sig(workdir, quant_table_tsv)
    get_sig_from_norm_by_gene(outDir + '/isoforms_sig_exp_change_norm_by_gene.tsv', workdir + '/counts.normbygene.tsv')

    # Convert count tables to dataframe and update isoform objects.
    genes, isoforms = separate_tables(quant_table_tsv, sFilter, samples, groups, workdir)

    # checks linear combination
    if len(combos) == 2:
        header = ['sample_id', 'condition']
        formulaMatrix = [[x, y] for x, y in zip(samples, groups)]
    elif len(set(batches)) > 1:
        header = ['sample_id', 'condition', 'batch']
        formulaMatrix = [[x, y, z] for x, y, z in zip(samples, groups, batches)]
    else:
        header = ['sample_id', ' condition']
        formulaMatrix = [[x, y] for x, y in zip(samples, groups)]

    formulaMatrixFile = workdir + "/formula_matrix.tsv"
    write_tsv(header, formulaMatrix, formulaMatrixFile)

    isoMatrixFile = workdir + "/filtered_iso_counts_ds2.tsv"
    geneMatrixFile = workdir + "/filtered_gene_counts_ds2.tsv"
    drimMatrixFile = workdir + "/filtered_iso_counts_drim.tsv"

    # DESeq2 genes & isoforms
    run_deseq2("genes_deseq2", workdir, groups, batches, geneMatrixFile, outDir, formulaMatrixFile)
    run_deseq2("isoforms_deseq2", workdir, groups, batches, isoMatrixFile, outDir, formulaMatrixFile)

    # DIRMSeq
    run_dirmseq("isoforms_drimseq", workdir, threads, groups, batches, drimMatrixFile, outDir, formulaMatrixFile)

def diffExp(counts_matrix=''):
    set_unix_path()
    parser = argparse.ArgumentParser()
#       parser.add_argument('diffExp')
    required = parser.add_argument_group('required named arguments')
    if not counts_matrix:
        required.add_argument('-q', '--counts_matrix', action='store',
                              type=str, required=True, help='Tab-delimited isoform count matrix from flair quantify module.')
    required.add_argument('-o', '--out_dir', action='store',
                          type=str, required=True, help='Output directory for tables and plots.')
    parser.add_argument('-t', '--threads', action='store',
                        type=int, required=False, default=4, help='Number of threads for parallel DRIMSeq.')
    parser.add_argument('-e', '--exp_thresh', action='store', type=int, required=False,
                        default=10, help='''Read count expression threshold. Isoforms in which
                        both conditions contain fewer than E reads are filtered out (Default E=10)''')
    parser.add_argument('-of', '--out_dir_force', action='store_true',
                        required=False, help='''Specify this argument to force overwriting of files in
                        an existing output directory''')
    args = parser.parse_args()
    if counts_matrix:
        args.counts_matrix = counts_matrix
        args.out_dir = args.out_dir + '.diffExp'

    if not os.path.exists(args.counts_matrix):
        raise FlairInputDataError('Counts matrix file path does not exist')

    calculate_sig(args)

if __name__ == "__main__":
    exit(diffExp())
