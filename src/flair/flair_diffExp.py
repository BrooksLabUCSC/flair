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
import errno
import csv
from collections import Counter
from statistics import median, mean
import pipettor

from flair import FlairError, FlairInputDataError
from flair.counts_matrix import (read_sample_columns, parse_sample_fields,
                                 condition_column_indexes, select_condition_pair)

os.environ['OPENBLAS_NUM_THREADS'] = '1'
import numpy as np  # noqa: E402


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


def multipletests(pvals, alpha=0.05):
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

    if pvals_corrected is not None:  # not necessary anymore
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


def do_mtc_ttest(filename, genetototcounts, ref_cols, test_cols):
    # imported here rather than at module scope; scipy takes 0.7s to import and
    # this is the only use of it, which would be paid by every flair command
    from scipy.stats import ttest_ind
    allids, allpval, alldeltas = [], [], []
    for line in open(filename):
        line = line.rstrip().split('\t')
        if line[0] != 'ids':
            id = line[0]
            gene = id.split('_')[-1]

            counts = [int(x) for x in line[1:]]
            wtcounts = [counts[i] for i in ref_cols]
            varcounts = [counts[i] for i in test_cols]

            # Compute median difference between variant and WT
            deltaval = median(varcounts) - median(wtcounts)
            genetot = genetototcounts[gene]
            wttot = mean([genetot[i] for i in ref_cols])
            vartot = mean([genetot[i] for i in test_cols])

            # Compute normalized usage difference, ignore if totals are zero
            deltausage = (mean(varcounts) / vartot if vartot > 0 else 0) - (mean(wtcounts) / wttot if wttot > 0 else 0)

            # Only test if median difference is large enough (|Δ| > 3)
            if abs(deltaval) > 3:
                # Two-sample t-test between WT and VAR counts
                # ranksums() wast too strict for small replicates
                pval = ttest_ind(wtcounts, varcounts).pvalue

                allids.append(id)
                allpval.append(pval)
                alldeltas.append(deltausage)

    if len(allpval) == 0:
        raise FlairInputDataError(f"no p-values with sufficient delta values from: {filename}")

    # Apply multiple-testing correction.  This is Holm-Sidak, not Benjamini-Hochberg
    # as this comment used to say: the adjusted values control the family-wise error
    # rate, not the false discovery rate

    corrpval = list(multipletests(allpval)[1])
    return allids, alldeltas, corrpval


def get_sig_from_norm_by_gene(outname, filename, ref_cols, test_cols):
    """
    This function runs t-tests with multiple testing correction on a file of isoforms counts normalized by gene
    This method essentially does differential isoform usage testing, but accounts for differences in gene expression
    This is better for detecting novel transcripts than DRIM-seq
    """

    genetototcounts = get_gene_to_counts(filename)
    allids, alldeltas, corrpval = do_mtc_ttest(filename, genetototcounts, ref_cols, test_cols)

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

def separate_tables(quant_table_tsv, thresh, samples, a_cols, b_cols, outDir):
    genes, isoforms = dict(), dict()
    duplicateID = 1

    for name, counts in quant_table_reader(quant_table_tsv):
        # FIXME: gene name parsing needs to be moved to a common module
        iso, gene = name, name.split("_")[-1]
        # if "-" in gene:
        #     gene = gene.split("-")[0]
        # m = iso.count("_")
        # if m > 1:
        #     iso = iso.replace("_", "", 1)

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

    # the two conditions' column indexes are chosen by name in calculate_sig; taking
    # them from the first and last column made the filter depend on column order, and
    # with an order like A,B,B,A it tested one condition twice
    g1Ind = np.asarray(a_cols)
    g2Ind = np.asarray(b_cols)

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

def run_deseq2(prefix, workdir, condition_a, condition_b, matrixFile, outDir, formulaMatrixFile):
    # no --batch: neither R script ever read it, and one arbitrary batch label would
    # not have said anything anyway.  Both take the batch column from the formula matrix
    stderr = f"{workdir}/{prefix}.txt"
    try:
        with open(stderr, "w") as stderr_fh:
            pipettor.run(["Rscript", diffExp_deseq2, "--condition_a", condition_a, "--condition_b", condition_b,
                          "--matrix", matrixFile, "--out_dir", outDir,
                          "--prefix", prefix, "--formula", formulaMatrixFile], stderr=stderr_fh)
    except pipettor.ProcessException as exc:
        raise FlairError(f'running {prefix} failed, please check {stderr} for details') from exc

def run_dirmseq(prefix, workdir, threads, condition_a, condition_b, matrixFile, outDir, formulaMatrixFile):
    stderr = f"{workdir}/{prefix}.txt"
    try:
        with open(stderr, "w") as stderr_fh:
            pipettor.run(["Rscript", diffExp_drimseq, "--threads", threads, "--condition_a", condition_a, "--condition_b", condition_b,
                          "--matrix", matrixFile, "--out_dir", outDir,
                          "--prefix", prefix, "--formula", formulaMatrixFile], stderr=stderr_fh)
    except pipettor.ProcessException as exc:
        raise FlairError(f'running {prefix} failed, please check {stderr} for details') from exc


def calculate_sig(*, counts_matrix, output, condition_a, condition_b, min_expression,  # noqa: C901 - FIXME: reduce complexity
                  threads, overwrite_output):
    outDir = output
    quant_table_tsv = counts_matrix
    sFilter = min_expression
    force_dir = overwrite_output

    # FIXME convert to just loading table upfront
    # Get sample data info
    header = read_sample_columns(quant_table_tsv)
    samples = ["%s_%s" % (h, num) for num, h in enumerate(header)]
    groups, batches = parse_sample_fields(header, quant_table_tsv)
    combos = set([(groups.index(x), batches.index(y)) for x, y in zip(groups, batches)])

    condition_a, condition_b = select_condition_pair(groups, condition_a, condition_b,
                                                     quant_table_tsv)
    a_cols = condition_column_indexes(groups, condition_a)
    b_cols = condition_column_indexes(groups, condition_b)

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
    # counts.normbygene.tsv keeps the counts matrix header, so the same column
    # indexes apply to it
    get_sig_from_norm_by_gene(outDir + '/isoforms_sig_exp_change_norm_by_gene.tsv',
                              workdir + '/counts.normbygene.tsv', a_cols, b_cols)

    # Convert count tables to dataframe and update isoform objects.
    genes, isoforms = separate_tables(quant_table_tsv, sFilter, samples, a_cols, b_cols, workdir)

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
    run_deseq2("genes_deseq2", workdir, condition_a, condition_b, geneMatrixFile, outDir, formulaMatrixFile)
    run_deseq2("isoforms_deseq2", workdir, condition_a, condition_b, isoMatrixFile, outDir, formulaMatrixFile)

    # DIRMSeq
    run_dirmseq("isoforms_drimseq", workdir, threads, condition_a, condition_b, drimMatrixFile, outDir, formulaMatrixFile)

def add_subparser(subparsers):
    desc = "Differential expression and differential usage analysis"
    parser = subparsers.add_parser('diffexp', help="Differential expression and usage analysis",
                                   description=desc)
    required = parser.add_argument_group('required named arguments')
    required.add_argument('--counts_matrix', type=str, required=True,
                          help='tab-delimited isoform count matrix from flair quantify')
    required.add_argument('-o', '--output', type=str, required=True,
                          help='output directory for tables and plots')
    parser.add_argument('-t', '--threads', type=int, default=4,
                        help='number of threads for parallel DRIMSeq (default: %(default)s)')
    parser.add_argument('--min_expression', type=int, default=10,
                        help='read count expression threshold; isoforms in which both conditions '
                             'contain fewer than this many reads are filtered out (default: %(default)s)')
    parser.add_argument('--condition_a', default='',
                        help='the reference condition, as named in the counts matrix columns; '
                             'fold changes are reported for --condition_b relative to this. '
                             'With neither this nor --condition_b given, the two conditions are '
                             'taken in sorted order')
    parser.add_argument('--condition_b', default='',
                        help='the condition compared against --condition_a')
    parser.add_argument('--overwrite_output', action='store_true',
                        help='overwrite files in an existing output directory')
    parser.set_defaults(entry=diffexp_cmd)

def diffexp_cmd(args):
    if not os.path.exists(args.counts_matrix):
        raise FlairInputDataError(f'counts matrix file does not exist: {args.counts_matrix}')
    calculate_sig(counts_matrix=args.counts_matrix, output=args.output,
                  condition_a=args.condition_a, condition_b=args.condition_b,
                  min_expression=args.min_expression, threads=args.threads,
                  overwrite_output=args.overwrite_output)
