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
from collections import Counter
from statistics import median, mean
import pipettor

from flair import FlairError, FlairInputDataError
from flair.conditions import condition_column_indexes, select_condition_pair
from flair.counts_matrix_tsv import read_counts_rows, CountsRow, write_counts_matrix
from flair.deseq2_counts_tsv import write_deseq2_counts
from flair.drimseq_counts_tsv import drimseq_counts_writer
from flair.formula_matrix_tsv import FormulaRow, write_formula_matrix
from flair.iso_usage_change_tsv import iso_usage_change_writer
from flair.output_compress import compressed_path, output_compression
from flair.sample_info_tsv import read_sample_info

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


def get_gene_to_counts(rows):
    genetototcounts = {}
    for row in rows:
        counts = [int(x) for x in row.counts]
        if row.gene_id not in genetototcounts:
            genetototcounts[row.gene_id] = [0 for x in range(len(counts))]
        genetototcounts[row.gene_id] = [genetototcounts[row.gene_id][x] + counts[x] for x in range(len(counts))]
    return genetototcounts


def row_ttest(ttest_ind, row, genetot, ref_cols, test_cols):
    """The t-test of one isoform, or None when the two conditions differ by too
    little to be worth testing."""
    counts = [int(x) for x in row.counts]
    wtcounts = [counts[i] for i in ref_cols]
    varcounts = [counts[i] for i in test_cols]

    # Compute median difference between variant and WT
    deltaval = median(varcounts) - median(wtcounts)
    wttot = mean([genetot[i] for i in ref_cols])
    vartot = mean([genetot[i] for i in test_cols])

    # Compute normalized usage difference, ignore if totals are zero
    deltausage = (mean(varcounts) / vartot if vartot > 0 else 0) - (mean(wtcounts) / wttot if wttot > 0 else 0)

    # Only test if median difference is large enough (|Δ| > 3)
    if abs(deltaval) <= 3:
        return None
    # Two-sample t-test between WT and VAR counts
    # ranksums() wast too strict for small replicates
    return deltausage, ttest_ind(wtcounts, varcounts).pvalue

def do_mtc_ttest(rows, genetototcounts, ref_cols, test_cols):
    # imported here rather than at module scope; scipy takes 0.7s to import and
    # this is the only use of it, which would be paid by every flair command
    from scipy.stats import ttest_ind
    allids, allpval, alldeltas = [], [], []
    for row in rows:
        tested = row_ttest(ttest_ind, row, genetototcounts[row.gene_id], ref_cols, test_cols)
        if tested is not None:
            deltausage, pval = tested
            allids.append((row.gene_id, row.isoform_id))
            allpval.append(pval)
            alldeltas.append(deltausage)

    if len(allpval) == 0:
        raise FlairInputDataError("no p-values with sufficient delta values from the counts matrix")

    # Apply multiple-testing correction.  This is Holm-Sidak, not Benjamini-Hochberg
    # as this comment used to say: the adjusted values control the family-wise error
    # rate, not the false discovery rate

    corrpval = list(multipletests(allpval)[1])
    return allids, alldeltas, corrpval


def get_sig_from_norm_by_gene(outname, norm_rows, ref_cols, test_cols):
    """
    This function runs t-tests with multiple testing correction on isoform counts normalized by gene
    This method essentially does differential isoform usage testing, but accounts for differences in gene expression
    This is better for detecting novel transcripts than DRIM-seq
    """

    genetototcounts = get_gene_to_counts(norm_rows)
    allids, alldeltas, corrpval = do_mtc_ttest(norm_rows, genetototcounts, ref_cols, test_cols)

    with iso_usage_change_writer(outname) as out:
        for i in range(len(allids)):
            gene_id, isoform_id = allids[i]
            if corrpval[i] < 0.05:
                out.writeColumns(gene_id=gene_id, isoform_id=isoform_id,
                                 delta_usage=round(alldeltas[i], 3),
                                 adj_pval=corrpval[i])

def add_counts_row(genes, isoforms, row, duplicateID):
    "add one counts row to its gene and to the isoform table"
    counts = np.asarray(row.counts, dtype=float)
    if row.gene_id not in genes:
        genes[row.gene_id] = Gene(row.gene_id, np.zeros(len(counts)))
    geneObj = genes[row.gene_id]
    geneObj.exp += counts

    iso = row.isoform_id
    if iso in isoforms:
        duplicateID += 1
        iso = iso + "-" + str(duplicateID)
    isoforms[iso] = Isoform(iso, geneObj, counts)
    return duplicateID

def separate_tables(counts_rows, thresh, samples, a_cols, b_cols, outDir):
    genes, isoforms = dict(), dict()
    duplicateID = 1

    for row in counts_rows:
        duplicateID = add_counts_row(genes, isoforms, row, duplicateID)

    # the two conditions' column indexes are chosen by name in calculate_sig; taking
    # them from the first and last column made the filter depend on column order, and
    # with an order like A,B,B,A it tested one condition twice
    g1Ind = np.asarray(a_cols)
    g2Ind = np.asarray(b_cols)

    # make gene table first
    geneIDs = np.asarray(list(genes.keys()))
    vals = np.asarray([genes[x].exp for x in geneIDs])
    if len(geneIDs) == 0:
        raise FlairInputDataError("no genes parsed from the counts matrix")

    # genes must be expressed in all samples of at least one group
    filteredRows = (np.min(vals[:, g1Ind], axis=1) > thresh) | (np.min(vals[:, g2Ind], axis=1) > thresh)
    filteredGeneVals = vals[filteredRows]
    filteredGeneIDs = geneIDs[filteredRows]
    write_deseq2_counts(outDir + "/filtered_gene_counts_ds2.tsv", samples,
                        filteredGeneIDs, filteredGeneVals)

    # now do isoforms
    isoformIDs = np.asarray(list(isoforms.keys()))
    vals = np.asarray([isoforms[x].exp for x in isoformIDs])
    filteredRows = (np.min(vals[:, g1Ind], axis=1) > thresh) | (np.min(vals[:, g2Ind], axis=1) > thresh)
    filteredIsoVals = vals[filteredRows]
    filteredIsoIDs = isoformIDs[filteredRows]

    write_deseq2_counts(outDir + "/filtered_iso_counts_ds2.tsv", samples,
                        filteredIsoIDs, filteredIsoVals)

    # the DRIMSeq table is every isoform, filtered or not; its serial irow column
    # keeps the index unique, which a repeated isoform id would not
    with drimseq_counts_writer(outDir + "/filtered_iso_counts_drim.tsv", samples) as writer:
        for iso, isoform in isoforms.items():
            writer.writeIsoform(isoform.parent.name, iso, isoform.exp)
    return genes, isoforms


def gene_sample_totals(counts_rows):
    "total counts of each gene in each sample"
    genetosampletotot = {}
    for row in counts_rows:
        counts = [float(x) for x in row.counts]
        if row.gene_id not in genetosampletotot:
            genetosampletotot[row.gene_id] = [0 for x in range(len(counts))]
        genetosampletotot[row.gene_id] = [genetosampletotot[row.gene_id][x] + counts[x]
                                          for x in range(len(counts))]
    return genetosampletotot

def norm_row_by_gene(row, genetot):
    "one isoform's counts scaled so that its gene totals the same in every sample"
    geneavg = sum(genetot) / len(genetot)
    counts = [(float(count) / genetot[i]) * geneavg if genetot[i] > 0 else 0
              for i, count in enumerate(row.counts)]
    return CountsRow(row.gene_id, row.isoform_id, [str(round(x)) for x in counts])

def calc_gene_norm_sig(workdir, counts_rows, sample_columns):
    """
    Counts normalized by gene, written for the record and returned for the t-test.
    This is not a standard normalization method, only used for downstream stats.
    """
    genetosampletotot = gene_sample_totals(counts_rows)
    norm_rows = [norm_row_by_gene(row, genetosampletotot[row.gene_id]) for row in counts_rows]
    write_counts_matrix(workdir + '/counts.normbygene.tsv', sample_columns, norm_rows)
    return norm_rows

def r_compress_option():
    "the R scripts compress their tables when this run does"
    return ['--compress'] if output_compression() else []

def run_deseq2(prefix, workdir, condition_a, condition_b, matrixFile, outDir, formulaMatrixFile):
    # no --batch: neither R script ever read it, and one arbitrary batch label would
    # not have said anything anyway.  Both take the batch column from the formula matrix
    stderr = f"{workdir}/{prefix}.txt"
    try:
        with open(stderr, "w") as stderr_fh:
            pipettor.run(["Rscript", diffExp_deseq2, "--condition_a", condition_a, "--condition_b", condition_b,
                          "--matrix", matrixFile, "--out_dir", outDir,
                          "--prefix", prefix, "--formula", formulaMatrixFile] + r_compress_option(),
                         stderr=stderr_fh)
    except pipettor.ProcessException as exc:
        raise FlairError(f'running {prefix} failed, please check {stderr} for details') from exc

def run_dirmseq(prefix, workdir, threads, condition_a, condition_b, matrixFile, outDir, formulaMatrixFile):
    stderr = f"{workdir}/{prefix}.txt"
    try:
        with open(stderr, "w") as stderr_fh:
            pipettor.run(["Rscript", diffExp_drimseq, "--threads", threads, "--condition_a", condition_a, "--condition_b", condition_b,
                          "--matrix", matrixFile, "--out_dir", outDir,
                          "--prefix", prefix, "--formula", formulaMatrixFile] + r_compress_option(),
                         stderr=stderr_fh)
    except pipettor.ProcessException as exc:
        raise FlairError(f'running {prefix} failed, please check {stderr} for details') from exc


def calculate_sig(*, counts_matrix, output, condition_a, condition_b, min_expression,  # noqa: C901 - FIXME: reduce complexity
                  threads, overwrite_output):
    outDir = output
    quant_table_tsv = counts_matrix
    sFilter = min_expression
    force_dir = overwrite_output

    counts_rows = read_counts_rows(quant_table_tsv)
    if len(counts_rows) == 0:
        raise FlairInputDataError(f"counts matrix {quant_table_tsv} has no isoform rows")

    # Get sample data info
    sample_infos = read_sample_info(quant_table_tsv)
    groups = [si.condition for si in sample_infos]
    batches = [si.batch for si in sample_infos]
    # the column number keeps these unique for R, which joins the formula matrix to
    # the counts matrix by them
    samples = ["%s_%s" % (si.sample_id, num) for num, si in enumerate(sample_infos)]
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

    # the normalized rows keep the counts matrix column order, so the same column
    # indexes apply to them
    norm_rows = calc_gene_norm_sig(workdir, counts_rows, samples)
    get_sig_from_norm_by_gene(compressed_path(outDir + '/isoforms_sig_exp_change_norm_by_gene.tsv'),
                              norm_rows, a_cols, b_cols)

    # Convert count tables to dataframe and update isoform objects.
    genes, isoforms = separate_tables(counts_rows, sFilter, samples, a_cols, b_cols, workdir)

    # DESeq2 models batch when the formula matrix has a batch column, so it is
    # written only when batch varies and is not a linear combination of condition.
    # The no-batch branches used to differ only in a stray leading space in the
    # condition column name
    with_batch = (len(combos) != 2) and (len(set(batches)) > 1)
    formula_rows = [FormulaRow(sample, group, batch)
                    for sample, group, batch in zip(samples, groups, batches)]

    formulaMatrixFile = workdir + "/formula_matrix.tsv"
    write_formula_matrix(formulaMatrixFile, formula_rows, with_batch=with_batch)

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
