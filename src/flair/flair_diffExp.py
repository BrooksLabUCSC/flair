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


import logging
import os
import os.path as osp
from collections import Counter
from statistics import median, mean
import pipettor

from flair import FlairInputDataError, FlairToolError
from flair.io_utils import keep_or_remove, make_run_temp_dir
from flair.run_checks import (check_option_ranges, count_records, precheck,
                              report_result_count)
from flair.conditions import condition_column_indexes, select_condition_pair
from flair.counts_matrix_tsv import (QUANTIFY_COUNTS, CountsRow, check_isoform_ids_unique,
                                     read_grouped_counts, write_counts_matrix)
from flair.deseq2_counts_tsv import write_deseq2_counts
from flair.drimseq_counts_tsv import drimseq_counts_writer
from flair.drimseq_input import check_feature_ids_unique
from flair.formula_matrix_tsv import FormulaRow, write_formula_matrix
from flair.iso_usage_change_tsv import iso_usage_change_writer
from flair.output_compress import compressed_path, output_compression
from flair.sample_info_tsv import read_sample_info

os.environ['OPENBLAS_NUM_THREADS'] = '1'
import numpy as np  # noqa: E402


# an isoform is tested only when the two conditions' median counts differ by more
# than this
MIN_MEDIAN_COUNT_DIFF = 3

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
    alphac_sidak = 1 - np.power((1. - alphaf), 1. / ntests)
    alphac_bonf = alphaf / float(ntests)

    alphac_sidak_all = 1 - np.power((1. - alphaf),
                                    1. / np.arange(ntests, 0, -1))
    notreject = pvals > alphac_sidak_all
    del alphac_sidak_all

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
    return reject_, pvals_corrected_, alphac_sidak, alphac_bonf


def get_gene_to_counts(rows):
    gene_to_tot_counts = {}
    for row in rows:
        counts = [int(x) for x in row.counts]
        if row.gene_id not in gene_to_tot_counts:
            gene_to_tot_counts[row.gene_id] = [0 for x in range(len(counts))]
        gene_to_tot_counts[row.gene_id] = [gene_to_tot_counts[row.gene_id][x] + counts[x] for x in range(len(counts))]
    return gene_to_tot_counts


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

    # Only test if median difference is large enough
    if abs(deltaval) <= MIN_MEDIAN_COUNT_DIFF:
        return None
    # Two-sample t-test between WT and VAR counts
    # ranksums() wast too strict for small replicates
    return deltausage, ttest_ind(wtcounts, varcounts).pvalue

def do_mtc_ttest(rows, gene_to_tot_counts, ref_cols, test_cols, unreported, counts_format):
    """the t-test of each row but those of unreported, whose counts are in their
    gene's total all the same"""
    # imported here rather than at module scope; scipy takes 0.7s to import and
    # this is the only use of it, which would be paid by every flair command
    from scipy.stats import ttest_ind
    allids, allpval, alldeltas = [], [], []
    for row in rows:
        if row.isoform_id in unreported:
            continue
        tested = row_ttest(ttest_ind, row, gene_to_tot_counts[row.gene_id], ref_cols, test_cols)
        if tested is not None:
            deltausage, pval = tested
            allids.append((row.gene_id, row.isoform_id))
            allpval.append(pval)
            alldeltas.append(deltausage)

    if len(allpval) == 0:
        raise FlairInputDataError(f"no {counts_format.row_noun}'s median count differs between the two "
                                  f"conditions by more than {MIN_MEDIAN_COUNT_DIFF} reads, so there is "
                                  "nothing to test; check that the conditions compared are the ones meant "
                                  "and that the matrix contains read counts rather than counts per million")

    # Apply multiple-testing correction.  This is Holm-Sidak, not Benjamini-Hochberg
    # as this comment used to read: the adjusted values control the family-wise error
    # rate, not the false discovery rate

    corrpval = list(multipletests(allpval)[1])
    return allids, alldeltas, corrpval


def get_sig_from_norm_by_gene(outname, norm_rows, ref_cols, test_cols, unreported, counts_format):
    """
    This function runs t-tests with multiple testing correction on isoform counts normalized by gene
    This method essentially does differential isoform usage testing, but accounts for differences in gene expression
    This is better for detecting novel transcripts than DRIM-seq
    """

    gene_to_tot_counts = get_gene_to_counts(norm_rows)
    allids, alldeltas, corrpval = do_mtc_ttest(norm_rows, gene_to_tot_counts, ref_cols,
                                               test_cols, unreported, counts_format)

    with iso_usage_change_writer(outname) as out:
        for i in range(len(allids)):
            gene_id, isoform_id = allids[i]
            if corrpval[i] < 0.05:
                out.writeColumns(gene_id=gene_id, isoform_id=isoform_id,
                                 delta_usage=round(alldeltas[i], 3),
                                 adj_pval=corrpval[i])

def add_counts_row(genes, isoforms, row):
    """add one counts row to its gene and to the isoform table.  The ids do not
    repeat, which check_isoform_ids_unique has reported if they do"""
    counts = np.asarray(row.counts, dtype=float)
    if row.gene_id not in genes:
        genes[row.gene_id] = Gene(row.gene_id, np.zeros(len(counts)))
    gene_obj = genes[row.gene_id]
    gene_obj.exp += counts
    isoforms[row.isoform_id] = Isoform(row.isoform_id, gene_obj, counts)

def separate_tables(counts_rows, thresh, samples, a_cols, b_cols, out_dir):
    genes, isoforms = dict(), dict()
    for row in counts_rows:
        add_counts_row(genes, isoforms, row)

    # the two conditions' column indexes are chosen by name in calculate_sig; taking
    # them from the first and last column made the filter depend on column order, and
    # with an order like A,B,B,A it tested one condition twice
    g1_ind = np.asarray(a_cols)
    g2_ind = np.asarray(b_cols)

    # make gene table first
    gene_ids = np.asarray(list(genes.keys()))
    vals = np.asarray([genes[x].exp for x in gene_ids])

    # a row is tested when the reads of at least one condition, summed over its
    # samples, reach the threshold, which is what --min_expression's help states.
    # It was the minimum over the samples of a condition, so one low sample in each
    # condition dropped the row however many reads the conditions held
    filtered_rows = (np.sum(vals[:, g1_ind], axis=1) >= thresh) | (np.sum(vals[:, g2_ind], axis=1) >= thresh)
    filtered_gene_vals = vals[filtered_rows]
    filtered_gene_ids = gene_ids[filtered_rows]
    write_deseq2_counts(out_dir + "/filtered_gene_counts_ds2.tsv", samples,
                        filtered_gene_ids, filtered_gene_vals)

    # now do isoforms
    isoform_ids = np.asarray(list(isoforms.keys()))
    vals = np.asarray([isoforms[x].exp for x in isoform_ids])
    filtered_rows = (np.sum(vals[:, g1_ind], axis=1) >= thresh) | (np.sum(vals[:, g2_ind], axis=1) >= thresh)
    filtered_iso_vals = vals[filtered_rows]
    filtered_iso_ids = isoform_ids[filtered_rows]

    write_deseq2_counts(out_dir + "/filtered_iso_counts_ds2.tsv", samples,
                        filtered_iso_ids, filtered_iso_vals)

    # the DRIMSeq table is every isoform, filtered or not; its serial irow column
    # keeps the index unique, which a repeated isoform id would not
    with drimseq_counts_writer(out_dir + "/filtered_iso_counts_drim.tsv", samples) as writer:
        for iso, isoform in isoforms.items():
            writer.writeIsoform(isoform.parent.name, iso, isoform.exp)
    return genes, isoforms


def gene_sample_totals(counts_rows):
    "total counts of each gene in each sample"
    gene_to_sample_to_tot = {}
    for row in counts_rows:
        counts = [float(x) for x in row.counts]
        if row.gene_id not in gene_to_sample_to_tot:
            gene_to_sample_to_tot[row.gene_id] = [0 for x in range(len(counts))]
        gene_to_sample_to_tot[row.gene_id] = [gene_to_sample_to_tot[row.gene_id][x] + counts[x]
                                              for x in range(len(counts))]
    return gene_to_sample_to_tot

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
    gene_to_sample_to_tot = gene_sample_totals(counts_rows)
    norm_rows = [norm_row_by_gene(row, gene_to_sample_to_tot[row.gene_id]) for row in counts_rows]
    write_counts_matrix(workdir + '/counts.normbygene.tsv', sample_columns, norm_rows)
    return norm_rows

def r_compress_option():
    "the R scripts compress their tables when this run does"
    return ['--compress'] if output_compression() else []

def run_deseq2(prefix, id_column, workdir, condition_a, condition_b, matrix_file, output, formula_matrix_file):
    # no --batch: neither R script ever read it, and one arbitrary batch label would
    # not have said anything anyway.  Both take the batch column from the formula matrix
    stderr = f"{workdir}/{prefix}.txt"
    try:
        with open(stderr, "w") as stderr_fh:
            pipettor.run(["Rscript", diffExp_deseq2, "--condition_a", condition_a, "--condition_b", condition_b,
                          "--matrix", matrix_file, "--out_prefix", output, "--workdir", workdir,
                          "--prefix", prefix, "--id_column", id_column,
                          "--formula", formula_matrix_file] + r_compress_option(),
                         stderr=stderr_fh)
    except pipettor.ProcessException as exc:
        raise FlairToolError(f'running {prefix} failed, please check {stderr} for details') from exc

def run_dirmseq(prefix, workdir, threads, condition_a, condition_b, matrix_file, output, formula_matrix_file,
                unreported_file):
    check_feature_ids_unique(matrix_file)
    stderr = f"{workdir}/{prefix}.txt"
    try:
        with open(stderr, "w") as stderr_fh:
            pipettor.run(["Rscript", diffExp_drimseq, "--threads", threads, "--condition_a", condition_a, "--condition_b", condition_b,
                          "--matrix", matrix_file, "--out_prefix", output, "--workdir", workdir,
                          "--prefix", prefix, "--formula", formula_matrix_file,
                          "--unreported", unreported_file] + r_compress_option(),
                         stderr=stderr_fh)
    except pipettor.ProcessException as exc:
        raise FlairToolError(f'running {prefix} failed, please check {stderr} for details') from exc


def calculate_sig(*, counts_matrix, output, condition_a, condition_b, min_expression,  # noqa: C901 - FIXME: reduce complexity
                  threads, temp_dir=None, keep_intermediate=False):
    quant_table_tsv = counts_matrix
    s_filter = min_expression

    grouped = read_grouped_counts(quant_table_tsv)
    counts_rows = grouped.rows
    if len(counts_rows) == 0:
        raise FlairInputDataError(f"counts matrix {quant_table_tsv} has no isoform rows")
    check_isoform_ids_unique(quant_table_tsv, counts_rows)

    # Get sample data info
    sample_infos = read_sample_info(quant_table_tsv)
    groups = [si.condition for si in sample_infos]
    batches = [si.batch for si in sample_infos]
    # R joins the formula matrix to the counts matrix by these, and the result
    # tables are written with them, so they are the matrix's own column names; a
    # matrix of two columns of one name is rejected when it is read
    samples = [si.sample_id for si in sample_infos]
    combos = set([(groups.index(x), batches.index(y)) for x, y in zip(groups, batches)])

    condition_a, condition_b = select_condition_pair(groups, condition_a, condition_b,
                                                     quant_table_tsv)
    a_cols = condition_column_indexes(groups, condition_a)
    b_cols = condition_column_indexes(groups, condition_b)

    group_counts = Counter(groups)
    if len(list(group_counts.keys())) != 2:
        raise FlairInputDataError(f"flair diffexp compares two conditions, but the sample info of "
                                  f"{quant_table_tsv} names {len(group_counts)}: "
                                  f"{', '.join(sorted(group_counts))}.  Cut the counts matrix and its "
                                  f"sample info down to the two conditions to compare")
    elif min(list(group_counts.values())) < 3:
        raise FlairInputDataError("flair diffexp needs at least 3 samples of each condition; "
                                  "with fewer, compare two samples with the diff_iso_usage program")
    elif set(groups).intersection(set(batches)):
        raise FlairInputDataError(f"a condition and a batch have the same name, "
                                  f"{', '.join(sorted(set(groups) & set(batches)))}, which the model "
                                  f"formula cannot tell apart; rename the batch in the sample info of "
                                  f"{quant_table_tsv}")
    elif sum([1 if x.isdigit() else 0 for x in groups]) > 0 or sum([1 if x.isdigit() else 0 for x in batches]) > 0:
        raise FlairInputDataError(f"condition and batch names must not be numbers, which are read as "
                                  f"numbers rather than as groups; name them in the sample info of "
                                  f"{quant_table_tsv}, batch1 rather than 1")

    # the unfiltered tables and the R scripts' other intermediates, in a directory
    # of this run's own, removed at the end unless they were asked for
    workdir = make_run_temp_dir(output, temp_dir)

    with keep_or_remove(workdir, keep_intermediate):
        # the normalized rows keep the counts matrix column order, so the same column
        # indexes apply to them
        norm_rows = calc_gene_norm_sig(workdir, counts_rows, samples)
        get_sig_from_norm_by_gene(compressed_path(output + '.isoforms_sig_exp_change_norm_by_gene.tsv'),
                                  norm_rows, a_cols, b_cols, grouped.unreported, grouped.format)

        # Convert count tables to dataframe and update isoform objects.
        genes, isoforms = separate_tables(counts_rows, s_filter, samples, a_cols, b_cols, workdir)

        # DESeq2 models batch when the formula matrix has a batch column, so it is
        # written only when batch varies and is not a linear combination of condition.
        # The no-batch branches used to differ only in a stray leading space in the
        # condition column name
        with_batch = (len(combos) != 2) and (len(set(batches)) > 1)
        formula_rows = [FormulaRow(sample, group, batch)
                        for sample, group, batch in zip(samples, groups, batches)]

        formula_matrix_file = workdir + "/formula_matrix.tsv"
        write_formula_matrix(formula_matrix_file, formula_rows, with_batch=with_batch)

        iso_matrix_file = workdir + "/filtered_iso_counts_ds2.tsv"
        gene_matrix_file = workdir + "/filtered_gene_counts_ds2.tsv"
        drim_matrix_file = workdir + "/filtered_iso_counts_drim.tsv"

        # DESeq2 genes & isoforms.  An event's reads are not the expression of a gene, so
        # of flair spliceevents counts only usage is tested
        expression_tested = grouped.format is QUANTIFY_COUNTS
        if expression_tested:
            run_deseq2("genes_deseq2", 'gene_id', workdir, condition_a, condition_b, gene_matrix_file,
                       output, formula_matrix_file)
            run_deseq2("isoforms_deseq2", 'isoform_id', workdir, condition_a, condition_b, iso_matrix_file,
                       output, formula_matrix_file)
        else:
            logging.info(f'{quant_table_tsv} is of {grouped.format.program}: testing usage only, '
                         'not expression with DESeq2')

        # DIRMSeq, its results of the rows that are not tested left out
        unreported_file = workdir + "/unreported_features.txt"
        with open(unreported_file, 'w') as fh:
            fh.writelines(feature_id + '\n' for feature_id in sorted(grouped.unreported))
        run_dirmseq("isoforms_drimseq", workdir, threads, condition_a, condition_b, drim_matrix_file, output,
                    formula_matrix_file, unreported_file)
        report_results(output, condition_a, condition_b, grouped.format, expression_tested)


# each filtered table a run writes, as (the noun for one row, its plural, the file)
def _result_tables(condition_a, condition_b, row_noun, expression_tested):
    tables = [(f'{row_noun} with differential usage', f'{row_noun}s with differential usage',
               f'isoforms_drimseq_{condition_a}_v_{condition_b}.tsv'),
              (f'{row_noun} whose usage changes by the t-test', f'{row_noun}s whose usage changes by the t-test',
               'isoforms_sig_exp_change_norm_by_gene.tsv')]
    if expression_tested:
        tables = [('differentially expressed gene', None, f'genes_deseq2_{condition_a}_v_{condition_b}.tsv'),
                  ('differentially expressed isoform', None,
                   f'isoforms_deseq2_{condition_a}_v_{condition_b}.tsv')] + tables
    return tables


def report_results(output, condition_a, condition_b, counts_format, expression_tested):
    """What each test reported, which is otherwise a table the user has to open: a
    test that reported nothing and a test that was not run look the same in the
    output directory."""
    for what, plural_what, name in _result_tables(condition_a, condition_b, counts_format.row_noun,
                                                  expression_tested):
        path = compressed_path(f'{output}.{name}')
        report_result_count(what, max(count_records(path) - 1, 0), plural_what=plural_what,
                            remedy='check that the conditions compared are the ones meant; the tables before '
                            'filtering are in the run directory, which --keep_intermediate keeps')


# the programs a run executes, checked before it starts
PROGRAMS = ('Rscript',)


def add_subparser(subparsers):
    desc = "Differential expression and differential usage analysis"
    parser = subparsers.add_parser('diffexp', help="Differential expression and usage analysis",
                                   description=desc)
    required = parser.add_argument_group('required named arguments')
    required.add_argument('--counts_matrix', type=str, required=True,
                          help='tab-delimited isoform count matrix from flair quantify, or the counts of flair '
                          'spliceevents (OUTPUT.spliceevents.counts.tsv), whose usage only is tested.  Conditions '
                          'and batches are read from the sample info TSV beside it')
    required.add_argument('-o', '--output', type=str, required=True,
                          help='output file name base: each table and plot is written as '
                               '<output>.<name>, as every other subcommand writes its outputs')
    parser.add_argument('-t', '--threads', type=int, default=12,
                        help='number of threads for parallel DRIMSeq (default: %(default)s)')
    parser.add_argument('--min_expression', type=int, default=10,
                        help='read count expression threshold: a gene or isoform is tested when the '
                             'reads of at least one condition, summed over its samples, reach this many, '
                             'and is left out when both conditions have fewer (default: %(default)s)')
    parser.add_argument('--condition_a', default='',
                        help='the reference condition, as named in the counts matrix columns; '
                             'fold changes are reported for --condition_b relative to this. '
                             'With neither this nor --condition_b given, the two conditions are '
                             'taken in sorted order')
    parser.add_argument('--condition_b', default='',
                        help='the condition compared against --condition_a')
    parser.add_argument('--temp_dir',
                        help='directory for temporary files; each run makes its own directory, named for the output, in it (default: $TMPDIR or the system temporary directory)')
    parser.add_argument('--keep_intermediate', action='store_true',
                        help="keep the run's own directory, which contains the tables before filtering and the "
                             "R scripts' diagnostics, and whose path the run reports")
    parser.set_defaults(entry=diffexp_cmd)

def diffexp_cmd(args):
    check_option_ranges(args, count=('min_expression',), positive=('threads',))
    precheck(files={'counts matrix': args.counts_matrix}, programs=PROGRAMS)
    calculate_sig(counts_matrix=args.counts_matrix, output=args.output,
                  condition_a=args.condition_a, condition_b=args.condition_b,
                  min_expression=args.min_expression, threads=args.threads,
                  temp_dir=args.temp_dir, keep_intermediate=args.keep_intermediate)
