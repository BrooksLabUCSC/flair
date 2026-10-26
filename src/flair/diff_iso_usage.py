#!/usr/bin/env python3

import argparse
from flair import FlairInputDataError
from flair.iso_usage_tsv import iso_usage_writer
from flair.counts_matrix_tsv import read_grouped_counts
from flair.pycbio.sys import cli


def build_parser():
    desc = """Calculates the usage of each isoform as a fraction of the total expression
    of the gene and compares this between samples.  It reads a flair quantify counts
    matrix, or the counts of flair spliceevents (each alternative of an event) or flair
    variantquant (the reads showing and not showing each allele)."""

    parser = argparse.ArgumentParser(prog='diff_iso_usage', description=desc)
    parser.add_argument('counts_matrix_tsv',
                        help='counts matrix TSV from flair quantify, OUTPUT.spliceevents.counts.tsv from '
                        'flair spliceevents, or OUTPUT.var.counts.tsv from flair variantquant')
    parser.add_argument('colname1',
                        help='the name of the column of the first sample')
    parser.add_argument('colname2',
                        help='the name of the column of the second sample')
    parser.add_argument('outfile',
                        help='output filename containing the p-value associated with differential '
                        'isoform usage for each isoform, and that p-value adjusted for the isoforms tested '
                        '(Benjamini-Hochberg)')
    parser.add_argument('--min_expression', type=int, default=0,
                        help="an isoform with fewer than this many reads in the two samples together is not tested, "
                        "its p-values NA, though its reads count toward its gene's other isoforms.  Filtering on "
                        "the total, not on each sample, leaves the adjusted p-values valid, and testing fewer "
                        "isoforms that could never be significant gives more power (default: %(default)s)")
    return parser

def sample_column_index(sample_columns, colname, counts_matrix_tsv):
    "index of a named sample column within a row's counts"
    if colname not in sample_columns:
        raise FlairInputDataError(
            f"{counts_matrix_tsv} has no sample column named {colname}; it has: "
            f"{' '.join(sample_columns)}")
    return sample_columns.index(colname)

def iso_psi(iso_counts, other_counts):
    "each sample's PSI for one isoform, and the change, NA where it cannot be taken"
    s1_psi, s2_psi, delta_psi = 'NA', 'NA', 'NA'
    if other_counts[0] + iso_counts[0] > 0:
        s1_psi = round(iso_counts[0] / (other_counts[0] + iso_counts[0]), 3)
    if other_counts[1] + iso_counts[1] > 0:
        s2_psi = round(iso_counts[1] / (other_counts[1] + iso_counts[1]), 3)
    if s1_psi != 'NA' and s2_psi != 'NA':
        delta_psi = round(s2_psi - s1_psi, 3)
    return [s1_psi, s2_psi, delta_psi]

def gene_total_counts(gene_counts):
    "the gene's counts in each sample, of all its isoforms"
    return [sum(counts[0] for counts in gene_counts.values()), sum(counts[1] for counts in gene_counts.values())]

def other_isoform_counts(gene_total, iso_counts):
    """the gene's counts without this isoform: its total less the isoform's, so a gene
    of many isoforms, or an event of many alternatives, is not summed again for each"""
    return [gene_total[0] - iso_counts[0], gene_total[1] - iso_counts[1]]

def gene_usage_rows(gene, gene_counts, colname1, colname2, unreported=frozenset(), min_expression=0):
    """one row per isoform of a gene, but those of unreported, which count toward the
    other isoforms of the rest.  An isoform of fewer than min_expression reads in the
    two samples together is not tested"""
    import scipy.stats as sps
    rows = []
    gene_total = gene_total_counts(gene_counts)
    for iso, iso_counts in gene_counts.items():
        if iso in unreported:
            continue
        other_counts = other_isoform_counts(gene_total, iso_counts)
        ctable = [iso_counts, other_counts]
        # an isoform with no expression of its gene in one sample cannot be tested
        if (iso_counts[0] + other_counts[0] == 0 or iso_counts[1] + other_counts[1] == 0
                or sum(iso_counts) == 0 or sum(other_counts) == 0):
            rows.append([gene, iso, 'NA'] + iso_counts + other_counts + ['NA', 'NA', 'NA'])
        elif sum(iso_counts) < min_expression:
            rows.append([gene, iso, 'NA'] + iso_counts + other_counts + iso_psi(iso_counts, other_counts))
        else:
            rows.append([gene, iso, sps.fisher_exact(ctable)[1]] + iso_counts +
                        other_counts + iso_psi(iso_counts, other_counts))
    return rows

def add_adjusted_pvalues(rows):
    """put after each row's p-value its Benjamini-Hochberg adjusted one, over the
    rows tested, NA for a row of none"""
    import scipy.stats as sps
    tested = [i for i, row in enumerate(rows) if row[2] != 'NA']
    adjusted = dict(zip(tested, map(float, sps.false_discovery_control([rows[i][2] for i in tested])))) if tested else {}
    for i, row in enumerate(rows):
        row.insert(3, adjusted.get(i, 'NA'))
    return rows

def diff_iso_usage(counts_matrix_tsv, colname1, colname2, outfilename,  # noqa: C901 - FIXME: reduce complexity
                   min_expression=0):
    grouped = read_grouped_counts(counts_matrix_tsv)
    col1 = sample_column_index(grouped.sample_columns, colname1, counts_matrix_tsv)
    col2 = sample_column_index(grouped.sample_columns, colname2, counts_matrix_tsv)

    counts = {}
    for row in grouped.rows:
        if row.gene_id not in counts:
            counts[row.gene_id] = {}
        counts[row.gene_id][row.isoform_id] = [row.counts[col1], row.counts[col2]]

    rows = [row for gene in sorted(counts.keys())
            for row in gene_usage_rows(gene, counts[gene], colname1, colname2, grouped.unreported, min_expression)]
    with iso_usage_writer(outfilename, colname1, colname2) as writer:
        for row in add_adjusted_pvalues(rows):
            writer.writeRow(row)


def main():
    args = cli.parseArgsWithLogging(build_parser())
    with cli.ErrorHandler():
        diff_iso_usage(args.counts_matrix_tsv, args.colname1, args.colname2, args.outfile,
                       min_expression=args.min_expression)


if __name__ == '__main__':
    main()
