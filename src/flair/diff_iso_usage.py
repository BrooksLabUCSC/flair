#!/usr/bin/env python3

import sys
import csv
import os
import scipy.stats as sps
from flair import FlairInputDataError
from flair.pycbio.sys import cli


def parse_args():
    desc = """Calculates the usage of each isoform as a fraction of the total expression
    of the gene and compares this between samples."""

    parser = cli.ArgumentParserExtras(description=desc)
    parser.add_argument('counts_matrix_tsv',
                        help='counts matrix TSV from flair-quantify')
    parser.add_argument('colname1',
                        help='the name of the column of the first sample')
    parser.add_argument('colname2',
                        help='the name of the column of the second sample')
    parser.add_argument('outfile',
                        help='output filename containing the p-value associated with differential '
                        'isoform usage for each isoform')
    return parser.parse_args()

def split_iso_gene(iso_gene):
    if '_chr' in iso_gene:
        iso = iso_gene[:iso_gene.rfind('_chr')]
        gene = iso_gene[iso_gene.rfind('_chr')+1:]
    elif '_XM' in iso_gene:
        iso = iso_gene[:iso_gene.rfind('_XM')]
        gene = iso_gene[iso_gene.rfind('_XM')+1:]
    elif '_XR' in iso_gene:
        iso = iso_gene[:iso_gene.rfind('_XR')]
        gene = iso_gene[iso_gene.rfind('_XR')+1:]
    elif '_NM' in iso_gene:
        iso = iso_gene[:iso_gene.rfind('_NM')]
        gene = iso_gene[iso_gene.rfind('_NM')+1:]
    elif '_NR' in iso_gene:
        iso = iso_gene[:iso_gene.rfind('_NR')]
        gene = iso_gene[iso_gene.rfind('_NR')+1:]
    else:
        iso = iso_gene[:iso_gene.rfind('_')]
        gene = iso_gene[iso_gene.rfind('_')+1:]
    return iso, gene


def diff_iso_usage(counts_matrix_tsv, colname1, colname2, outfilename):
    counts_matrix_fh = open(counts_matrix_tsv)
    header = counts_matrix_fh.readline().rstrip().split('\t')

    if colname1 in header:
        col1 = header.index(colname1)
    else:
        raise FlairInputDataError('Could not find {} in {}\n'.format(colname1, ' '.join(header)))
    if colname2 in header:
        col2 = header.index(colname2)
    else:
        raise FlairInputDataError('Could not find {} in {}\n'.format(colname2, ' '.join(header)))

    counts = {}
    for line in counts_matrix_fh:
        line = line.rstrip().split('\t')
        iso_gene, count1, count2 = line[0], float(line[col1]), float(line[col2])
        if '_' not in iso_gene:
            raise FlairInputDataError('Incorrect isoform names: Please run identify_annotated_gene first so that \n'
                             'isoforms can be grouped by their parent genes\n')
        iso, gene = split_iso_gene(iso_gene)
        if gene not in counts:
            counts[gene] = {}
        counts[gene][iso] = [count1, count2]

    with open(outfilename, 'wt') as outfile:
        writer = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep)
        writer.writerow(['geneID', 'isoID', 'fisher_pval', 'this_iso_sample1_count', 'this_iso_sample2_count', 'other_isos_sample1_count', 'other_isos_sample2_count'])
        geneordered = sorted(counts.keys())
        for gene in geneordered:
            generes = []
            for iso in counts[gene]:
                thesecounts = counts[gene][iso]
                othercounts = [0, 0]
                ctable = [thesecounts, othercounts]
                for iso_ in counts[gene]:
                    if iso_ == iso:
                        continue
                    othercounts[0] += counts[gene][iso_][0]
                    othercounts[1] += counts[gene][iso_][1]
                if thesecounts[0] + othercounts[0] == 0 or thesecounts[1] + othercounts[1] == 0:  # do not test this isoform if no gene exp in one sample
                    continue
                ctable[1] = othercounts
                if ctable[0][0] + ctable[1][0] == 0 or ctable[0][1] + ctable[1][1] == 0 or not sum(ctable[1]):
                    continue
                generes += [[gene, iso, sps.fisher_exact(ctable)[1]] + ctable[0] + ctable[1]]
            if not generes:
                writer.writerow([gene, iso, 'NA'] + ctable[0] + ctable[1])
                continue

            for res in generes:
                writer.writerow(res)


def main():
    args = parse_args()
    with cli.ErrorHandler():
        diff_iso_usage(args.counts_matrix_tsv, args.colname1, args.colname2, args.outfile)

if __name__ == '__main__':
    main()
