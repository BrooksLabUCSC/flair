#!/usr/bin/env python3
import sys
import csv
import os
from flair import FlairInputDataError

def main():
    try:
        isoforms = sys.argv[1]
        s = float(sys.argv[2])
        outfilename = sys.argv[3]
    except:
        raise FlairInputDataError('usage: filter_isoforms_by_proportion_of_gene_expr.py isoforms support_percentage outfile')

    if s >= 1:
        raise FlairInputDataError('Support percentage should be a decimal e.g. 0.1 for 10%')

    filter_isoforms_by_proportion_of_gene_expr(isoforms=isoforms, outfilename=outfilename,
                                        support=s,)


def split_iso_gene(iso_gene):
    if '_chr' in iso_gene:
        splitchar = '_chr'
    elif '_XM' in iso_gene:
        splitchar = '_XM'
    elif '_XR' in iso_gene:
        splitchar = '_XR'
    elif '_NM' in iso_gene:
        splitchar = '_NM'
    elif '_NR' in iso_gene:
        splitchar = '_NR'
    elif '_R2_' in iso_gene:
        splitchar = '_R2_'
    elif '_NC_' in iso_gene:
        splitchar = '_NC_'
    else:
        splitchar = '_'
    iso = iso_gene[:iso_gene.rfind(splitchar)]
    gene = iso_gene[iso_gene.rfind(splitchar)+1:]
    return iso, gene

def filter_isoforms_by_proportion_of_gene_expr(isoforms, outfilename, support):
    genes = dict()
    isoFH = open(isoforms, 'r')
    for line in isoFH:
        line = line.rstrip().split('\t')
        name = line[3]
        iso, gene = split_iso_gene(name)
        if gene not in genes:
            genes[gene] = set()
        genes[gene].add(tuple(line))

    with open(outfilename, 'wt') as outfile:
        writer = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep)
        for gene in genes:
            gene_total = float(sum([float(iso[-1]) for iso in genes[gene]]))
            if gene_total == 0:
                continue
            for iso in genes[gene]:
                if float(iso[-1])/gene_total >= support:
                    writer.writerow(iso)

if __name__ == "__main__":
    main()
