#!/usr/bin/env python3
import sys
import csv
import os
import argparse

from flair.gtf_to_bed import get_iso_info

def main():
    parser = argparse.ArgumentParser(description='''identifies the most likely gene id associated with
            each isoform and renames the isoform''')
    parser.add_argument('bed', type=str,
            action='store', help='isoforms in bed format')
    parser.add_argument('gtf', type=str,
            action='store', help='annotated isoform gtf')
    parser.add_argument('outfilename', type=str,
            action='store', help='Name of output file')
    parser.add_argument('--proportion', action='store', default=0.8, dest='proportion_annotated_covered',
            type=float, help='''proportion should be a decimal < 1 specifying the % of an annotated single-exon
            gene a FLAIR isoform has to cover (default=0.8)''')
    parser.add_argument('--annotation_reliant', action='store_true', dest='annotation_reliant',
            help='name all isoforms with -* starting with -0')
    parser.add_argument('--gene_only', action='store_true', dest='gene_only',
            help='only append gene name to read name')
    parser.add_argument('--field_name', action='store', dest='field_name', default='gene_id',
            help='field name to use for gene id, e.g. gene_type or gene_name (default: gene_id)')
    args = parser.parse_args()

    identify_gene_isoform(gtf=args.gtf, field_name=args.field_name, outfilename=args.outfilename,
                          query=args.bed,
                          proportion_annotated_covered=args.proportion_annotated_covered,
                          gene_only=args.gene_only, annotation_reliant=args.annotation_reliant)


def get_junctions(line):
    junctions = set()
    starts = [int(n) + 1 for n in line[20].split(',')[:-1]]
    sizes = [int(n) - 1 for n in line[18].split(',')[:-1]]  # for indexing pupropses
    if len(starts) == 1:
        return
    for b in range(len(starts)-1): # block
        junctions.add((starts[b]+sizes[b], starts[b+1]))
    return junctions


def get_junctions_bed12(line):
    junctions = set()
    chrstart = int(line[1])
    starts = [int(n) + chrstart + 1 for n in line[11].split(',')[:-1]]
    sizes = [int(n) - 1 for n in line[10].split(',')[:-1]]
    if len(starts) == 1:
        return
    for b in range(len(starts)-1): # block
        junctions.add((starts[b]+sizes[b], starts[b+1]))
    return junctions


def bin_search(query, data):
    """ Query is a coordinate interval. Binary search for the query in sorted data,
    which is a list of coordinates. Finishes when an overlapping value of query and
    data exists and returns the index in data. """
    i = int(round(len(data)/2))  # binary search prep
    lower, upper = 0, len(data)
    while True:
        if upper - lower < 2:  # stop condition but not necessarily found
            break
        if data[i][1] < query[0]:
            lower = i
            i = int(round((i + upper)/2))
        elif data[i][0] > query[1]:
            upper = i
            i = int(round((lower + i)/2))
        else:  # found
            break
    return i


def overlapping_bases(coords0, coords1):
    """ complete coverage of coords0 by coords1, and coords0 can be tol larger.
    if coords0 is contained by coords1, then return the number of
    overlapping basepairs """
    if coords0[1] > coords1[0] and coords1[1] > coords0[0]:
        return min(coords1[1], coords0[1]) - max(coords1[0], coords0[0])
    return


def update_tn_dicts(chrom, junctions, prev_transcript, prev_exon, junc_to_tn,
        tn_to_juncs, all_se):
    if chrom not in junc_to_tn:
        junc_to_tn[chrom] = {}
        tn_to_juncs[chrom] = {}
        all_se[chrom] = []
    if not junctions or len(junctions) == 0:
        all_se[chrom].append(prev_exon)
    else:
        tn_to_juncs[chrom][prev_transcript] = junctions
        for j in junctions:
            if j not in junc_to_tn[chrom]:
                junc_to_tn[chrom][j] = set()
            junc_to_tn[chrom][j].add(prev_transcript)
    return junc_to_tn, tn_to_juncs, all_se


def update_gene_dicts(chrom, j, gene, junctions, gene_unique_juncs, junc_to_gene):
    junctions.add(j)
    if gene not in gene_unique_juncs:
        gene_unique_juncs[gene] = set()
    gene_unique_juncs[gene].add(j)
    if j not in junc_to_gene[chrom]:
        junc_to_gene[chrom][j] = set()
    junc_to_gene[chrom][j].add(gene)
    return junctions, gene_unique_juncs, junc_to_gene


def identify_gene_isoform(gtf, outfilename, query, field_name='gene_id', proportion_annotated_covered=0.8,
                          gene_only=False, annotation_reliant=False):
    junc_to_tn = {}  # matches intron to transcript; chrom: {intron: [transcripts], ... }
    tn_to_juncs = {}  # matches transcript to intron; i.e. chrom: {transcript_name: (junction1, junction2), ... }
    all_se = {}  # all single exon genes
    junc_to_gene = {}  # matches a splice junction (i.e. an intron) to gene name
    gene_unique_juncs = {}  # matches a gene to its set of unique splice junctions

    if gtf:
        iso_to_info, iso_to_exons, iso_to_cds = get_iso_info(gtf, adjustpos=False)

        for transcript in iso_to_info:
            chrom, strand, gene = iso_to_info[transcript]
            exons = sorted(iso_to_exons[transcript])
            if chrom not in junc_to_gene:
                junc_to_gene[chrom] = {}
            junctions = set()
            for i in range(len(exons)-1):
                junctions, gene_unique_juncs, junc_to_gene = update_gene_dicts(chrom, (exons[i][1], exons[i+1][0]), gene,
                                                                               junctions,gene_unique_juncs, junc_to_gene)

            junc_to_tn, tn_to_juncs, all_se = update_tn_dicts(chrom, junctions, transcript, (exons[-1][0], exons[-1][1], gene),
                                                              junc_to_tn, tn_to_juncs, all_se)

        for chrom in all_se:
            all_se[chrom] = sorted(list(all_se[chrom]), key=lambda x: x[0])

    name_counts = {}  # to avoid redundant names
    with open(outfilename, 'wt') as outfile:
        writer = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep)
        for line in open(query):
            line = line.rstrip().split('\t')
            junctions = get_junctions_bed12(line)
            chrom, name, start, end = line[0], line[3], int(line[1]), int(line[2])
            if ';' in name:
                name = name[:name.find(';')]

            if chrom not in junc_to_tn:  # chrom not in reference file
                if name not in name_counts:
                    name_counts[name] = 0
                else:
                    name_counts[name] += 1
                    name = name + '-' + str(name_counts[name])
                noref = chrom + ':' + str(start)[:-3] + '000'
                newname = name + '_' + noref
                line[3] = newname
                writer.writerow(line)
                continue

            gene_hits = {}
            se_gene_tiebreaker = {}
            if not junctions:
                exon = (start, end)
                i = bin_search(exon, all_se[chrom])
                for e in all_se[chrom][i-2:i+2]:
                    overlap = overlapping_bases(exon, e)
                    if overlap:
                        proportion = float(overlap)/(exon[1]-exon[0])  # base coverage of long-read isoform by the annotated isoform
                        proportion2 = float(overlap)/(e[1]-e[0])  # base coverage of the annotated isoform by the long-read isoform
                        if proportion > 0.5 and proportion2 > proportion_annotated_covered:
                            if e[2] in gene_hits: # gene name
                                if proportion <= gene_hits[e[2]]:
                                    continue
                            gene_hits[e[2]] = proportion
                            se_gene_tiebreaker[e[2]] = proportion2
            else:
                for j in junctions:
                    if j in junc_to_gene[chrom]:
                        for gene in junc_to_gene[chrom][j]:
                            if gene not in gene_hits:
                                gene_hits[gene] = 0
                            gene_hits[gene] += 1  # gene name, number of junctions this isoform shares with this gene

            if not gene_hits:  # gene name will just be a chromosome locus
                gene = chrom + ':' + str(start)[:-3] + '000'
            else:  # gene name will be whichever gene the entry has more shared junctions with
                genes = sorted(gene_hits.items(), key=lambda x: x[1])  # sort by number of junctions shared with gene
                if len(genes) > 1 and genes[-1][1] == genes[-2][1]: # tie, break by gene size
                    genes = sorted(genes, key=lambda x: x[0])
                    genes = sorted(genes, key=lambda x: x[1])
                    if not junctions:
                        g = genes[-1], se_gene_tiebreaker[genes[-1][0]]
                        for i in reversed(range(len(genes)-1)):
                            if genes[i][1] == g[0][1]:
                                if se_gene_tiebreaker[genes[i][0]] > g[1]:
                                    g = genes[i], se_gene_tiebreaker[genes[i][0]]
                            else:
                                break
                        genes[-1] = g[0]
                    else:
                        g = genes[-1], len(gene_unique_juncs[genes[-1][0]])
                        for i in reversed(range(len(genes)-1)):
                            if genes[i][1] == g[0][1]:
                                if len(gene_unique_juncs[genes[i][0]]) < g[1]:
                                    g = genes[i], len(gene_unique_juncs[genes[i][0]])
                            else:
                                break
                        genes[-1] = g[0]
                gene = genes[-1][0]

            transcript = ''
            if junctions:
                matches = set()
                for j in junctions:
                    if j in junc_to_tn[chrom]:
                        matches.update(junc_to_tn[chrom][j])
                for t in sorted(list(matches)):
                    if tn_to_juncs[chrom][t] == junctions:
                        transcript = t  # annotated transcript identified
                        break
            name = transcript if transcript and not gene_only else name

            if gene_only:
                newname = name + '_' + gene
            elif name not in name_counts:
                name_counts[name] = 0
                if annotation_reliant:
                    newname = name + '-0_' + gene
                else:
                    newname = name + '_' + gene
            else:
                name_counts[name] += 1
                newname = name + '-' + str(name_counts[name]) + '_' + gene

            line[3] = newname
            line[8] = "20,47,181" if transcript else "232,142,23" ##blue if annotated, orange if novel
            if line[9] == '1': line[8] = "242,208,17" #yellow if monoexon
            writer.writerow(line)

if __name__ == "__main__":
    main()
