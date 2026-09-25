#!/usr/bin/env python3

########################################################################
# File: predictProductivity.py
#  executable: predictProductivity.py
# Purpose:
#
#
# Author: Cameron M. Soulette
# History:      cms 10/09/2018 Created
#                   06/15/26 Rewritten by Colette Felton
#
########################################################################


# import pipettor
# import os
# import argparse
# from flair import FlairInputDataError
# from flair.gtf_io import gtf_record_parser, GtfAttrsSet
# from flair.pycbio.hgdata.bed import Bed, BedReader
# from flair.pycbio.sys import fileOps
# from flair.bed_to_gtf import bed_to_gtf
# from flair.flair_bed import FlairBed
from flair.isoform_data import get_reverse_complement, translate as translate_seq

STOP_CODON_SEQS = set(['TAA', 'TGA', 'TAG'])
MAX_DIST_FROM_EXON_EDGE_FOR_PTC = 55

def get_annot_start_codons(transcript, gene_to_cds_starts):
    my_annot_starts = set()
    for gene in transcript.gene.gene_desc:
        if gene in gene_to_cds_starts:
            my_annot_starts.update(gene_to_cds_starts[gene])
    return my_annot_starts

def identify_start_exon_index(my_exons, annot_start):
    start_exon_index = None
    for i, e in enumerate(my_exons):
        if e.start <= annot_start <= e.end:
            start_exon_index = i
            break
    return start_exon_index

def calc_transcript_rel_start_pos(annot_start, exon_sizes, my_exons, start_exon_index, strand):
    if strand == '+':
        rel_start = sum(exon_sizes[:start_exon_index]) + (annot_start - my_exons[start_exon_index].start)
    else:
        rel_start = sum(exon_sizes[start_exon_index + 1:]) + (my_exons[start_exon_index].end - annot_start)
    return rel_start

def calc_stop_codon_pos(seq_from_start):
    """Offset of the first stop codon, and whether one was found.  The loop variable
    used to be returned even when the loop never ran, which is an UnboundLocalError on
    an empty sequence, as happens when the annotated start codon is at the last base."""
    stop_codon_pos = 0
    for stop_codon_pos in range(0, len(seq_from_start), 3):
        if seq_from_start[stop_codon_pos:stop_codon_pos + 3] in STOP_CODON_SEQS:
            return True, stop_codon_pos
    return False, stop_codon_pos

def calc_ptc(exon_sizes, orf_end_pos, ref_transcript_id, transcript_to_nmd_except):
    is_ptc = True
    if ref_transcript_id in transcript_to_nmd_except and transcript_to_nmd_except[ref_transcript_id]:
        is_ptc = False
    elif orf_end_pos > sum(exon_sizes[:-2]) and orf_end_pos > sum(exon_sizes[:-1]) - MAX_DIST_FROM_EXON_EDGE_FOR_PTC:
        is_ptc = False
    return "PTC" if is_ptc else "PRO"

def calc_genomic_end_pos(my_exons, exon_sizes, orf_end_pos, strand):
    curr_start = 0
    genomic_end_pos = None
    for i in range(len(my_exons)):
        if curr_start <= orf_end_pos < curr_start + exon_sizes[i]:
            if strand == '+':
                genomic_end_pos = my_exons[i].start + (orf_end_pos - curr_start)
            else:
                genomic_end_pos = my_exons[i].end - (orf_end_pos - curr_start)
        curr_start += exon_sizes[i]
    return genomic_end_pos

def identify_prod_from_start(orfs, annot_start, rel_start, my_exons, exon_sizes, my_seq, strand, ref_transcript_id,
                             transcript_to_nmd_except):
    five_UTR, seq_from_start = my_seq[:rel_start], my_seq[rel_start:].upper()
    stop_reached, stop_codon_pos = calc_stop_codon_pos(seq_from_start)
    if not stop_reached:
        transcript_end = my_exons[-1].end if strand == "+" else my_exons[0].start
        orfs.append(["NST", annot_start, transcript_end, len(five_UTR) + stop_codon_pos - rel_start, rel_start, None])
    else:
        orf_end_pos = len(five_UTR) + stop_codon_pos + 3
        # order exon sizes from 5' to 3'
        if strand == '-':
            exon_sizes = exon_sizes[::-1]
            my_exons = my_exons[::-1]
        ptc = calc_ptc(exon_sizes, orf_end_pos, ref_transcript_id, transcript_to_nmd_except)
        genomic_end_pos = calc_genomic_end_pos(my_exons, exon_sizes, orf_end_pos, strand)
        aaseq = None if ptc == 'PTC' else translate(seq_from_start[:orf_end_pos - rel_start])
        orfs.append([ptc, annot_start, genomic_end_pos, orf_end_pos - rel_start, rel_start, aaseq])

def identify_best_orf_from_starts(transcript, my_annot_starts, my_seq, transcript_to_nmd_except):
    my_exons = sorted(transcript.exons)  # sorted from left to right on genome
    exon_sizes = [x.end - x.start for x in my_exons]
    orfs = []
    for annot_start in my_annot_starts:
        start_exon_index = identify_start_exon_index(my_exons, annot_start)
        if start_exon_index is not None:
            rel_start = calc_transcript_rel_start_pos(annot_start, exon_sizes, my_exons, start_exon_index, transcript.strand)
            identify_prod_from_start(orfs, annot_start, rel_start, my_exons, exon_sizes, my_seq, transcript.strand,
                                     transcript.ref_transcript_id, transcript_to_nmd_except)
    if len(orfs) == 0:
        orfs.append(["NGO", sorted(my_exons)[0].start, sorted(my_exons)[0].start, 0, 0, None])
    orfs.sort(key=lambda x: x[3], reverse=True)
    return orfs[0]

def translate_from_bed(bed_transcript, genome):
    if bed_transcript.productivity == 'PRO':
        seq_to_translate = ''
        my_exons = sorted(bed_transcript.blocks)  # sorted from left to right on genome
        for e in my_exons:
            if e.start <= bed_transcript.thickStart <= e.end:
                seq_to_translate += genome.fetch(bed_transcript.chrom, bed_transcript.thickStart, e.end)
            elif bed_transcript.thickStart < e.start and e.end < bed_transcript.thickEnd:
                seq_to_translate += genome.fetch(bed_transcript.chrom, e.start, e.end)
            elif e.start <= bed_transcript.thickEnd <= e.end:
                seq_to_translate += genome.fetch(bed_transcript.chrom, e.start, bed_transcript.thickEnd)
        if bed_transcript.strand == '-':
            seq_to_translate = get_reverse_complement(seq_to_translate)
        # print(seq_to_translate[-10:], translate(seq_to_translate)[-10:])
        return translate(seq_to_translate)
    return None

def predict_prod_temp(transcript, start_codon_count, gene_to_cds_starts, transcript_to_nmd_except, genome):
    my_prod = None
    thickStart, thickEnd = None, None
    my_aaseq = None
    if start_codon_count > 0:  # annotations exist and contain start codons
        my_annot_starts = get_annot_start_codons(transcript, gene_to_cds_starts)
        my_orf = identify_best_orf_from_starts(transcript, my_annot_starts, transcript.get_sequence(genome), transcript_to_nmd_except)
        if transcript.strand == '+':
            thickStart, thickEnd = my_orf[1], my_orf[2]
        else:
            thickStart, thickEnd = my_orf[2], my_orf[1]
        my_prod = my_orf[0]
        my_aaseq = my_orf[5]
    return thickStart, thickEnd, my_prod, my_aaseq

def translate(seq):
    """Protein for a coding sequence, without the trailing stop.  The table and the
    codon lookup live in isoform_data; this had its own identical copy."""
    return translate_seq(seq)[:-1]
