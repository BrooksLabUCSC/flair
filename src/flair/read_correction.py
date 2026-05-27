"""Shared filter + junction-correction + grouping for FLAIR pipelines."""

import logging

from flair.isoform_data import Junc, ReadRec
from flair.read_processing import should_process_read, add_corrected_read_to_groups


def _correct_and_group_read(read, *, read_to_annot_transcript, annots,
                            junction_corrector, sj_to_ends, genome,
                            keep_single_exon):
    """Correct a single read's splice junctions and add it to sj_to_ends groups.

    Spliced and single-exon reads are fundamentally different:
    - Spliced: junctions corrected from annotation or intron support, strand from correction
    - Single-exon: no correction, strand resolved later in group_se_by_overlap

    keep_single_exon=False drops single-exon reads (both reads without juncs
    and reads matching annotated single-exon transcripts).
    """
    readrec = ReadRec.from_read(read, genome=genome)

    # annotated spliced: correct junctions and strand from annotation
    if read.query_name in read_to_annot_transcript:
        # FIXME more id assumptions
        tid, startindex, startdist, endindex, enddist = read_to_annot_transcript[read.query_name]
        transcript = '_'.join(tid.split('_')[:-1])
        gene = tid.split('_')[-1]
        exons = annots.transcript_to_exons[(transcript, gene)]
        annot_juncs = [(exons[x].end, exons[x + 1].start) for x in range(len(exons) - 1)]
        if len(annot_juncs) > 0:
            newstart = annot_juncs[startindex][0] - startdist
            newend = annot_juncs[endindex][1] + enddist
            juncs = tuple([Junc(x[0], x[1]) for x in annot_juncs[startindex:endindex + 1]])
            readrec.correct_from_annotation(newstart, newend, annots.gene_to_strand[gene], juncs)
            add_corrected_read_to_groups(readrec, sj_to_ends)
            return

    # unannotated spliced: correct junctions and strand from intron support
    if readrec.juncs:
        if junction_corrector.correct_readrec(readrec):
            add_corrected_read_to_groups(readrec, sj_to_ends)
        else:
            logging.debug(f"read dropped: junction correction failed: {readrec.name}")
        return

    # single-exon: no correction, strand resolved later in group_se_by_overlap
    if keep_single_exon:
        add_corrected_read_to_groups(readrec, sj_to_ends)
    else:
        logging.debug(f"read dropped: single-exon: {readrec.name}")


def filter_correct_group_reads(*, bam_file, region, read_to_annot_transcript,
                               annots, junction_corrector, genome,
                               quality, keep_sup, sj_to_ends,
                               allow_secondary=False, allow_outside_range=False,
                               keep_single_exon=True):
    """Filter reads, correct splice junctions, and group by junction chain.
    sj_to_ends is mutated in place."""
    for read in bam_file.fetch(region.name, region.start, region.end):
        if should_process_read(read, region, quality, keep_sup,
                               allow_secondary, allow_outside_range):
            _correct_and_group_read(read,
                                    read_to_annot_transcript=read_to_annot_transcript,
                                    annots=annots,
                                    junction_corrector=junction_corrector,
                                    sj_to_ends=sj_to_ends,
                                    genome=genome,
                                    keep_single_exon=keep_single_exon)
