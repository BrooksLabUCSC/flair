"""Shared filter + junction-correction + grouping for FLAIR pipelines."""

import logging

from flair.annotation_precheck import splice_sites_cleanly_aligned
from flair.iso_gene_id import split_iso_gene
from flair.isoform_data import Junc, ReadRec
from flair.read_processing import should_process_read, add_corrected_read_to_groups


def _add_corrected_spliced_read(read, readrec, aligned_juncs, genome, sj_to_ends):
    """Flag whether correction moved the read's junctions from its alignment's
    introns and, if not, whether the alignment is clean around them, then group it"""
    juncs = tuple((j.start, j.end) for j in sorted(readrec.juncs))
    readrec.junctions_moved = juncs != tuple((j.start, j.end) for j in sorted(aligned_juncs))
    if not readrec.junctions_moved and genome is not None:
        readrec.clean_splice_sites = splice_sites_cleanly_aligned(read, juncs, genome)
    add_corrected_read_to_groups(readrec, sj_to_ends)


def _correct_and_group_read(read, *, read_to_annot_transcript, annots,
                            junction_corrector, sj_to_ends, genome,
                            keep_single_exon, trust_strand, check_motifs):
    """Correct a single read's splice junctions and add it to sj_to_ends groups.

    Spliced and single-exon reads are fundamentally different:
    - Spliced: junctions corrected from annotation or intron support, strand from correction
    - Single-exon: no correction, strand resolved later in group_se_by_overlap

    Given the genome and check_motifs, a spliced read whose junctions are only
    supported by introns that are weak strand evidence (unannotated, without a
    GT-AG motif) gets UNKNOWN_STRAND, to be resolved by gene identification.

    With trust_strand, a spliced read keeps the strand of its alignment: an
    annotated transcript match on the other strand is not used, and junctions are
    corrected only from introns on the read's strand.

    keep_single_exon=False drops single-exon reads (both reads without juncs
    and reads matching annotated single-exon transcripts).
    """
    readrec = ReadRec.from_read(read, genome=genome)
    aligned_juncs = readrec.juncs

    # annotated spliced: correct junctions and strand from annotation.  A read spliced
    # inside an annotated exon, as at an intron in a 3' UTR, aligns to the transcript
    # with that intron as a deletion, which allow_UTR_indels tolerates in a terminal
    # exon, and the match gives it fewer junctions than its alignment has.  Such a
    # read is corrected from intron support, keeping a real intron, and only failing
    # that, from the annotation, the junction more likely an alignment artifact.  A
    # short exon the genome alignment missed gives the match more junctions, which
    # is what it corrects
    if read.query_name in read_to_annot_transcript:
        tid, startindex, startdist, endindex, enddist = read_to_annot_transcript[read.query_name]
        # the composite id heuristic, shared with bed_to_gtf and the rest, rather than
        # a bare split on the last '_', which lands in the wrong place for any gene id
        # containing one
        transcript, gene = split_iso_gene(tid)
        exons = annots.transcript_to_exons[(transcript, gene)]
        annot_juncs = [(exons[x].end, exons[x + 1].start) for x in range(len(exons) - 1)]
        if len(annot_juncs) > 0 and not (trust_strand and annots.gene_to_strand[gene] != readrec.strand):
            newstart = annot_juncs[startindex][0] - startdist
            newend = annot_juncs[endindex][1] + enddist
            juncs = tuple([Junc(x[0], x[1]) for x in annot_juncs[startindex:endindex + 1]])
            from_introns = (len(juncs) < len(readrec.juncs)
                            and junction_corrector.correct_readrec(readrec, trust_strand, genome if check_motifs else None))
            if from_introns:
                logging.debug(f"annotation match would drop read junctions, corrected from intron support: {readrec.name}")
            else:
                readrec.correct_from_annotation(newstart, newend, annots.gene_to_strand[gene], juncs)
            _add_corrected_spliced_read(read, readrec, aligned_juncs, genome, sj_to_ends)
            return

    # unannotated spliced: correct junctions and strand from intron support
    if readrec.juncs:
        if junction_corrector.correct_readrec(readrec, trust_strand, genome if check_motifs else None):
            _add_corrected_spliced_read(read, readrec, aligned_juncs, genome, sj_to_ends)
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
                               keep_single_exon=True, trust_strand=False, check_motifs=True):
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
                                    keep_single_exon=keep_single_exon,
                                    trust_strand=trust_strand,
                                    check_motifs=check_motifs)
