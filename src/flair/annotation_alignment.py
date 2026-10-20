"""
Aligning reads to transcripts, shared by flair transcriptome and flair spliceevents:
the reference of a region's annotated transcripts, with each gene's terminal exons
extended to its transcripts' furthest ends, and the alignment of reads to it, whose
good matches correct the reads' junctions and ends.
"""
import logging
from flair import thread_share
from flair.annotation_precheck import reads_skipping_annotation_alignment
from flair.count_sam_transcripts import run_count_sam_transcripts
from flair.isoform_data import Exon, get_bed_exons_from_exons, get_sequence_for_exons
from flair.terminal_exon_ends import TerminalExonEnds


####
# terminal exon normalization
####
def _exon_bounds(exons):
    return [(e.start, e.end) for e in exons]

def max_terminal_exons_ends_from_annots(annots):
    max_terminal_exons_ends = TerminalExonEnds()
    for transcript_id, gene_id, strand in annots.transcripts:
        exons = annots.transcript_to_exons[(transcript_id, gene_id)]
        max_terminal_exons_ends.add_transcript(gene_id, strand, _exon_bounds(exons))
    return max_terminal_exons_ends

####
# transcriptome reference
####
def normalize_gene_terminal_exons(max_terminal_exons_ends, gene_id, strand, exons, chrom_len):
    "updates terminal exons ends"
    start, end = max_terminal_exons_ends.normalized_ends(gene_id, strand, _exon_bounds(exons), chrom_len)
    exons[0] = Exon(start, exons[0].end)
    exons[-1] = Exon(exons[-1].start, end)
    return exons

def generate_transcriptome_reference_transcript(strand, transcript_to_strand, transcript_id, gene_id, annots, max_terminal_exons_ends,
                                                transcript_to_new_exons, chrom, genome, annot_bed_fh, annot_fa_fh, annot_uniqueseq_fh):
    transcript_to_strand[(transcript_id, gene_id)] = strand
    exons = list(annots.transcript_to_exons[(transcript_id, gene_id)])
    assert isinstance(exons[0], Exon)  # FIXME tmp debugging
    if len(exons) > 1:
        normalize_gene_terminal_exons(max_terminal_exons_ends, gene_id, strand, exons, genome.get_reference_length(chrom))
        transcript_to_new_exons[(transcript_id, gene_id)] = tuple(exons)
    exons = tuple(exons)
    start, end = exons[0].start, exons[-1].end

    # FIXME: duplicated code
    exon_starts, exon_sizes = get_bed_exons_from_exons(exons, start)
    # FIXME: duplicated use BED class,
    bed_line = [chrom, start, end, transcript_id + '_' + gene_id, '.', strand, start, end, '0', len(exons),
                ','.join([str(x) for x in exon_sizes]), ','.join([str(x) for x in exon_starts])]
    trans_seq = get_sequence_for_exons(genome, chrom, strand, exons)
    annot_bed_fh.write('\t'.join([str(x) for x in bed_line]) + '\n')
    annot_fa_fh.write('>' + transcript_id + '_' + gene_id + '\n')
    annot_fa_fh.write(''.join(trans_seq) + '\n')

def generate_transcriptome_reference_guts(annots, chrom, genome, annot_bed_fh, annot_fa_fh, annot_uniqueseq_fh):
    transcript_to_strand = {}
    transcript_to_new_exons = {}
    max_terminal_exons_ends = max_terminal_exons_ends_from_annots(annots)

    for transcript_id, gene_id, strand in annots.transcripts:
        generate_transcriptome_reference_transcript(strand, transcript_to_strand, transcript_id, gene_id, annots, max_terminal_exons_ends,
                                                    transcript_to_new_exons, chrom, genome, annot_bed_fh, annot_fa_fh, annot_uniqueseq_fh)
    return transcript_to_strand, transcript_to_new_exons

def generate_transcriptome_reference(temp_prefix, annots, chrom, genome):
    with (open(temp_prefix + '.annotated_transcripts.bed', 'w') as annot_bed_fh,
          open(temp_prefix + '.annotated_transcripts.fa', 'w') as annot_fa_fh,
          open(temp_prefix + '.annotated_transcripts_uniquebound.txt', 'w') as annot_uniqueseq_fh):
        return generate_transcriptome_reference_guts(annots, chrom, genome, annot_bed_fh, annot_fa_fh, annot_uniqueseq_fh)


def _write_reads_subset(in_fasta, out_fasta, skip):
    "write the reads of a FASTA not named in skip; returns the number written"
    nwritten, keep = 0, False
    with open(in_fasta) as in_fh, open(out_fasta, 'w') as out_fh:
        for line in in_fh:
            if line.startswith('>'):
                keep = line[1:].split()[0] not in skip
                nwritten += keep
            if keep:
                out_fh.write(line)
    return nwritten

def identify_good_match_to_annot(temp_prefix, region, annots, genome, bam_file, *, no_align_to_annot, directRNA,
                                 check_splice, trust_strand, fusion_breakpoints, windows=None):
    """read -> (annotated transcript, start splice site index and distance, end splice
    site index and distance) of the reads in region that align well to an annotated
    transcript, of the reads in temp_prefix.reads.fasta, their genomic clipping in
    temp_prefix.reads.genomicclipping.txt.  Reads whose correction the alignment can't
    change aren't aligned.  Given windows, ReadWindows, only the reads in them are
    checked for needing the alignment, as only they are in the fasta."""
    chrom = region.name
    read_to_transcript = {}
    if not no_align_to_annot and len(annots.transcripts) > 0:
        # reads whose correction aligning to the annotation can't change go straight
        # to junction correction
        skip = reads_skipping_annotation_alignment(bam_file, region, annots, genome, windows)
        annot_reads = temp_prefix + '.annotreads.fasta'
        nreads = _write_reads_subset(temp_prefix + '.reads.fasta', annot_reads, skip)
        logging.info(f'aligning {nreads} reads to the annotated transcripts; {len(skip)} need no alignment')
        if nreads == 0:
            return read_to_transcript
        # this part generates the fasta file for the annotation
        transcript_to_strand, transcript_to_new_exons = \
            generate_transcriptome_reference(temp_prefix, annots, chrom, genome)
        # FIXME: make a TSV
        clipping_file = temp_prefix + '.reads.genomicclipping.txt'
        transcriptome_align_and_count(annot_reads,
                                      temp_prefix + '.annotated_transcripts.fa',
                                      temp_prefix + '.annotated_transcripts.bed',
                                      temp_prefix + '.matchannot.counts.txt',
                                      None, True,
                                      clipping_file,
                                      temp_prefix + '.annotated_transcripts_uniquebound.txt',
                                      directRNA=directRNA, no_stringent=False, check_splice=check_splice,
                                      trust_strand=trust_strand, fusion_breakpoints=fusion_breakpoints)
        for line in open(temp_prefix + '.matchannot.ends.tsv'):
            line = line.rstrip().split('\t')
            read, transcript = line[:2]
            start_sj_index, start_sj_dist, start_tend_dist, end_sj_index, end_sj_dist, end_tend_dist = [int(x) if x != 'None' else None for x in line[2:]]
            # is not None: the value was converted on the line above, so the old
            # comparison with the string 'None' was always true
            if start_sj_index is not None and end_sj_index is not None:  # not a single exon transcript
                read_to_transcript[read] = (transcript, start_sj_index, start_sj_dist, end_sj_index, end_sj_dist)
    return read_to_transcript


# the secondary alignments of each read to the transcripts that minimap2 reports for
# count_sam_transcripts to choose from (-N)
SECONDARY_ALIGNMENTS = 4


def transcriptome_align_and_count(input_reads, align_ref_fasta, ref_bed, output_name, map_file, is_annot, clipping_file, unique_bound,  # noqa: C901 - FIXME: reduce complexity
                                  *, directRNA, no_stringent, check_splice, trust_strand, fusion_breakpoints,
                                  secondary_alignments=SECONDARY_ALIGNMENTS):
    """Align reads to transcripts with minimap2 and assign each to its best
    transcript with count_sam_transcripts: to the annotated transcripts (is_annot),
    or to flair transcriptome's firstpass isoforms.  minimap2 reports
    secondary_alignments other alignments of each read, which must reach the
    transcript it fits among many alike, as a truncated read among the isoforms
    containing its own"""
    # minimap (results are piped into count_sam_transcripts.py)
    # '--split-prefix', 'minimap2transcriptomeindex', doesn't work with MD tag
    if isinstance(input_reads, str):
        input_reads = [input_reads]
    mm2_cmd = ['minimap2', '-a', '-N', str(secondary_alignments), '--MD'] + [align_ref_fasta] + input_reads

    # FIXME add in step to filter out chimeric reads here
    # FIXME really need to go in and check on how count_sam_transcripts is working
    trimmedreads = clipping_file or None
    generate_map = map_file or None
    output_endpos = output_name.split('.counts.txt')[0] + '.ends.tsv'  # if (args.output_endpos or is_annot) else None)
    stringent = (not is_annot) and (not no_stringent)
    unique_bound_path = unique_bound if unique_bound and (not no_stringent or is_annot) else None

    # minimap2 borrows threads left idle in the thread budget while it runs
    with thread_share.program_threads() as mm2_threads:
        run_count_sam_transcripts(
            mm2_cmd=mm2_cmd[:1] + ['-t', str(mm2_threads)] + mm2_cmd[1:],
            output=output_name,
            trimmedreads=trimmedreads,
            generate_map=generate_map,
            output_endpos=output_endpos,
            stringent=stringent,
            allow_UTR_indels=not directRNA,
            transcript_fasta=align_ref_fasta,
            check_splice=check_splice,
            isoforms=ref_bed,
            trust_ends=False,
            unique_bound=unique_bound_path,
            fusion_breakpoints=fusion_breakpoints,
            # reads come from samtools fasta in their sequenced orientation, and
            # transcript sequences are sense, so an antisense read aligns reversed
            stranded=trust_strand,
            # in the annotation pass, an alignment clipped where the read runs past the
            # transcript's end would cut the read's corrected ends short; without its
            # match the read keeps its own ends through junction correction.  The final
            # realignment is to padded transcripts, and an alignment clipped past the
            # read's genomic clipping is a read whose sequence the transcript doesn't
            # have; it also gives the reported ends, which would stop at the
            # transcript's end
            no_extra_clipping=True)
