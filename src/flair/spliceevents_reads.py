"""Reads corrected to junction chains, per sample and region: their genes, the
fragments of first and last exons among their unspliced reads, and the files
handing them to event calling.  Used by flair spliceevents"""

import os
import os.path as osp
import shutil
import pysam
import logging
from collections import namedtuple
from flair.read_processing import generate_genomic_alignment_read_to_clipping_file
from flair.read_correction import filter_correct_group_reads
from flair.annotation_alignment import identify_good_match_to_annot
from flair.junction_correct import UNKNOWN_STRAND
from flair.gene_identification import get_gene_name_firstpass, get_gene_name_unknown_strand
from flair.interval_index import IntervalIndex
from flair.spliceevents_events import lies_left


def correct_sample_reads(args, region, work_dir, bamfile_name, genome, annots, junction_corrector, read_starts=None):
    """The sample's reads in region, corrected as flair transcriptome corrects them,
    grouped by junction chain: aligned to the annotated transcripts, then their
    junctions corrected, the unspliced reads under no junctions.  Given read_starts,
    a SeqRange, only the reads aligned to start in it.  The alignment's
    intermediates are written in work_dir, removed after"""
    os.makedirs(work_dir, exist_ok=True)
    temp_prefix = osp.join(work_dir, 'reads')
    # the reads as samtools fasta writes them: primary alignments, in the read's own orientation
    with pysam.AlignmentFile(bamfile_name, 'rb') as bam, open(temp_prefix + '.reads.fasta', 'w') as out:
        for read in bam.fetch(region.name, region.start, region.end):
            if not (read.is_secondary or read.is_supplementary) \
                    and (read_starts is None or read_starts.start <= read.reference_start < read_starts.end):
                out.write(f'>{read.query_name}\n{read.get_forward_sequence()}\n')
    sj_to_ends = {}
    with pysam.AlignmentFile(bamfile_name, 'rb') as bam_file:
        generate_genomic_alignment_read_to_clipping_file(temp_prefix, bam_file, region.name, region.start, region.end)
        read_to_transcript = identify_good_match_to_annot(
            temp_prefix, region, annots, genome, bam_file, no_align_to_annot=args.no_align_to_annot,
            directRNA=False, check_splice=True, trust_strand=args.trust_strand, fusion_breakpoints=None)
        filter_correct_group_reads(bam_file=bam_file, region=region,
                                   read_to_annot_transcript=read_to_transcript,
                                   annots=annots, junction_corrector=junction_corrector,
                                   genome=genome,
                                   quality=args.quality, keep_sup=args.keep_supplementary,
                                   sj_to_ends=sj_to_ends,
                                   allow_outside_range=True,
                                   keep_single_exon=True, trust_strand=args.trust_strand,
                                   check_motifs=not args.trust_junctions, read_starts=read_starts)
    shutil.rmtree(work_dir)
    return sj_to_ends


# a read's corrected ends, all event calling needs of it beyond its junction chain
ReadEnds = namedtuple('ReadEnds', ('start', 'end'))


def read_chains_path(partition, sample):
    "the file of a sample's corrected reads in a partition, handed from its correction to the region's event calling"
    return f'{partition.file_prefix}_{sample}.read_chains.tsv'


# a sample's corrected junction chain: its annotated genes, 5' to 3', or () for none,
# its junctions, strand, and its reads' ends.  Or its unspliced reads that are
# fragments of one gene, with no junctions, and fragment 'last' for 3' fragments,
# 'first' for 5' (AnnotatedFragmentExons)
ReadChain = namedtuple('ReadChain', ('genes', 'juncs', 'strand', 'ends', 'fragment'), defaults=('',))


def write_read_chains(path, chrom, chains):
    """Write a sample's corrected reads, ReadChains, a line per junction chain: its
    annotated genes, chromosome, junctions, strand, the ends of each of its reads,
    0-based start-end pairs, and, for fragments, without junctions, which kind.  It's
    written whole or not at all, as a rerun with --keep_intermediate skips the samples
    whose file exists"""
    with open(path + '.tmp', 'w') as out:
        for chain in chains:
            out.write('\t'.join([','.join(chain.genes), chrom, ','.join(f'{s}-{e}' for s, e in chain.juncs), chain.strand,
                                 ','.join(f'{r.start}-{r.end}' for r in chain.ends), chain.fragment]) + '\n')
    os.replace(path + '.tmp', path)


def _ranges(field):
    return tuple(tuple(int(x) for x in pair.split('-')) for pair in field.split(',')) if field else ()


def read_read_chains(path):
    """a sample's corrected reads from write_read_chains, ReadChains; of an earlier
    version's files, without fragments, too"""
    chains = []
    with open(path) as fh:
        for line in fh:
            genes, _, juncs, strand, ends, *fragment = line.rstrip('\n').split('\t')
            chains.append(ReadChain(tuple(genes.split(',')) if genes else (), _ranges(juncs), strand,
                                    [ReadEnds(*pair) for pair in _ranges(ends)], *fragment))
    return chains


def annotated_genes(isoform, annots):
    """A corrected junction chain's annotated genes, 5' to 3', and its strand, as
    flair transcriptome finds them (gene_identification): () if it has none, and a
    strand of None for a chain of unknown strand that no gene gives one"""
    if isoform.strand == UNKNOWN_STRAND:
        genes, _, strand = get_gene_name_unknown_strand(isoform, annots)
        if genes is None:
            return (), None
    else:
        genes, _ = get_gene_name_firstpass(isoform, annots)
        if genes is None:
            return (), isoform.strand
        strand = annots.gene_to_strand[genes[0]]
    # a readthrough's genes in the order it reads them
    return tuple(sorted(genes, key=lambda g: annots.gene_to_exons[g][0][0], reverse=(strand == '-'))), strand


def fragment_chains(reads, fragment_exons, trust_strand):
    """the unspliced reads that are fragments of one gene (AnnotatedFragmentExons), as a ReadChain
    per gene and kind"""
    by_fragment = {}
    for read in reads:
        found = fragment_exons.fragment(read.start, read.end, read.strand if trust_strand else None)
        if found is not None:
            by_fragment.setdefault(found, []).append(ReadEnds(read.start, read.end))
    return [ReadChain((gene,), (), strand, ends, kind) for (kind, gene, strand), ends in sorted(by_fragment.items())]


def get_juncs_single_sample(args, region, path, bamfile_name, genome, junction_corrector, annots, terminal_exons,
                            read_starts=None):
    """correct a sample's reads in region, give their junction chains annotated genes,
    and write them to path, with the unspliced reads that are fragments of one gene's
    first or last exons, terminal_exons; the chains of no annotated gene get novel
    genes in event calling, over all samples"""
    sj_to_ends = correct_sample_reads(args, region, path + '.correction', bamfile_name, genome, annots, junction_corrector,
                                      read_starts)
    fragment_exons = AnnotatedFragmentExons(terminal_exons, annots)
    chains = []
    for (_, juncs), isoform in sj_to_ends.items():
        if len(juncs) > 0:
            genes, strand = annotated_genes(isoform, annots)
            if strand is None:
                logging.debug(f'reads dropped: strand unknown and not given by a gene: {juncs}')
                continue
            chains.append(ReadChain(genes, tuple((j.start, j.end) for j in juncs), strand, isoform.reads))
        else:
            chains.extend(fragment_chains(isoform.reads, fragment_exons, args.trust_strand))
    write_read_chains(path, region.name, chains)


# an exon's first or last exon by kind: 'first' or 'last', strand, start, end, gene
AnnotatedTerminalExon = namedtuple('AnnotatedTerminalExon', ('kind', 'strand', 'start', 'end', 'gene'))


def annotated_terminal_exons(transcripts):
    """chrom -> the first and last exons, AnnotatedTerminalExons, of the spliced transcripts, as
    of GENCODE basic, whose first and last exons are real ones"""
    exons_by_chrom = {}
    for trans in transcripts:
        exons = sorted((exon.start, exon.end) for exon in trans.exons)
        if len(exons) > 1:
            first, last = (exons[0], exons[-1]) if trans.strand == '+' else (exons[-1], exons[0])
            terminal = exons_by_chrom.setdefault(trans.chrom, set())
            terminal.add(AnnotatedTerminalExon('first', trans.strand, *first, trans.gene_id))
            terminal.add(AnnotatedTerminalExon('last', trans.strand, *last, trans.gene_id))
    return {chrom: sorted(exons) for chrom, exons in exons_by_chrom.items()}


def annotated_terminal_splice_sites(terminal_exons):
    """chrom -> strand -> the splice sites, sorted, of terminal_exons' first exons, and
    of their last exons, which mark alternative first and last exons.  By strand, not
    gene, so they don't depend on the reads' genes"""
    afe_ss, ale_ss = {}, {}
    for chrom, exons in terminal_exons.items():
        for sites in (afe_ss, ale_ss):
            sites[chrom] = {'+': set(), '-': set()}
        for exon in exons:
            # the splice site of a first exon is its 3' end, of a last exon its 5' end
            exon_left = lies_left(exon.kind == 'first', exon.strand)
            site = exon.end if exon_left else exon.start
            (afe_ss if exon.kind == 'first' else ale_ss)[chrom][exon.strand].add(site)
    return tuple({chrom: {strand: sorted(strand_sites) for strand, strand_sites in by_strand.items()}
                  for chrom, by_strand in sites.items()} for sites in (afe_ss, ale_ss))


# how far past its splice site, into the intron, the end of an unspliced read in a
# first or last exon may be, as flair transcriptome's --end_window
FRAGMENT_SPLICE_SITE_WINDOW = 100


class AnnotatedFragmentExons:
    """The annotated first and last exons of a region, and the exons of all its genes,
    to tell the unspliced reads that are fragments of a gene's transcripts: in a last
    exon, a 3' fragment, whose end is a transcript's end, or in a first exon, a 5'
    fragment, whose start is a transcript's start"""
    def __init__(self, terminal_exons, annots):
        self.terminal = IntervalIndex()
        for exon in terminal_exons:
            self.terminal.add(exon.start, exon.end, exon)
        # each gene's exons, merged, so the bases a read shares with a gene count once
        self.exons = IntervalIndex()
        for gene, exons in annots.gene_to_exons.items():
            strand = annots.gene_to_strand[gene]
            merged = []
            for start, end in sorted(exons):
                if merged and start <= merged[-1][1]:
                    merged[-1][1] = max(merged[-1][1], end)
                else:
                    merged.append([start, end])
            for start, end in merged:
                self.exons.add(start, end, (gene, strand, start, end))

    def fragment(self, start, end, strand=None):
        """('last' or 'first', gene, strand) of an unspliced read, or None: of the genes
        whose exons overlap the read, the one sharing the most bases with it, if that's
        one gene, and the read's end at its splice site, its 5' end in a last exon or
        its 3' end in a first, lies in the exon or within FRAGMENT_SPLICE_SITE_WINDOW
        past the site, while its other end may run past the transcript's end.  Given
        strand, the read's trusted strand, only the genes on it; otherwise the genes of
        both strands, as an unspliced read's strand isn't known"""
        shared = {}
        for gene, gene_strand, exon_start, exon_end in self.exons.overlap(start, end):
            if strand in (None, gene_strand):
                key = (gene, gene_strand)
                shared[key] = shared.get(key, 0) + min(end, exon_end) - max(start, exon_start)
        if len(shared) == 0:
            return None
        most = max(shared.values())
        best = [key for key, bases in shared.items() if bases == most]
        if len(best) != 1:
            return None
        gene, gene_strand = best[0]
        found = set()
        for exon in self.terminal.overlap(start, end, slack=FRAGMENT_SPLICE_SITE_WINDOW):
            if exon.gene != gene:
                continue
            # a first exon is upstream of its splice site, a last exon downstream
            if lies_left(exon.kind == 'first', exon.strand):  # the splice site is the exon's end
                at_site = exon.start < end <= exon.end + FRAGMENT_SPLICE_SITE_WINDOW
            else:
                at_site = exon.start - FRAGMENT_SPLICE_SITE_WINDOW <= start < exon.end
            if at_site:
                found.add(exon.kind)
        if len(found) != 1:
            return None
        return found.pop(), gene, gene_strand
