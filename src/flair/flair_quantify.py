#! /usr/bin/env python3

import os
import pysam
import shutil
import itertools
import logging
from bisect import bisect_right
from flair import FlairInputDataError
from flair.counts_matrix import SampleInfo, sample_info_path, write_sample_info
from flair.io_utils import make_run_temp_dir
from flair import thread_share
from flair.pycbio.hgdata.bed import BedReader, BedBlock
from flair.flair_bed import FlairBed
from flair.terminal_exon_ends import GeneTerminalExonEnds
from flair.count_sam_transcripts import run_count_sam_transcripts, TRUST_ENDS_WINDOW
from flair.annotation_precheck import read_introns, splice_sites_cleanly_aligned
from flair.interval_index import IntervalIndex
from flair.identify_prelim_fusions import needs_transcriptome
from flair.counts_to_tpm import counts_to_tpm
from flair.isoform_data import get_reverse_complement
import multiprocessing as mp


os.environ['OPENBLAS_NUM_THREADS'] = '1'

# batches of loci per thread, for the pool to balance them
BATCHES_PER_THREAD = 8
# a read clipped this much at an end may be a fusion's, as flair fusion's default --min_fragment_size
FUSION_MIN_CLIP = 40

# the direct assignments' lookups over all genes (direct_assignment_chains,
# single_exon_index), and per sample the reads counted for fusions (count_fusions),
# set by quantify before the batch tasks fork, which share them
_chains = {}
_single_exon_index = {}
_fusion_reads = []

def add_subparser(subparsers):
    desc = "Quantify the expression level of isoforms across samples"
    parser = subparsers.add_parser('quantify', help="Quantify isoform expression",
                                   description=desc)
    required = parser.add_argument_group('required named arguments')
    required.add_argument('--manifest', type=str, required=True,
                          help='Tab delimited file containing sample id, condition, batch, and the path to '
                               "that sample's reads aligned to the genome as a sorted, indexed BAM")
    required.add_argument('-g', '--genome', type=str, required=True,
                          help='FastA of genome')
    required.add_argument('--isoform_bed', type=str, required=True,
                          help='isoform .bed file')
    parser.add_argument('-o', '--output', type=str, default='flair.quantify',
                        help='output file name base for FLAIR quantify (default: %(default)s)')
    parser.add_argument('-t', '--threads', type=int, default=4,
                        help='minimap2 number of threads (default: %(default)s)')
    parser.add_argument('--tpm', action='store_true',
                        help='convert counts matrix to transcripts per million and output as a separate file named <output>.tpm.tsv')
    parser.add_argument('--trust_ends', action='store_true',
                        help='specify if reads are generated from a long read method with minimal fragmentation')
    parser.add_argument('--generate_map', action='store_true',
                        help='create read-to-isoform assignment files for each sample')
    parser.add_argument('--with_gene', action='store_true',
                        help='output lines with isoform_gene')
    parser.add_argument('--normalize_ends', action='store_true',
                        help="normalize transcript ends; recommended when not using --trust_ends and "
                             "different transcript ends do not matter")
    parser.add_argument('--temp_dir',
                        help='directory for temporary files; each run makes its own directory, named for the '
                             'output, in it, with the reads and isoforms of each batch of genes realigned '
                             '(default: $TMPDIR or the system temporary directory)')
    parser.add_argument('--keep_intermediate', action='store_true',
                        help="keep intermediate and temporary files for debugging, in the run's directory in "
                             '--temp_dir, which must be given')
    # parser.add_argument('--stringent', action='store_true',
    #                     help="""Supporting reads must cover 80 percent of their isoform and extend at least 25 nt into the
    #                     first and last exons. If those exons are themselves shorter than 25 nt, the requirement becomes
    #                     'must start within 4 nt from the start" or "must end within 4 nt from the end" """)
    # parser.add_argument('--check_splice', action='store_true',
    #                     help="""enforce coverage of 4 out of 6 bp around each splice site and no
    #                     insertions greater than 3 bp at the splice site""")
    parser.set_defaults(entry=quantify_cmd)

def quantify_cmd(args):
    if args.keep_intermediate and args.temp_dir is None:
        raise FlairInputDataError('--keep_intermediate requires --temp_dir, the directory to keep them in')
    check_input_files(manifest=args.manifest, genome=args.genome, isoform_bed=args.isoform_bed)
    quantify(manifest=args.manifest, genome=args.genome, isoform_bed=args.isoform_bed,
             output=args.output, threads=args.threads,
             tpm=args.tpm, trust_ends=args.trust_ends, generate_map=args.generate_map,
             with_gene=args.with_gene, normalize_ends=args.normalize_ends,
             temp_dir=args.temp_dir, keep_intermediate=args.keep_intermediate)

def check_input_files(*, manifest, genome, isoform_bed):
    if isoform_bed.endswith('.psl'):
        raise FlairInputDataError('FLAIR no longer accepts PSL input, convert it with psl_to_bed: ' + isoform_bed)
    for what, path in (('isoform models bed', isoform_bed), ('genome fasta', genome),
                       ('manifest', manifest)):
        if not os.path.exists(path):
            raise FlairInputDataError(f'{what} file does not exist: {path}')

def load_manifest(manifest):
    """the manifest rows.  The id, condition and batch fields may contain anything:
    they are written to the sample info file as their own columns rather than being
    joined into the counts column names."""
    sample_data = []
    for line in open(manifest):
        cols = line.rstrip().split('\t')
        if len(cols) < 4:
            raise FlairInputDataError(f'expected 4 tab-separated columns in {manifest}, got {len(cols)}')

        sample, group, batch, readFile = cols
        if not os.path.exists(readFile):
            raise FlairInputDataError(f'reads file does not exist: {readFile}')

        sample_data.append((sample, group, batch, readFile))
    return sample_data

class GeneData():
    def __init__(self, gene_id, chrom, strand):
        self.gene_id = gene_id
        self.chrom = chrom
        self.strand = strand
        self.left_bound = None
        self.right_bound = None
        self.isoform_beds = []
        # self.isoform_fas = []
        self.terminal_exon_ends = GeneTerminalExonEnds(gene_id)

    def add_isoform_bed(self, isoform):
        self.isoform_beds.append(isoform)
        if self.left_bound is None or isoform.chromStart < self.left_bound:
            self.left_bound = isoform.chromStart
        if self.right_bound is None or isoform.chromEnd > self.right_bound:
            self.right_bound = isoform.chromEnd
        if len(isoform.blocks) > 1:
            self.terminal_exon_ends.add(isoform.blocks[0].start, isoform.blocks[0].end,
                                        isoform.blocks[-1].start, isoform.blocks[-1].end)

def load_isoform_data(isoform_bed):
    """The isoforms by gene, gene id -> GeneData, and fusion isoforms, name -> their
    pieces in pos_in_fusion order.  A fusion's pieces share its name and gene id, on
    any chromosomes and strands, so they aren't in a gene (count_fusions)."""
    gene_data, fusions = {}, {}
    names = {}
    for bed in BedReader(isoform_bed, bedClass=FlairBed):
        if bed.pos_in_fusion is not None:
            fusions.setdefault(bed.name, []).append(bed)
            continue
        if bed.gene_id not in gene_data:
            gene_data[bed.gene_id] = GeneData(bed.gene_id, bed.chrom, bed.strand)
        elif bed.strand != gene_data[bed.gene_id].strand:
            raise FlairInputDataError(f"gene id '{bed.gene_id}' has isoforms on both strands, "
                                      f"'{gene_data[bed.gene_id].strand}' and '{bed.strand}'; give each strand its own "
                                      "gene id")
        # isoforms are realigned to and counted by name, a batch of genes at a time
        if (bed.chrom, bed.name) in names:
            raise FlairInputDataError(f"isoform name '{bed.name}' is used more than once on {bed.chrom}")
        names[(bed.chrom, bed.name)] = bed.gene_id
        gene_data[bed.gene_id].add_isoform_bed(bed)
    other_names = {name for _, name in names}
    for name, pieces in fusions.items():
        pieces.sort(key=lambda bed: bed.pos_in_fusion)
        if [bed.pos_in_fusion for bed in pieces] != list(range(len(pieces))) or len(pieces) < 2:
            raise FlairInputDataError(f"fusion isoform '{name}' needs pieces with pos_in_fusion 0, 1, ...; "
                                      f"it has {[bed.pos_in_fusion for bed in pieces]}")
        if name in other_names:
            raise FlairInputDataError(f"isoform name '{name}' is both a fusion and another isoform")
    return gene_data, fusions

def normalize_isoform_ends(gene_data, genome_file):
    """With --normalize_ends, give each spliced isoform the ends its gene's isoforms
    with terminal exons like its own share (GeneTerminalExonEnds), before batching,
    so the batches' loci and the direct assignments see the isoforms as realigned to"""
    with pysam.FastaFile(genome_file) as genome:
        for gene_info in gene_data.values():
            for bed in gene_info.isoform_beds:
                if len(bed.blocks) > 1:
                    bed.chromStart, bed.chromEnd = gene_info.terminal_exon_ends.normalized_ends(
                        bed.blocks[0].end, bed.blocks[-1].start, genome.get_reference_length(bed.chrom))
                    bed.blocks[0] = BedBlock(bed.chromStart, bed.blocks[0].end)
                    bed.blocks[-1] = BedBlock(bed.blocks[-1].start, bed.chromEnd)
            gene_info.left_bound = min(bed.chromStart for bed in gene_info.isoform_beds)
            gene_info.right_bound = max(bed.chromEnd for bed in gene_info.isoform_beds)

def make_batches(gene_data, num_batches):
    """The genes in batches realigned together, as dicts of chrom, start and end of
    their loci, the loci as (start, end), the end of the previous batch's last locus
    on the chromosome (prev_end, None for the first), gene ids, and number of
    isoforms.  Genes whose spans overlap, a locus, are always in one batch, and
    consecutive loci of a chromosome are grouped up to an even share of the isoforms
    over num_batches.  A batch takes the reads overlapping its loci (batch_reads),
    as quantify always fetched a gene's reads, but not those starting before
    prev_end, which overlap the previous batch's last locus and are its reads, so each
    read is considered by one batch."""
    by_chrom = {}
    for gene_info in gene_data.values():
        by_chrom.setdefault(gene_info.chrom, []).append(gene_info)
    per_batch = max(1, -(-sum(len(g.isoform_beds) for g in gene_data.values()) // max(1, num_batches)))
    batches = []
    for chrom, genes in by_chrom.items():
        genes.sort(key=lambda g: (g.left_bound, g.right_bound))
        loci = []   # [start, end, gene ids, isoforms]
        for gene_info in genes:
            if loci and gene_info.left_bound < loci[-1][1]:
                loci[-1][1] = max(loci[-1][1], gene_info.right_bound)
                loci[-1][2].append(gene_info.gene_id)
                loci[-1][3] += len(gene_info.isoform_beds)
            else:
                loci.append([gene_info.left_bound, gene_info.right_bound, [gene_info.gene_id], len(gene_info.isoform_beds)])
        chrom_batches = []
        for locus in loci:
            if chrom_batches and chrom_batches[-1]['isoforms'] < per_batch:
                batch = chrom_batches[-1]
                batch['end'] = locus[1]
                batch['loci'].append((locus[0], locus[1]))
                batch['gene_ids'].extend(locus[2])
                batch['isoforms'] += locus[3]
            else:
                chrom_batches.append({'chrom': chrom, 'start': locus[0], 'end': locus[1], 'loci': [(locus[0], locus[1])],
                                      'prev_end': chrom_batches[-1]['end'] if chrom_batches else None,
                                      'gene_ids': list(locus[2]), 'isoforms': locus[3]})
        batches.extend(chrom_batches)
    return batches

def _overlaps_loci(loci, start, end):
    "does start-end overlap any of loci, sorted (start, end) not overlapping each other"
    i = bisect_right(loci, (end, -1)) - 1
    return i >= 0 and loci[i][1] > start

def direct_assignment_chains(isoform_beds):
    """The isoforms reads can be assigned to without the realignment, by chromosome and
    junction chain: the spliced isoforms, but not a chain more than one isoform has,
    where the realignment chooses between them by their ends.  Fusions' pieces
    aren't genes' isoforms (load_isoform_data)."""
    by_chain = {}
    for bed in isoform_beds:
        if len(bed.blocks) > 1:
            by_chain.setdefault((bed.chrom, tuple((gap.start, gap.end) for gap in bed.getGaps())), []).append(bed)
    return {chain: beds[0] for chain, beds in by_chain.items() if len(beds) == 1}

def direct_isoform(read, chains, genome, trust_ends):
    """The isoform a read is assigned to without the realignment, or None: its
    introns are the chain of one isoform (direct_assignment_chains), aligned cleanly
    at their splice sites, as flair transcriptome's direct assignment requires of
    reads.  Unlike there, the read needn't cover the isoform's unique sequence
    past the exons of an isoform whose chain contains its own: a read is
    full-length for the one isoform with exactly its junctions, so that is its
    isoform, as in the realignment (count_batch).  With trust_ends, its
    ends must also be near the isoform's, as count_sam_transcripts requires.  A
    read with a supplementary alignment, perhaps a fusion's, is realigned."""
    if read.has_tag('SA'):
        return None
    introns = read_introns(read)
    bed = chains.get((read.reference_name, introns))
    if bed is None or not splice_sites_cleanly_aligned(read, introns, genome):
        return None
    if trust_ends and (abs(read.reference_start - bed.chromStart) > TRUST_ENDS_WINDOW
                       or abs(read.reference_end - bed.chromEnd) > TRUST_ENDS_WINDOW):
        return None
    return bed

def is_realigned_to(bed):
    """Is an isoform in the realignment's reference: spliced isoforms, but not
    single-exon isoforms, whose reads are the unspliced reads (single_exon_isoform),
    as in flair transcriptome, so a spliced read isn't counted for one with its
    intron as a deletion"""
    return len(bed.blocks) > 1

def single_exon_index(gene_data):
    """The single-exon isoforms of every gene, by chromosome, as (start, end, name,
    gene id), that unspliced reads are assigned to (single_exon_isoform)"""
    index = {}
    for gene_info in gene_data.values():
        for bed in gene_info.isoform_beds:
            if not is_realigned_to(bed):
                index.setdefault(bed.chrom, IntervalIndex()).add(
                    bed.chromStart, bed.chromEnd, (bed.chromStart, bed.chromEnd, bed.name, gene_info.gene_id))
    for chrom_index in index.values():
        chrom_index.build()
    return index

def single_exon_isoform(read, se_index):
    """The single-exon isoform an unspliced read is assigned to, as (start, end, name,
    gene id), or None, without the realignment, as flair transcriptome assigns its
    unspliced reads: of the single-exon isoforms the read covers more than half of,
    the one whose ends are nearest the read's, then the leftmost.  Strand isn't
    known.  The isoforms of every gene are candidates, whatever batch the read is in."""
    chrom_index = se_index.get(read.reference_name)
    if chrom_index is None:
        return None
    start, end = read.reference_start, read.reference_end
    covered = [iso for iso in chrom_index.overlap(start, end)
               if min(iso[1], end) - max(iso[0], start) > (iso[1] - iso[0]) / 2]
    return min(covered, key=lambda iso: (abs(start - iso[0]) + abs(end - iso[1]),) + iso, default=None)

def _genomic_clipping(read):
    "the read's soft and hard clipping at the two ends of its genome alignment"
    clipping_types = (pysam.CIGAR_OPS.CSOFT_CLIP, pysam.CIGAR_OPS.CHARD_CLIP)
    cigar = read.cigartuples
    return (cigar[0][1] if cigar[0][0] in clipping_types else 0,
            cigar[-1][1] if cigar[-1][0] in clipping_types else 0)

def batch_reads(bam_file, batch, chains, genome, trust_ends, se_index, fasta, clipping, fusion_reads=frozenset()):
    """Assign the batch's primary reads (make_batches), returning isoform name -> read
    names for those assigned without the realignment, and the number of reads
    written to fasta, with their genome alignment's clipping to clipping, to be
    realigned.  Reads matching a spliced isoform's chain go to it (direct_isoform),
    and unspliced reads, never realigned, to a single-exon isoform
    (single_exon_isoform), in this batch or not.  Reads counted for a fusion
    (count_fusions), fusion_reads, are skipped."""
    direct, num_realign = {}, 0
    loci, prev_end = batch['loci'], batch['prev_end']
    with open(fasta, 'w') as fa_fh, open(clipping, 'w') as clip_fh:
        for read in bam_file.fetch(batch['chrom'], batch['start'], batch['end']):
            if ((prev_end is not None and read.reference_start < prev_end)
                    or not _overlaps_loci(loci, read.reference_start, read.reference_end)
                    or read.is_secondary or read.is_supplementary or read.query_sequence is None
                    or read.query_name in fusion_reads):
                continue
            if not read_introns(read):
                iso = single_exon_isoform(read, se_index)
                if iso is not None:
                    direct.setdefault(iso[2], []).append(read.query_name)
                continue
            bed = direct_isoform(read, chains, genome, trust_ends)
            if bed is not None:
                direct.setdefault(bed.name, []).append(read.query_name)
            else:
                fa_fh.write('>' + read.query_name + '\n' + read.get_forward_sequence() + '\n')
                clip_fh.write(read.query_name + '\t' + '\t'.join(str(c) for c in _genomic_clipping(read)) + '\n')
                num_realign += 1
    return direct, num_realign

def read_count_sam_results(counts_file, read_map_file, assigned):
    """Add count_sam_transcripts' assignments, from its read map, or its counts without
    one, to assigned, isoform name -> read names (or counts)"""
    if read_map_file is not None:
        for line in open(read_map_file):
            name, reads = line.rstrip('\n').split('\t')
            assigned.setdefault(name, []).extend(reads.split(',') if reads else [])
    else:
        for line in open(counts_file):
            name, count = line.rstrip('\n').split('\t')
            assigned[name] = assigned.get(name, 0) + int(count)

def write_batch_reference(batch_beds, genome, temp_prefix):
    "the isoforms realigned to (is_realigned_to), as bed and fasta, and whether there are any"
    any_reference = False
    with open(temp_prefix + 'isoforms.bed', 'w') as fh_bed, open(temp_prefix + 'isoforms.fa', 'w') as fh_fa:
        for bed in batch_beds:
            if not is_realigned_to(bed):
                continue
            any_reference = True
            bed.write(fh_bed)
            seq = ''.join(genome.fetch(bed.chrom, block.start, block.end) for block in bed.blocks)
            if bed.strand == '-':
                seq = get_reverse_complement(seq)
            fh_fa.write('>' + bed.name + '\n' + seq + '\n')
    return any_reference

def count_batch(batch, batch_beds, sample_data, temp_dir, genome_file, trust_ends, generate_map, keep_intermediate):
    """Count the reads of one batch of genes in every sample.  Returns a list, per
    sample, of isoform name -> read names with generate_map, otherwise -> counts.
    Only minimap2's input and count_sam_transcripts' output are written, in the
    batch's own directory, removed after unless keep_intermediate."""
    temp_prefix = temp_dir + f"{batch['chrom']}_{batch['start']}/"
    os.makedirs(temp_prefix)
    results = []
    try:
        with pysam.FastaFile(genome_file) as genome:
            has_reference = write_batch_reference(batch_beds, genome, temp_prefix)
            for i, (sample, group, sample_batch, bamfile) in enumerate(sample_data):
                prefix = temp_prefix + sample
                with pysam.AlignmentFile(bamfile, 'rb') as bam_file:
                    direct, num_realign = batch_reads(bam_file, batch, _chains, genome, trust_ends, _single_exon_index,
                                                      prefix + '.reads.fasta', prefix + '.reads.genomicclipping.txt',
                                                      _fusion_reads[i] if _fusion_reads else frozenset())
                assigned = {name: (reads if generate_map else len(reads)) for name, reads in direct.items()}
                logging.debug(f"{batch['chrom']}:{batch['start']} {sample}: "
                              f"{sum(len(r) for r in direct.values())} reads assigned without realignment, "
                              f"{num_realign} realigned")
                if num_realign > 0 and has_reference:
                    # a read counts for an isoform it is full-length for, covering all its junctions
                    # and into its first and last exons.  Without unique sequence bounds, a read
                    # whose junctions are a contiguous run of a longer isoform's is that shorter
                    # isoform's, not given to neither: it is not full-length for the longer one.
                    # The batch's isoforms of all its genes are realigned to together, so a read
                    # fitting isoforms of two overlapping genes is ranked against both, once.
                    read_map_file = prefix + '.isoform.read.map.txt' if generate_map else None
                    mm2_cmd = ['minimap2', '-a', '-N', '4', '--MD', temp_prefix + 'isoforms.fa', prefix + '.reads.fasta']
                    # minimap2 borrows threads left idle in the thread budget while it runs
                    with thread_share.program_threads() as mm2_threads:
                        run_count_sam_transcripts(
                            mm2_cmd=mm2_cmd[:1] + ['-t', str(mm2_threads)] + mm2_cmd[1:],
                            output=prefix + '.isoform.counts.txt',
                            trimmedreads=prefix + '.reads.genomicclipping.txt',
                            generate_map=read_map_file,
                            stringent=True,
                            allow_UTR_indels=True,  # is_annot,
                            transcript_fasta=temp_prefix + 'isoforms.fa',
                            check_splice=True,
                            isoforms=temp_prefix + 'isoforms.bed',
                            trust_ends=trust_ends,
                        )
                    read_count_sam_results(prefix + '.isoform.counts.txt', read_map_file, assigned)
                results.append(assigned)
    finally:
        if not keep_intermediate:
            shutil.rmtree(temp_prefix, ignore_errors=True)
    return results

def fused_transcript(pieces, genome):
    """A fusion isoform's transcript, its pieces' exons in pos_in_fusion order, each
    piece in its own strand's orientation, as its sequence and exon sizes"""
    seqs, sizes = [], []
    for piece in pieces:
        for block in (piece.blocks if piece.strand == '+' else piece.blocks[::-1]):
            seq = genome.fetch(piece.chrom, block.start, block.end)
            seqs.append(seq if piece.strand == '+' else get_reverse_complement(seq))
            sizes.append(block.end - block.start)
    return ''.join(seqs), sizes

def write_fusion_reference(fusions, genome, temp_prefix):
    """The fusion isoforms' fused transcripts as fasta, and as a bed of plus-strand
    records, each on a sequence of its own named for it, exon after exon.  A fusion's
    junctions are then the transcript's junctions to count_sam_transcripts, so a read
    counts for a fusion only full-length across all of them, with each splice site,
    those of the fusion junctions too, checked as for any isoform."""
    with open(temp_prefix + 'fusions.bed', 'w') as fh_bed, open(temp_prefix + 'fusions.fa', 'w') as fh_fa:
        for name, pieces in fusions.items():
            seq, sizes = fused_transcript(pieces, genome)
            fh_fa.write('>' + name + '\n' + seq + '\n')
            starts = [sum(sizes[:i]) for i in range(len(sizes))]
            fh_bed.write('\t'.join([name, '0', str(len(seq)), name, '0', '+', '0', '0', '0', str(len(sizes)),
                                    ','.join(map(str, sizes)) + ',', ','.join(map(str, starts)) + ',']) + '\n')

def fusion_candidate_reads(bam_file, fusions, fasta, clipping):
    """Write the primary alignments overlapping a fusion's piece of the reads flair
    fusion takes as possibly chimeric (needs_transcriptome): those with a supplementary
    alignment (an SA tag), or FUSION_MIN_CLIP or more clipped at an end, as a short
    partner piece the genome aligner didn't place, once each, to fasta, and their
    clipping as none: a read of a fused transcript aligns to it end to end.  Returns
    the number written."""
    regions = {}
    for pieces in fusions.values():
        for piece in pieces:
            regions.setdefault(piece.chrom, []).append((piece.chromStart, piece.chromEnd))
    seen = set()
    with open(fasta, 'w') as fa_fh, open(clipping, 'w') as clip_fh:
        for chrom, spans in regions.items():
            if chrom not in bam_file.references:
                continue
            spans.sort()
            merged = [list(spans[0])]
            for start, end in spans[1:]:
                if start <= merged[-1][1]:
                    merged[-1][1] = max(merged[-1][1], end)
                else:
                    merged.append([start, end])
            for start, end in merged:
                for read in bam_file.fetch(chrom, start, end):
                    if (read.is_secondary or read.is_supplementary or read.query_sequence is None
                            or read.query_name in seen or not needs_transcriptome(read, FUSION_MIN_CLIP)):
                        continue
                    seen.add(read.query_name)
                    fa_fh.write('>' + read.query_name + '\n' + read.get_forward_sequence() + '\n')
                    clip_fh.write(read.query_name + '\t0\t0\n')
    return len(seen)

def count_fusions(fusions, sample_data, temp_dir, genome_file, threads, trust_ends, generate_map, keep_intermediate):
    """Count the fusion isoforms' reads in every sample, before the genes' (count_batch),
    by realigning the chimeric reads (fusion_candidate_reads) to the fused transcripts
    (write_fusion_reference) with count_sam_transcripts, as the genes' reads are.
    Returns a list, per sample, of isoform name -> read names with generate_map,
    otherwise -> counts, and of the sets of the reads counted, which the genes' batches
    skip, so a fusion's read isn't counted for its partners' isoforms too."""
    results, counted = [], []
    if not fusions:
        return [{} for _ in sample_data], [frozenset() for _ in sample_data]
    temp_prefix = temp_dir + 'fusions/'
    os.makedirs(temp_prefix)
    try:
        with pysam.FastaFile(genome_file) as genome:
            write_fusion_reference(fusions, genome, temp_prefix)
        for sample, group, sample_batch, bamfile in sample_data:
            prefix = temp_prefix + sample
            assigned = {}
            with pysam.AlignmentFile(bamfile, 'rb') as bam_file:
                num_reads = fusion_candidate_reads(bam_file, fusions, prefix + '.reads.fasta',
                                                   prefix + '.reads.genomicclipping.txt')
            logging.info(f'{sample}: {num_reads} chimeric reads realigned to {len(fusions)} fusion isoforms')
            if num_reads > 0:
                run_count_sam_transcripts(
                    mm2_cmd=['minimap2', '-a', '-N', '4', '--MD', '-t', str(threads), temp_prefix + 'fusions.fa',
                             prefix + '.reads.fasta'],
                    output=prefix + '.isoform.counts.txt',
                    trimmedreads=prefix + '.reads.genomicclipping.txt',
                    generate_map=prefix + '.isoform.read.map.txt',
                    stringent=True,
                    allow_UTR_indels=True,
                    transcript_fasta=temp_prefix + 'fusions.fa',
                    check_splice=True,
                    isoforms=temp_prefix + 'fusions.bed',
                    trust_ends=trust_ends,
                )
                read_count_sam_results(prefix + '.isoform.counts.txt', prefix + '.isoform.read.map.txt', assigned)
            counted.append(frozenset(read for reads in assigned.values() for read in reads))
            results.append(assigned if generate_map else {name: len(reads) for name, reads in assigned.items()})
    finally:
        if not keep_intermediate:
            shutil.rmtree(temp_prefix, ignore_errors=True)
    return results, counted

def _count_batch_task(packed):
    batch = packed[0]
    thread_share.task_start(batch['isoforms'])
    try:
        return count_batch(*packed)
    finally:
        thread_share.task_done()

def write_outputs(results_iter, sample_data, gene_data, fusions, output, with_gene, generate_map):
    """Write the counts, and with generate_map each sample's read map, from the fusions'
    and batches' results as they arrive, in genomic order.  Genes' isoforms are in the
    bed's order in the counts, then a line per fusion isoform, and read map lines as
    they were assigned."""
    counts = {}
    map_fhs = [open(f'{output}.{sample}.read.map.txt', 'w') for sample, *_ in sample_data] if generate_map else []
    try:
        for batch_results in results_iter:
            for i, assigned in enumerate(batch_results):
                for name, reads in assigned.items():
                    n = len(reads) if generate_map else reads
                    counts.setdefault(name, [0] * len(sample_data))[i] += n
                    if generate_map:
                        map_fhs[i].write(name + '\t' + ','.join(reads) + '\n')
    finally:
        for fh in map_fhs:
            fh.close()
    with open(output + '.counts.tsv', 'w') as fh:
        # just the sample id; the condition and batch of each column are in the
        # sample info file written beside this one
        fh.write('\t'.join(['ids'] + [x[0] for x in sample_data]) + '\n')
        isoform_ids = [(bed.name, gene_id) for gene_id, gene_info in gene_data.items() for bed in gene_info.isoform_beds]
        isoform_ids += [(name, pieces[0].gene_id) for name, pieces in fusions.items()]
        for name, gene_id in isoform_ids:
            iso_id = name + '_' + gene_id if with_gene else name
            fh.write('\t'.join([iso_id] + [str(x) for x in counts.get(name, [0] * len(sample_data))]) + '\n')

def quantify(*, manifest, genome, isoform_bed, output, threads, tpm,
             trust_ends, generate_map, with_gene, normalize_ends, temp_dir=None, keep_intermediate=False):
    logging.info('loading isos')
    sample_data = load_manifest(manifest)
    gene_data, fusions = load_isoform_data(isoform_bed)
    if normalize_ends:
        normalize_isoform_ends(gene_data, genome)
    # built before the batch tasks fork, which share them
    global _chains, _single_exon_index, _fusion_reads
    all_beds = [bed for gene_info in gene_data.values() for bed in gene_info.isoform_beds]
    _chains = direct_assignment_chains(all_beds)
    _single_exon_index = single_exon_index(gene_data)
    batches = make_batches(gene_data, threads * BATCHES_PER_THREAD)

    # a directory of this run's own, so no earlier run's files are read
    run_temp_dir = make_run_temp_dir(output, temp_dir)
    logging.info(f'Re-aligning reads to transcriptome in {len(batches)} batches of genes'
                 f'{f" and {len(fusions)} fusion isoforms" if fusions else ""}. '
                 f'Writing temporary files into {run_temp_dir}')
    try:
        fusion_results, _fusion_reads = count_fusions(fusions, sample_data, run_temp_dir, genome, threads, trust_ends,
                                                      generate_map, keep_intermediate)
        packed = [(batch, [bed for gene_id in batch['gene_ids'] for bed in gene_data[gene_id].isoform_beds],
                   sample_data, run_temp_dir, genome, trust_ends, generate_map, keep_intermediate)
                  for batch in batches]
        # minimap2 borrows threads left idle as the pool drains, more for batches of more isoforms
        thread_share.init(threads, len(packed), max((batch['isoforms'] for batch in batches), default=1))
        if threads == 1:
            write_outputs(itertools.chain([fusion_results], map(_count_batch_task, packed)),
                          sample_data, gene_data, fusions, output, with_gene, generate_map)
        else:
            mp.set_start_method('fork', force=True)
            with mp.Pool(threads) as pool:
                write_outputs(itertools.chain([fusion_results], pool.imap(_count_batch_task, packed)),
                              sample_data, gene_data, fusions, output, with_gene, generate_map)
        write_sample_info(sample_info_path(output + '.counts.tsv'),
                          [SampleInfo(sample, condition, batch)
                           for sample, condition, batch, _ in sample_data])
        if tpm:
            counts_to_tpm(open(output + '.counts.tsv'), output + '.tpm.tsv')
    finally:
        # also on failure, so runs don't leave their files in $TMPDIR
        if keep_intermediate:
            logging.info(f'intermediate files kept in {run_temp_dir}')
        else:
            shutil.rmtree(run_temp_dir, ignore_errors=True)
