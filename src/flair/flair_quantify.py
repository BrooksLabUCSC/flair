#! /usr/bin/env python3

import os
import pysam
from shutil import rmtree
import logging
from flair import FlairInputDataError
from flair.counts_matrix import SampleInfo, sample_info_path, write_sample_info
from flair.io_utils import make_temp_dir
from flair import thread_share
from flair.pycbio.hgdata.bed import BedReader, BedBlock
from flair.flair_bed import FlairBed
from flair.terminal_exon_ends import GeneTerminalExonEnds
from flair.read_processing import generate_genomic_alignment_read_to_clipping_file
from flair.count_sam_transcripts import run_count_sam_transcripts, TRUST_ENDS_WINDOW
from flair.annotation_precheck import read_introns, splice_sites_cleanly_aligned
from flair.interval_index import IntervalIndex
from flair.counts_to_tpm import counts_to_tpm
from flair.isoform_data import get_reverse_complement
import multiprocessing as mp


os.environ['OPENBLAS_NUM_THREADS'] = '1'

# the single-exon isoforms of every gene (single_exon_index), set by quantify
_single_exon_index = {}

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
    # parser.add_argument('--stringent', action='store_true',
    #                     help="""Supporting reads must cover 80 percent of their isoform and extend at least 25 nt into the
    #                     first and last exons. If those exons are themselves shorter than 25 nt, the requirement becomes
    #                     'must start within 4 nt from the start" or "must end within 4 nt from the end" """)
    # parser.add_argument('--check_splice', action='store_true',
    #                     help="""enforce coverage of 4 out of 6 bp around each splice site and no
    #                     insertions greater than 3 bp at the splice site""")
    parser.set_defaults(entry=quantify_cmd)

def quantify_cmd(args):
    check_input_files(manifest=args.manifest, genome=args.genome, isoform_bed=args.isoform_bed)
    quantify(manifest=args.manifest, genome=args.genome, isoform_bed=args.isoform_bed,
             output=args.output, threads=args.threads,
             tpm=args.tpm, trust_ends=args.trust_ends, generate_map=args.generate_map,
             with_gene=args.with_gene, normalize_ends=args.normalize_ends)

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
    gene_data = {}
    for bed in BedReader(isoform_bed, bedClass=FlairBed):
        if bed.gene_id not in gene_data:
            gene_data[bed.gene_id] = GeneData(bed.gene_id, bed.chrom, bed.strand)
        elif bed.strand != gene_data[bed.gene_id].strand:
            raise FlairInputDataError(f"gene id '{bed.gene_id}' has isoforms on both strands, "
                                      f"'{gene_data[bed.gene_id].strand}' and '{bed.strand}'; give each strand its own "
                                      "gene id")
        gene_data[bed.gene_id].add_isoform_bed(bed)
    return gene_data

def write_combined_counts(sample_data, gene_data, temp_dir, output, with_gene):
    with open(output + '.counts.tsv', 'w') as fh:
        # just the sample id; the condition and batch of each column are in the
        # sample info file written beside this one
        fh.write('\t'.join(['ids'] + [x[0] for x in sample_data]) + '\n')
        for gene_id in gene_data:
            gene_info = gene_data[gene_id]
            iso_to_counts = {x.name: [0] * len(sample_data) for x in gene_info.isoform_beds}
            for i, (sample, group, batch, bamfile) in enumerate(sample_data):
                for line in open(temp_dir + gene_id + '/' + sample + '.isoform.counts.txt'):
                    line = line.rstrip().split('\t')
                    iso_to_counts[line[0]][i] = int(line[1])
            for iso in iso_to_counts:
                iso_id = iso + '_' + gene_id if with_gene else iso
                fh.write('\t'.join([iso_id] + [str(x) for x in iso_to_counts[iso]]) + '\n')

def write_map_out(sample_data, gene_data, temp_dir, output, generate_map):
    if generate_map:
        for sample, group, batch, bamfile in sample_data:
            with open(f'{output}.{sample}.read.map.txt', 'w') as fh:
                for gene_id in gene_data:
                    for line in open(temp_dir + gene_id + '/' + sample + '.isoform.read.map.txt'):
                        fh.write(line)

def direct_assignment_chains(isoform_beds):
    """The isoforms reads can be assigned to without the realignment, by junction
    chain: the gene's spliced isoforms, but not a fusion's pieces, whose reads need
    their partner too, nor a chain more than one isoform has, where the realignment
    chooses between them by their ends"""
    by_chain = {}
    for bed in isoform_beds:
        if len(bed.blocks) > 1 and bed.pos_in_fusion is None and bed.transcript_class != 'fusion':
            by_chain.setdefault(tuple((gap.start, gap.end) for gap in bed.getGaps()), []).append(bed)
    return {chain: beds[0] for chain, beds in by_chain.items() if len(beds) == 1}

def direct_isoform(read, chains, genome, trust_ends):
    """The isoform a read is assigned to without the realignment, or None: its
    introns are the chain of one isoform (direct_assignment_chains), aligned cleanly
    at their splice sites, as flair transcriptome's direct assignment requires of
    reads.  Unlike there, the read needn't cover the isoform's unique sequence
    past the exons of an isoform whose chain contains its own: a read is
    full-length for the one isoform with exactly its junctions, so that is its
    isoform, as in the realignment (get_counts_for_sample).  With trust_ends, its
    ends must also be near the isoform's, as count_sam_transcripts requires.  A
    read with a supplementary alignment, perhaps a fusion's, is realigned."""
    if read.has_tag('SA'):
        return None
    introns = read_introns(read)
    bed = chains.get(introns)
    if bed is None or not splice_sites_cleanly_aligned(read, introns, genome):
        return None
    if trust_ends and (abs(read.reference_start - bed.chromStart) > TRUST_ENDS_WINDOW
                       or abs(read.reference_end - bed.chromEnd) > TRUST_ENDS_WINDOW):
        return None
    return bed

def is_realigned_to(bed):
    """Is an isoform in the realignment's reference: spliced isoforms and fusions'
    pieces, but not other single-exon isoforms, whose reads are the unspliced reads
    (single_exon_isoform), as in flair transcriptome, so a spliced read isn't
    counted for one with its intron as a deletion"""
    return len(bed.blocks) > 1 or bed.pos_in_fusion is not None or bed.transcript_class == 'fusion'

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
    known.  The isoforms of every gene are candidates, so each gene the read is
    fetched for finds the same one, and the read is counted once."""
    chrom_index = se_index.get(read.reference_name)
    if chrom_index is None:
        return None
    start, end = read.reference_start, read.reference_end
    covered = [iso for iso in chrom_index.overlap(start, end)
               if min(iso[1], end) - max(iso[0], start) > (iso[1] - iso[0]) / 2]
    return min(covered, key=lambda iso: (abs(start - iso[0]) + abs(end - iso[1]),) + iso, default=None)

def write_reads_to_realign(bam_file, gene_info, chains, genome, trust_ends, fasta, se_index):
    """Write the gene's primary reads to fasta to be realigned, as samtools fasta would,
    except the reads assigned without the realignment, which are returned as
    isoform name -> read names, with the number of reads written: reads matching a
    spliced isoform (direct_isoform), and unspliced reads, which are never realigned,
    assigned to a single-exon isoform (single_exon_isoform), if it is this gene's"""
    direct, num_realign = {}, 0
    with open(fasta, 'w') as fh:
        # samtools' region, with the 1-based start
        for read in bam_file.fetch(region=f'{gene_info.chrom}:{gene_info.left_bound}-{gene_info.right_bound}'):
            if read.is_secondary or read.is_supplementary or read.query_sequence is None:
                continue
            if not read_introns(read):
                iso = single_exon_isoform(read, se_index)
                if iso is not None and iso[3] == gene_info.gene_id:
                    direct.setdefault(iso[2], []).append(read.query_name)
                continue
            bed = direct_isoform(read, chains, genome, trust_ends)
            if bed is not None:
                direct.setdefault(bed.name, []).append(read.query_name)
            else:
                fh.write('>' + read.query_name + '\n' + read.get_forward_sequence() + '\n')
                num_realign += 1
    return direct, num_realign

def add_direct_assignments(direct, counts_file, read_map_file):
    """Add the reads assigned without the realignment, isoform name -> read names, to
    count_sam_transcripts' counts and read map"""
    counts = {}
    for line in open(counts_file):
        name, count = line.rstrip('\n').split('\t')
        counts[name] = int(count)
    for name, reads in direct.items():
        counts[name] = counts.get(name, 0) + len(reads)
    with open(counts_file, 'w') as fh:
        for name, count in counts.items():
            fh.write(name + '\t' + str(count) + '\n')
    if read_map_file is not None:
        read_map = {}
        for line in open(read_map_file):
            name, reads = line.rstrip('\n').split('\t')
            read_map[name] = reads.split(',') if reads else []
        for name, reads in direct.items():
            read_map.setdefault(name, []).extend(reads)
        with open(read_map_file, 'w') as fh:
            for name, reads in read_map.items():
                fh.write(name + '\t' + ','.join(reads) + '\n')

def get_counts_for_sample(sample, bamfile, temp_prefix, gene_info, generate_map, trust_ends, chains, genome,
                          se_index, has_realign_reference):
    temp_prefix_sample = temp_prefix + sample
    # with, not a bare open: this runs once per gene per sample, so the descriptors
    # grew without bound in a worker handling many genes
    with pysam.AlignmentFile(bamfile, 'rb') as bam_file:
        direct, num_realign = write_reads_to_realign(bam_file, gene_info, chains, genome, trust_ends,
                                                     temp_prefix_sample + '.reads.fasta', se_index)
        generate_genomic_alignment_read_to_clipping_file(temp_prefix_sample, bam_file, gene_info.chrom, gene_info.left_bound, gene_info.right_bound)
    logging.debug(f'{gene_info.gene_id} {sample}: {sum(len(r) for r in direct.values())} reads assigned without '
                  f'realignment, {num_realign} realigned')
    read_map_file = temp_prefix_sample + '.isoform.read.map.txt' if generate_map else None
    if num_realign == 0 or not has_realign_reference:
        for path in (temp_prefix_sample + '.isoform.counts.txt', read_map_file):
            if path is not None:
                open(path, 'w').close()
    else:
        mm2_cmd = ['minimap2', '-a', '-N', '4', '--MD', temp_prefix + 'isoforms.fa', temp_prefix_sample + '.reads.fasta']
        # a read counts for an isoform it is full-length for, covering all its junctions
        # and into its first and last exons.  Without unique sequence bounds, a read
        # whose junctions are a contiguous run of a longer isoform's is that shorter
        # isoform's, not given to neither: it is not full-length for the longer one.
        # minimap2 borrows threads left idle in the thread budget while it runs
        with thread_share.program_threads() as mm2_threads:
            run_count_sam_transcripts(
                mm2_cmd=mm2_cmd[:1] + ['-t', str(mm2_threads)] + mm2_cmd[1:],
                output=temp_prefix_sample + '.isoform.counts.txt',
                trimmedreads=temp_prefix_sample + '.reads.genomicclipping.txt',
                generate_map=read_map_file,
                stringent=True,
                allow_UTR_indels=True,  # is_annot,
                transcript_fasta=temp_prefix + 'isoforms.fa',
                check_splice=True,
                isoforms=temp_prefix + 'isoforms.bed',
                trust_ends=trust_ends,
            )
    add_direct_assignments(direct, temp_prefix_sample + '.isoform.counts.txt', read_map_file)

def get_counts_for_gene(input):
    temp_dir, gene_id, gene_info, sample_data, generate_map, trust_ends, genome_file, norm_ends = input
    logging.debug(f'realigning reads to {gene_id}')
    temp_prefix = temp_dir + gene_id + '/'
    os.makedirs(temp_prefix)
    with pysam.FastaFile(genome_file) as genome:
        has_realign_reference = False
        with open(temp_prefix + 'isoforms.bed', 'w') as fh_bed, open(temp_prefix + 'isoforms.fa', 'w') as fh_fa:
            for bed in gene_info.isoform_beds:
                if not is_realigned_to(bed):
                    continue
                has_realign_reference = True
                if norm_ends and len(bed.blocks) > 1:
                    bed.chromStart, bed.chromEnd = gene_info.terminal_exon_ends.normalized_ends(bed.blocks[0].end, bed.blocks[-1].start,
                                                                                                genome.get_reference_length(bed.chrom))
                    bed.blocks[0] = BedBlock(bed.chromStart, bed.blocks[0].end)
                    bed.blocks[-1] = BedBlock(bed.blocks[-1].start, bed.chromEnd)
                bed.write(fh_bed)
                seq = []
                for block in bed.blocks:
                    seq.append(genome.fetch(bed.chrom, block.start, block.end))
                seq = ''.join(seq)
                if bed.strand == '-':
                    seq = get_reverse_complement(seq)
                fh_fa.write('>' + bed.name + '\n' + seq + '\n')

        # the isoforms with their ends as realigned to, normalized or not
        chains = direct_assignment_chains(gene_info.isoform_beds)
        for sample, group, batch, bamfile in sample_data:
            get_counts_for_sample(sample, bamfile, temp_prefix, gene_info, generate_map, trust_ends, chains, genome,
                                  _single_exon_index, has_realign_reference)

def _count_gene_task(input):
    thread_share.task_start()
    try:
        get_counts_for_gene(input)
    finally:
        thread_share.task_done()

def quantify(*, manifest, genome, isoform_bed, output, threads, tpm,
             trust_ends, generate_map, with_gene, normalize_ends):
    logging.info('loading isos')
    temp_dir = make_temp_dir(output)
    sample_data = load_manifest(manifest)
    gene_data = load_isoform_data(isoform_bed)
    # built before the gene tasks fork, which share it
    global _single_exon_index
    _single_exon_index = single_exon_index(gene_data)

    logging.info(f'Re-aligning reads to transcriptome. Writing temporary files into {temp_dir}')

    packed = []
    for gene_id in gene_data:
        packed.append((temp_dir, gene_id, gene_data[gene_id], sample_data, generate_map, trust_ends, genome, normalize_ends))

    # minimap2 borrows threads left idle as the pool drains
    thread_share.init(threads, len(packed))
    if threads == 1:
        for p in packed:
            _count_gene_task(p)
    else:
        mp.set_start_method('fork', force=True)
        with mp.Pool(threads) as pool:
            pool.map(_count_gene_task, packed)

    logging.info('writing quantify output')
    write_combined_counts(sample_data, gene_data, temp_dir, output, with_gene)
    write_sample_info(sample_info_path(output + '.counts.tsv'),
                      [SampleInfo(sample, condition, batch)
                       for sample, condition, batch, _ in sample_data])

    write_map_out(sample_data, gene_data, temp_dir, output, generate_map)

    if tpm:
        counts_to_tpm(open(output + '.counts.tsv'), output + '.tpm.tsv')

    rmtree(temp_dir)
