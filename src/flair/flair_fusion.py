#! /usr/bin/env python3

import os
import glob
import pipettor
import pysam
import logging
from collections import namedtuple
from flair.gtf_to_bed import gtf_data_to_bed
from flair.convert_synthetic_to_genome_bed import convert_synthetic_isos, get_paralog_ref, paralog_groups_for_genes
from flair.identify_prelim_fusions import id_chimeras, scan_genome_alignment
from flair import FlairInputDataError
from flair.gtf_io import gtf_data_parser, GtfAttrsSet, TRANSCRIPT_EXON_FEATURES
from flair.annotation_data import AnnotsByChrom
from flair.make_synthetic_fusion_reference import make_synthetic_reference
from flair.read_processing import get_sequence_from_bed
from flair.pycbio.hgdata.bed import Bed, BedReader
from flair.flair_bed import FlairBed
from flair.bed_to_gtf import bed_to_gtf
from flair.isoform_data import make_big_bed
from flair.fusion_table import write_fusion_table, write_table

def add_subparser(subparsers):
    desc = ("Identify gene fusions and generate a fusion transcriptome.  Reads must be stranded, each "
            "oriented 5' to 3' as its RNA: the 5' and 3' genes of a fusion come from the order of a read's "
            "pieces, so a reversed read gives the fusion backwards, which is then rejected, and poorly "
            "stranded data will miss fusions.")
    parser = subparsers.add_parser('fusion', help="Identify gene fusions", description=desc)
    required = parser.add_argument_group('required named arguments')
    required.add_argument('-g', '--genome', type=str, required=True,
                          help='FastA of reference genome')
    required.add_argument('-f', '--gtf', type=str, required=True,
                          help='GTF annotation file, used for renaming FLAIR isoforms to annotated isoforms and adjusting TSS/TESs')
    required.add_argument('-b', '--genome_aligned_bam', type=str, required=True,
                          help='bam file of reads aligned to the genome with minimap2 -ax splice.  All reads, '
                               'not only the chimeric ones: reads are also aligned to annotated transcripts, which '
                               'finds partner pieces the genome alignment leaves clipped.  The chimeric reads taken '
                               'from it are kept as OUTPUT.genomealigned.chim.bam and '
                               'OUTPUT.transcriptomealigned.chim.bam, and reused if they exist: with a different '
                               'bam, use a new -o or remove them first')
    required.add_argument('--sample_name', required=True,
                          help='name of sample, will be added as metadata to output files')
    parser.add_argument('-o', '--output',
                        help='output file name base for FLAIR isoforms, defaults to --sample_name.  The chimeric '
                             'reads bams OUTPUT.genomealigned.chim.bam and OUTPUT.transcriptomealigned.chim.bam '
                             'from an earlier run with this -o are reused, whatever -b and --min_fragment_size are')
    parser.add_argument('-t', '--threads', type=int, default=4,
                        help='number of threads for minimap2 and for flair transcriptome on the synthetic '
                             'fusion genome (default: %(default)s)')
    parser.add_argument('--min_fragment_size', type=int, default=40,
                        help='minimum size of alignment kept, used in minimap -s. More important '
                             'when doing downstream fusion detection (default: %(default)s)')
    parser.add_argument('--min_support', type=float, default=3.0,
                        help='minimum number of supporting reads for a fusion (default: %(default)s)')
    parser.add_argument('--max_loci', type=int, default=2,
                        help='max loci detected in fusion. Set higher for detection of 3-gene+ fusions (default: %(default)s)')
    parser.add_argument('--max_dist_to_tss', type=int, default=15000,
                        help="maximum allowed distance of a fusion's 5' end to an annotated transcript start "
                             "of its 5' gene upstream of it.  With fragmented or 5' truncated reads, which start "
                             "far from any TSS, set this very high or to -1, which does not check this "
                             "(default: %(default)s)")
    parser.add_argument('--min_dist_between_breakpoints', type=int, default=100000,
                        help='minimum allowed distance between breakpoints when they are on the same '
                             'strand; removes read-through transcripts.  A fusion from a deletion between genes '
                             'closer than the maximum intron length of the genome alignment (minimap2 -G) aligns '
                             'as one spliced read and is not found: to find fusions between close genes, align '
                             'with a smaller minimap2 -G and lower this to match (default: %(default)s)')
    parser.add_argument('--max_partners', type=int, default=1000,
                        help='maximum number of fusion partners of a gene in one role (5\' gene, 3\' gene or '
                             'between them) before its fusions are rejected as promiscuous.  A gene with many '
                             'partners can be an alignment or library artifact, but in highly rearranged genomes, '
                             'such as cancer cell lines with amplified regions, genes have many real partners: in '
                             'MCF7, SULF2 is the 5\' gene of 5 fusions and the 3\' gene of ARFGEF2--SULF2.  The '
                             'high default keeps them; lower it for genomes without substantially rearranged genomes '
                             '(default: %(default)s)')
    parser.add_argument('--quality', type=int, default=10,
                        help='minimum MAPQ of every piece of a chimeric read in the genome alignment for the '
                             'read to count toward a fusion. Applies to the genome alignments only, not to the '
                             'alignment of reads to annotated transcripts, where a read of a gene with several '
                             'overlapping isoforms has a low MAPQ regardless (default: %(default)s)')
    parser.add_argument('--frac_support', type=float, default=0.05,
                        help='minimum fraction of a fusion\'s reads for a fusion isoform to be reported; '
                             'a fusion with one dominant isoform hides its minor ones below this. Set to 0 '
                             'for max recall (default: %(default)s)')
    parser.add_argument('--no_transcriptome_alignment', action='store_true',
                        help='find fusions from the chimeric reads of the genome alignment only, without '
                             'aligning reads to annotated transcripts.  About half the time and a third of the '
                             'memory, and -b can then hold only the chimeric reads, but loses fusion reads whose '
                             'partner piece the genome alignment leaves clipped, particularly short 5\' exons')
    parser.add_argument('--keep_intermediate', action='store_true',
                        help='keep the intermediate files: the annotated transcripts, the preliminary fusion '
                             'reads, and the synthetic fusion genome with its annotation, breakpoints, '
                             'alignments and flair transcriptome output, where the fusion isoforms can be '
                             'viewed.  The chimeric reads bams and OUTPUT.prelimfusions.bed are always kept')
    parser.set_defaults(entry=fusion_cmd)

def fusion_cmd(args):
    for what, path in (('genome fasta', args.genome), ('GTF annotation', args.gtf)):
        if not os.path.exists(path):
            raise FlairInputDataError(f'{what} file does not exist: {path}')
    detectfusions(genome=args.genome, gtf=args.gtf,
                  genome_aligned_bam=args.genome_aligned_bam, sample_name=args.sample_name,
                  output=args.output if args.output is not None else args.sample_name,
                  threads=args.threads, min_fragment_size=args.min_fragment_size,
                  min_support=args.min_support, max_loci=args.max_loci,
                  # -1 turns the TSS distance check off
                  max_dist_to_tss=None if args.max_dist_to_tss == -1 else args.max_dist_to_tss,
                  min_dist_between_breakpoints=args.min_dist_between_breakpoints,
                  frac_support=args.frac_support, quality=args.quality, max_partners=args.max_partners,
                  transcriptome_alignment=not args.no_transcriptome_alignment, keep_intermediate=args.keep_intermediate)

def report_nofusions(outputprefix):
    logging.info('no fusions detected. Exiting')
    for file in [outputprefix + '.fusions.isoforms.bed', outputprefix + '.fusions.isoforms.fa',
                 outputprefix + '.fusions.isoform.counts.txt', outputprefix + '.fusions.isoform.read.map.txt']:
        f = open(file, 'w')
        f.close()
    with open(outputprefix + '.fusions.aaseq.tsv', 'w') as fh:
        fh.write('aaseq_id\taaseq\n')
    write_table([], outputprefix + '.fusions.tsv')

def align_to_synth_genome(genome, reads, output, additional_options):
    mm2_cmd = ['minimap2', '-ax', 'splice'] + additional_options + [genome, reads]
    samtools_filter_cmd = ('samtools', 'view', '-F', '2048', '-hb', '-')
    samtools_sort_cmd = ('samtools', 'sort', '-o', output, '-')
    samtools_index_cmd = ('samtools', 'index', output)
    pipettor.run([mm2_cmd, samtools_filter_cmd, samtools_sort_cmd])
    pipettor.run([samtools_index_cmd])


CANONICAL_MOTIFS = {"GT/AG", "GC/AG", "AT/AC"}
# a non-canonical breakpoint junction must carry this share of the reads crossing the breakpoint
NONCANONICAL_MIN_SHARE = 0.5


def count_breakpoint_reads(bam, fusiontobp):
    """fusion contig -> the number of reads aligned across its breakpoint"""
    crossing = {f: 0 for f in fusiontobp}
    with pysam.AlignmentFile(bam, 'rb') as aligns:
        for a in aligns:
            if a.is_mapped and not a.is_secondary and not a.is_supplementary:
                breakpoint = fusiontobp.get(a.reference_name)
                if breakpoint is not None and a.reference_start < breakpoint < a.reference_end:
                    crossing[a.reference_name] += 1
    return crossing


def select_junctions(junctions, fusiontobp, crossing):
    """The rows of the junctions (intron-prospector BED6 records on the synthetic
    genome) to pass flair transcriptome: canonical ones with 2 or more reads, and for
    a fusion with no such junction across its breakpoint, the junctions across it with
    2 or more reads that carry NONCANONICAL_MIN_SHARE of the reads crossing it.  A
    breakpoint inside an exon has no splice site, so its junction needn't be
    canonical, but the reads should agree on it: in MCF7, TAOK3--RPS2's 21 reads
    crossed its breakpoint at a dozen non-canonical junctions, and 3 that happened to
    agree made a fusion isoform."""
    junctions = [bed for bed in junctions if bed.score >= 2]
    good_sj, has_canonical_bp = [], set()
    for bed in junctions:
        # for synthetic alignment, all junctions should be '+'
        if bed.name.split('_')[-1] in CANONICAL_MOTIFS and bed.strand == '+':
            good_sj.append(bed.toRow())
            breakpoint = fusiontobp.get(bed.chrom)
            if breakpoint is not None and bed.chromStart < breakpoint < bed.chromEnd:
                has_canonical_bp.add(bed.chrom)
    for bed in junctions:
        breakpoint = fusiontobp.get(bed.chrom)
        if breakpoint is not None and bed.chrom not in has_canonical_bp and bed.chromStart < breakpoint < bed.chromEnd \
                and bed.score >= NONCANONICAL_MIN_SHARE * crossing.get(bed.chrom, 0):
            good_sj.append(bed.toRow())
    return good_sj


Annotation = namedtuple('Annotation', ['genetoinfo', 'genetoname', 'chrom_to_gene_pos', 'gene_to_all_exons',
                                       'juncs_to_gene', 'intronLocs', 'intronToGenome', 'gene_to_paralogs'])


def align_reads_to_transcripts(genome, gtf_data, output, reads_fa, nreads, threads, min_fragment_size):
    "the chimeric alignments of reads_fa to the annotated transcripts, in OUTPUT.transcriptomealigned.chim.bam"
    logging.info('aligning to transcriptome and getting chimeric reads')
    transcriptchimbam = output + '.transcriptomealigned.chim.bam'
    annotated_bed = output + '.annotated_transcripts.bed'
    gtf_data_to_bed(annotated_bed, gtf_data, include_gene=True, name_sep='|')
    annotated_fa = output + '.annotated_transcripts.fa'
    get_sequence_from_bed(genome, output + '.annotated_transcripts.bed', annotated_fa)

    logging.info(f'{nreads} reads to align to transcripts')
    mm2_cmd = ('minimap2', '-a', '-s', str(min_fragment_size), '-t', str(threads), '--secondary=no',
               annotated_fa, reads_fa)

    filter_cmd = ('samtools', 'view', '-hF', '0x104', '-e', '[SA] != ""')
    pipettor.run([mm2_cmd, filter_cmd], stdout=transcriptchimbam)
    os.remove(reads_fa)
    pipettor.run([('samtools', 'sort', '-o', output + '.transcriptomealigned.chim.sorted.bam', transcriptchimbam)])
    pipettor.run([('mv', output + '.transcriptomealigned.chim.sorted.bam', transcriptchimbam)])
    pysam.index(transcriptchimbam)


def make_chimeric_bams(genome, gtf_data, genome_aligned_bam, output, threads, min_fragment_size, transcriptome_alignment):
    """The chimeric reads of the genome alignment and of the alignment to annotated
    transcripts, each unless an earlier run left them, from one pass over the genome
    alignment.  Returns the names of the reads aligned to transcripts, None if not."""
    genomechimbam = output + '.genomealigned.chim.bam'
    make_chim = not os.path.exists(genomechimbam)
    align_to_transcripts = transcriptome_alignment and not os.path.exists(output + '.transcriptomealigned.chim.bam')
    transcriptome_reads_fa = output + '.transcriptome_reads.fa'
    transcriptome_reads = None
    if make_chim or transcriptome_alignment:
        logging.info('getting chimeric reads from genome')
        transcriptome_reads = scan_genome_alignment(genome_aligned_bam, min_fragment_size,
                                                    chim_bam=genomechimbam if make_chim else None,
                                                    reads_fa=transcriptome_reads_fa if align_to_transcripts else None)
        if make_chim:
            pysam.index(genomechimbam)
    if align_to_transcripts:
        align_reads_to_transcripts(genome, gtf_data, output, transcriptome_reads_fa, len(transcriptome_reads),
                                   threads, min_fragment_size)
    return transcriptome_reads


def parse_gtf(gtf_data):
    """From the parsed annotation (GtfData): gene -> [chrom, start, end, strand,
    transcript 5' ends], gene -> transcript -> exons, gene -> name, chrom -> genes'
    (start, end, strand, gene), and chrom -> an empty junction dict.  A gene spans its
    transcripts on the chromosome of its first one"""
    genetoinfo, genetoexons, genetoname = {}, {}, {}
    chrom_to_gene_pos, juncs_to_gene = {}, {}
    for transcript in gtf_data.transcripts:
        gene_id = transcript.gene_id
        if gene_id not in genetoinfo:
            genetoinfo[gene_id] = [transcript.chrom, transcript.start, transcript.end, transcript.strand, []]
            genetoname[gene_id] = transcript.gene_name if transcript.gene_name else gene_id
            juncs_to_gene.setdefault(transcript.chrom, {})
        info = genetoinfo[gene_id]
        if transcript.chrom != info[0]:
            continue
        info[1], info[2] = min(info[1], transcript.start), max(info[2], transcript.end)
        info[4].append(transcript.start if transcript.strand == '+' else transcript.end)
        if len(transcript.exons) > 0:
            genetoexons.setdefault(gene_id, {})[transcript.transcript_id] = [(e.start, e.end) for e in transcript.exons]
    for gene_id, (chrom, start, end, strand, _) in genetoinfo.items():
        chrom_to_gene_pos.setdefault(chrom, []).append((start, end, strand, gene_id))
    for chrom in chrom_to_gene_pos:
        chrom_to_gene_pos[chrom] = sorted(chrom_to_gene_pos[chrom])
    return genetoinfo, genetoexons, genetoname, chrom_to_gene_pos, juncs_to_gene


def merge_exons(exons):
    "sorted exons merged where they overlap"
    newexons = []
    laststart, lastend = -1, -1
    for s, e in exons:
        if s > lastend:
            if laststart != -1:
                newexons.append((laststart, lastend))
            laststart, lastend = s, e
        else:
            lastend = max(lastend, e)
    if laststart != -1:
        newexons.append((laststart, lastend))
    return newexons


def gene_exons_and_junctions(genetoexons, genetoinfo, juncs_to_gene):
    """gene -> its transcripts' exons merged, and each annotated junction -> its genes,
    filled into juncs_to_gene by chrom"""
    gene_to_all_exons = {}
    for gene in genetoexons:
        allexons = set()
        chrom = genetoinfo[gene][0]
        for t in genetoexons[gene]:
            exons = sorted(list(genetoexons[gene][t]))
            allexons.update(set(exons))
            for x in range(len(exons) - 1):
                juncs_to_gene[chrom].setdefault((exons[x][1], exons[x + 1][0]), set()).add(gene)
        gene_to_all_exons[gene] = merge_exons(sorted(list(allexons)))
    return gene_to_all_exons


def intron_locs(exons, strand):
    """transcript positions of a transcript's splice sites, with its two ends, each with
    the genomic interval around it: the intron, or 500 bases past an end"""
    myexons = sorted(exons)
    mylocs = [[0, myexons[0][0] - 500, myexons[0][0]]]  # add start of transcript
    runningtot = 0
    for i in range(len(myexons) - 1):
        runningtot += myexons[i][1] - myexons[i][0]  # add size of last exon
        mylocs.append([runningtot, myexons[i][1], myexons[i + 1][0]])  # add intron
    runningtot += myexons[-1][1] - myexons[-1][0]
    mylocs.append([runningtot, myexons[-1][1], myexons[-1][1] + 500])
    if strand == '-':
        mylocs = [[runningtot - mylocs[x][0], mylocs[x][1], mylocs[x][2]] for x in range(len(mylocs))]
    return sorted([x[0] for x in mylocs]), {x[0]: (x[1], x[2]) for x in mylocs}


def read_annotation(gtf_data):
    "the annotation (GtfData) as id_chimeras uses it"
    genetoinfo, genetoexons, genetoname, chrom_to_gene_pos, juncs_to_gene = parse_gtf(gtf_data)
    logging.info('continuing to parse annot')
    gene_to_all_exons = gene_exons_and_junctions(genetoexons, genetoinfo, juncs_to_gene)
    intronLocs, intronToGenome = {}, {}
    for g in genetoexons:
        strand = genetoinfo[g][3]
        for t in genetoexons[g]:
            intronLocs[t], intronToGenome[t] = intron_locs(genetoexons[g][t], strand)
    gene_to_paralogs = paralog_groups_for_genes(
        get_paralog_ref(os.path.realpath(__file__).split('flair_fusion')[0] + 'dgd_Hsa_all_v71.tsv'), genetoinfo)
    return Annotation(genetoinfo, genetoname, chrom_to_gene_pos, gene_to_all_exons, juncs_to_gene,
                      intronLocs, intronToGenome, gene_to_paralogs)


def write_prelim_fusions(combchim, bedfile):
    "write the preliminary fusions' loci to bedfile, and return their reads"
    fusionreads = set()
    with open(bedfile, 'w') as bedout:
        for f in combchim:
            fusionreads.update(combchim[f]['reads'])
            genes = f.split('__')
            startdiststr = ','.join([str(x) for x in combchim[f]['disttostart']])
            for g in genes:
                strand = '+' if combchim[f][g][1] < combchim[f][g][2] else '-'
                positions = sorted(combchim[f][g][1:])
                Bed(combchim[f][g][0], positions[0], positions[1],
                    name=g + '__' + '--'.join(genes), score=len(combchim[f]['reads']), strand=strand,
                    extraCols=[startdiststr]).write(bedout)
    return fusionreads


def write_fusion_reads(fusionreads, bams, fasta):
    "the sequences of fusionreads, from the first of the chimeric read bams holding each, to fasta"
    seenreads = set()
    with open(fasta, 'w') as out_fa:
        for bam in bams:
            with pysam.AlignmentFile(bam, 'rb') as bamfile:
                for a in bamfile:
                    if not a.is_secondary and not a.is_supplementary and a.query_name in fusionreads \
                            and a.query_name not in seenreads:
                        out_fa.write('>' + a.query_name + '\n' + a.get_forward_sequence() + '\n')
                        seenreads.add(a.query_name)


def make_synthetic_genome(gtf, gtf_data, genome, output):
    "the synthetic fusion genome and its annotation; False if it holds no fusion"
    logging.info('generating synthetic reference')
    make_synthetic_reference(gtf_data, gtf, genome, output + '.prelimfusions.bed', output)
    if os.path.getsize(output + '-syntheticFusionGenome.fa') == 0:
        return False
    pipettor.run([['samtools', 'faidx', output + '-syntheticFusionGenome.fa']])
    return True


def is_primary(a):
    return a.is_mapped and not a.is_secondary and not a.is_supplementary


def choose_synthetic_alignments(withsplice_bam, nosplice_bam):
    """read -> its alignment to the synthetic genome: with splicing, or without when that
    matches more bases with fewer insertions"""
    rname_to_read = {}
    with pysam.AlignmentFile(withsplice_bam, 'rb') as bamfile:
        for a in bamfile:
            if is_primary(a):
                rname_to_read[a.query_name] = a
    with pysam.AlignmentFile(nosplice_bam, 'rb') as bamfile:
        for a in bamfile:
            if not is_primary(a):
                continue
            if a.query_name not in rname_to_read:
                rname_to_read[a.query_name] = a
                continue
            nosplice_match, nosplice_insert = a.get_cigar_stats()[0][0], a.get_cigar_stats()[0][1]
            withsplice_stats = rname_to_read[a.query_name].get_cigar_stats()[0]
            if nosplice_match > withsplice_stats[0] and nosplice_insert < withsplice_stats[1]:
                rname_to_read[a.query_name] = a
    return rname_to_read


def align_to_synthetic_genome(output, freadsname, min_fragment_size, threads):
    "OUTPUT.syntheticAligned.bam: each fusion read's better alignment to the synthetic genome"
    logging.info('aligning to synthetic fusion genome')
    synth_fa = output + '-syntheticFusionGenome.fa'
    withsplice, nosplice = output + '.syntheticAligned.withsplice.bam', output + '.syntheticAligned.nosplice.bam'
    align_to_synth_genome(synth_fa, freadsname, nosplice,
                          ['-s', str(min_fragment_size), '-t', str(threads), '-un', '--secondary=no', '-G', '1000k'])
    align_to_synth_genome(synth_fa, freadsname, withsplice,
                          ['-s', str(min_fragment_size), '-t', str(threads), '--secondary=no', '-G', '1000k'])
    rname_to_read = choose_synthetic_alignments(withsplice, nosplice)
    with pysam.AlignmentFile(withsplice, 'rb') as template, \
            pysam.AlignmentFile(output + '.syntheticAligned.unsorted.bam', 'wb', template=template) as outbam:
        for a in rname_to_read:
            outbam.write(rname_to_read[a])
    pipettor.run([('samtools', 'sort', '-o', output + '.syntheticAligned.bam', output + '.syntheticAligned.unsorted.bam')])
    pipettor.run([('samtools', 'index', output + '.syntheticAligned.bam')])
    pipettor.run([('rm', withsplice, withsplice + '.bai', nosplice, nosplice + '.bai',
                   output + '.syntheticAligned.unsorted.bam')])


def write_synthetic_junctions(output):
    "OUTPUT.syntheticAligned.SJ.bed: the junctions of the synthetic alignments to pass flair transcriptome"
    logging.info('getting ss')
    ipcmd = ('intron-prospector', f'--genome-fasta={output}-syntheticFusionGenome.fa', f'--intron-bed6={output}.syntheticAligned.IPSJ.bed', '-C', '0.0', '--sj-filter=all', f'{output}.syntheticAligned.bam')
    pipettor.run([ipcmd])

    fusiontobp = {}
    for bed in BedReader(f'{output}-syntheticBreakpointLoc.bed', numStdCols=3):
        fusiontobp[bed.chrom] = bed.chromStart

    crossing = count_breakpoint_reads(f'{output}.syntheticAligned.bam', fusiontobp)
    good_sj = select_junctions(BedReader(f'{output}.syntheticAligned.IPSJ.bed', numStdCols=6), fusiontobp, crossing)

    with open(f'{output}.syntheticAligned.SJ.bed', 'w') as out:
        for line in good_sj:
            out.write('\t'.join(line) + '\n')


def run_synthetic_transcriptome(output, sample_name, threads, frac_support):
    "flair transcriptome on the reads aligned to the synthetic fusion genome"
    # NOTE: only the chimeric reads supporting a fusion (fusionreads above) are
    # aligned to the synthetic genome, and this transcriptome call would throw out
    # any non-chimeric read anyway: with --fusion_breakpoints, count_sam_transcripts
    # only accepts an alignment that covers a fusion breakpoint (check_fusionbp).
    # Realigning all reads through flair fusion would need that changed as well.
    transcriptome_command = ['flair', 'transcriptome',
                             '--genome_aligned_bam', output + '.syntheticAligned.bam',
                             '--genome', output + '-syntheticFusionGenome.fa',
                             '--gtf', output + '-syntheticReferenceAnno.gtf',
                             '--ss_window', '8',
                             '--generate_map',
                             '--quality', '0',
                             '--sjc_support', '2',
                             '--end_window', '300',
                             '--frac_support', str(frac_support),
                             # '--no_check_splice',
                             # '--no_stringent',
                             '--no_align_to_annot',
                             '--fusion_breakpoints', output + '-syntheticBreakpointLoc.bed',
                             '--sample_name', sample_name,
                             '--output', output + '.syntheticAligned.flair',
                             '--threads', str(threads)]
    # only include junction if any where found
    junc_bed = output + '.syntheticAligned.SJ.bed'
    if os.path.exists(junc_bed) and (os.path.getsize(junc_bed) > 0):
        transcriptome_command.extend(['--junction_bed', junc_bed])

    logging.info('generating fusion transcriptome')
    logging.info(' '.join(transcriptome_command))
    pipettor.run(transcriptome_command)


def filter_fasta(in_fa, out_fa, names):
    "the records of in_fa named in names, to out_fa"
    with open(out_fa, 'w') as out:
        good = False
        for line in open(in_fa):
            if line[0] == '>':
                good = line[1:].rstrip() in names
            if good:
                out.write(line)


def filter_by_first_column(in_tsv, out_tsv, names, header=False):
    "the lines of in_tsv whose first column is in names, to out_tsv, after its header line if it has one"
    with open(out_tsv, 'w') as out:
        for i, line in enumerate(open(in_tsv)):
            if (header and i == 0) or line.split('\t', 1)[0] in names:
                out.write(line)


def write_final_fusions(output, genome, min_dist_between_breakpoints, genetoname, gtf_data):
    """the fusion isoforms in genome coordinates, as flair transcriptome writes its
    isoforms: bed, sequences, counts, reads, proteins, GTF and bigBed"""
    logging.info('converting coordinates from synthetic to genomic')
    iso_genes = convert_synthetic_isos(output + '.syntheticAligned.flair.isoforms.bed', output + '-syntheticBreakpointLoc.bed',
                                       output + '.fusions.isoforms.bed', min_dist_between_breakpoints, AnnotsByChrom(gtf_data))

    fusion_beds = list(BedReader(output + '.fusions.isoforms.bed', fixScores=True, bedClass=FlairBed))
    goodisos = {bed.name for bed in fusion_beds}
    filter_fasta(output + '.syntheticAligned.flair.isoforms.fa', output + '.fusions.isoforms.fa', goodisos)
    filter_by_first_column(output + '.syntheticAligned.flair.isoform.counts.txt',
                           output + '.fusions.isoform.counts.txt', goodisos)
    filter_by_first_column(output + '.syntheticAligned.flair.isoform.read.map.txt',
                           output + '.fusions.isoform.read.map.txt', goodisos)
    filter_by_first_column(output + '.syntheticAligned.flair.aaseq.tsv', output + '.fusions.aaseq.tsv',
                           {bed.aaseq_id for bed in fusion_beds if bed.aaseq_id is not None}, header=True)

    write_fusion_table(output, iso_genes, genetoname)

    bed_to_gtf(output + '.fusions.isoforms.bed', output + '.fusions.isoforms.gtf', is_flair_bed=True)

    with pysam.FastaFile(genome) as genome_fa:
        make_big_bed(genome_fa, output + '.syntheticAligned.flair.chrom.sizes', output + '.fusions.isoforms')


INTERMEDIATE_SUFFIXES = ('.annotated_transcripts.bed', '.annotated_transcripts.fa', '.fusionreads.prelim.fa',
                         '-syntheticFusionGenome.fa', '-syntheticFusionGenome.fa.fai', '-syntheticReferenceAnno.gtf',
                         '-syntheticBreakpointLoc.bed', '.syntheticAligned.bam', '.syntheticAligned.bam.bai',
                         '.syntheticAligned.IPSJ.bed', '.syntheticAligned.SJ.bed')


def remove_intermediates(output):
    "remove the files --keep_intermediate keeps"
    for filename in [output + suffix for suffix in INTERMEDIATE_SUFFIXES] \
            + glob.glob(glob.escape(output) + '.syntheticAligned.flair*'):
        if os.path.exists(filename):
            os.remove(filename)


def detectfusions(*, genome, gtf, genome_aligned_bam, sample_name, output, threads,
                  min_fragment_size, min_support, max_loci, max_dist_to_tss,
                  min_dist_between_breakpoints, keep_intermediate, frac_support=0.05, quality=10,
                  transcriptome_alignment=True, max_partners=1000):
    # fusion detection relies on stranded reads, as the help and docs say: a read's pieces
    # are ordered 5' to 3' by their position in the read
    # parsed once for everything that needs the annotation.  Without CDS records,
    # which nothing here uses, it takes a third less memory and time
    logging.info('reading gtf')
    gtf_data = gtf_data_parser(gtf, attrs=GtfAttrsSet.FLAIR, include_features=TRANSCRIPT_EXON_FEATURES)
    transcriptome_reads = make_chimeric_bams(genome, gtf_data, genome_aligned_bam, output, threads, min_fragment_size,
                                             transcriptome_alignment)
    annot = read_annotation(gtf_data)
    logging.info('finding preliminary fusions')
    # only the reads that needed aligning to transcripts, as a bam left by an earlier
    # run may hold others.  A read aligned to the genome in one piece is left out, and
    # with it a fusion from a deletion between genes closer than the genome aligner's
    # maximum intron length, which aligns to the genome as one spliced read
    prelim_fusions = id_chimeras(output + '.genomealigned.chim.bam', annot, min_support, max_loci=max_loci,
                                 max_dist_to_tss=max_dist_to_tss, max_partners=max_partners,
                                 transcriptome_bam=output + '.transcriptomealigned.chim.bam' if transcriptome_alignment else None,
                                 min_mapq=quality, transcriptome_reads=transcriptome_reads)
    fusionreads = write_prelim_fusions(prelim_fusions, output + '.prelimfusions.bed')
    if os.path.getsize(output + '.prelimfusions.bed') > 0:
        find_fusion_isoforms(fusionreads, gtf, gtf_data, genome, output, sample_name, threads, min_fragment_size,
                             min_dist_between_breakpoints, frac_support, transcriptome_alignment, annot.genetoname)
    else:
        report_nofusions(output)
    if not keep_intermediate:
        remove_intermediates(output)


def find_fusion_isoforms(fusionreads, gtf, gtf_data, genome, output, sample_name, threads, min_fragment_size,
                         min_dist_between_breakpoints, frac_support, transcriptome_alignment, genetoname):
    "the fusion isoforms of the preliminary fusions' reads, from the synthetic fusion genome"
    logging.info('obtaining fusion reads')
    freadsname = output + '.fusionreads.prelim.fa'
    chim_bams = [output + '.genomealigned.chim.bam']
    if transcriptome_alignment:
        chim_bams.append(output + '.transcriptomealigned.chim.bam')
    write_fusion_reads(fusionreads, chim_bams, freadsname)

    if not make_synthetic_genome(gtf, gtf_data, genome, output):
        report_nofusions(output)
        return
    align_to_synthetic_genome(output, freadsname, min_fragment_size, threads)
    write_synthetic_junctions(output)
    run_synthetic_transcriptome(output, sample_name, threads, frac_support)
    write_final_fusions(output, genome, min_dist_between_breakpoints, genetoname, gtf_data)
