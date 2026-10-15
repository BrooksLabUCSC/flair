#! /usr/bin/env python3

import os
import sys
import glob
import pipettor
import pysam
import logging
from flair.gtf_to_bed import gtf_to_bed
from flair.convert_synthetic_to_genome_bed import convert_synthetic_isos, get_paralog_ref, paralog_groups_for_genes
from flair.identify_prelim_fusions import id_chimeras, chimeric_read_names, genome_unsplit_reads
from flair import FlairInputDataError
from flair.gtf_io import gtf_record_parser, GtfAttrsSet
from flair.read_processing import get_sequence_from_bed
from flair.pycbio.hgdata.bed import Bed, BedReader
from flair.bed_to_gtf import bed_to_gtf
from flair.isoform_data import make_big_bed

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
                             'aligning reads to annotated transcripts.  Substantially faster, and -b can then '
                             'hold only the chimeric reads, but loses fusion reads whose partner piece the '
                             'genome alignment leaves clipped, particularly short 5\' exons')
    parser.add_argument('--keep_intermediate', action='store_true',
                        help='keep intermediate and temporary files for debugging purposes')
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

def def_value():
    return set()

def report_nofusions(outputprefix):
    logging.info('no fusions detected. Exiting')
    for file in [outputprefix + '.fusions.isoforms.bed', outputprefix + '.fusions.isoforms.fa']:
        f = open(file, 'w')
        f.close()

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


def detectfusions(*, genome, gtf, genome_aligned_bam, sample_name, output, threads,  # noqa: C901 - FIXME: reduce complexity
                  min_fragment_size, min_support, max_loci, max_dist_to_tss,
                  min_dist_between_breakpoints, keep_intermediate, frac_support=0.05, quality=10,
                  transcriptome_alignment=True, max_partners=1000):
    path = os.path.dirname(os.path.realpath(__file__)) + '/'

    # fusion detection relies on stranded reads, as the help and docs say: a read's pieces
    # are ordered 5' to 3' by their position in the read

    # Processing the gtf file so many times is really inefficient, how can we resolve this??
    genomechimbam = output + '.genomealigned.chim.bam'
    transcriptchimbam = output + '.transcriptomealigned.chim.bam'

    if not os.path.exists(genomechimbam):
        logging.info('getting chimeric reads from genome')
        infile = pysam.AlignmentFile(genome_aligned_bam, 'rb')
        outfile = pysam.AlignmentFile(genomechimbam, 'wb', template=infile)
        for align in infile:
            if align.is_mapped and not align.is_secondary:
                if align.has_tag('SA'):
                    outfile.write(align)
        infile.close()
        outfile.close()
        pysam.index(genomechimbam)

    if transcriptome_alignment and not os.path.exists(transcriptchimbam):
        logging.info('aligning to transcriptome and getting chimeric reads')
        annotated_bed = output + '.annotated_transcripts.bed'
        gtf_to_bed(annotated_bed, gtf, include_gene=True, name_sep='|')
        annotated_fa = output + '.annotated_transcripts.fa'
        get_sequence_from_bed(genome, output + '.annotated_transcripts.bed', annotated_fa)

        fa_cmd = ('samtools', 'fasta', genome_aligned_bam)
        mm2_cmd = ('minimap2', '-a', '-s', str(min_fragment_size), '-t', str(threads), '--secondary=no',
                   annotated_fa, '-')

        filter_cmd = ('samtools', 'view', '-hF', '0x104', '-e', '[SA] != ""')
        pipettor.run([fa_cmd, mm2_cmd, filter_cmd], stdout=output + '.transcriptomealigned.chim.bam')
        pipettor.run([('samtools', 'sort', '-o', output + '.transcriptomealigned.chim.sorted.bam', transcriptchimbam)])
        pipettor.run([('mv', output + '.transcriptomealigned.chim.sorted.bam', transcriptchimbam)])
        pysam.index(transcriptchimbam)

    logging.info('reading gtf')
    genetoinfo, genetoexons, genetoname = {}, {}, {}
    chrom_to_gene_pos = {}
    gene_to_all_exons, juncs_to_gene = {}, {}
    for rec in gtf_record_parser(gtf, include_features={'gene', 'exon', 'transcript'}, attrs=GtfAttrsSet.ALL):
        gene_id = rec.gene_id
        if rec.feature == 'gene':
            genetoinfo[gene_id] = [rec.chrom, rec.start, rec.end, rec.strand, []]
            genetoname[gene_id] = rec.gene_name if rec.gene_name else gene_id
            if rec.chrom not in chrom_to_gene_pos:
                chrom_to_gene_pos[rec.chrom] = []
                juncs_to_gene[rec.chrom] = {}
            chrom_to_gene_pos[rec.chrom].append((rec.start, rec.end, rec.strand, gene_id))

        elif rec.feature == 'exon':
            if gene_id not in genetoexons:
                genetoexons[gene_id] = {}
            if rec.transcript_id not in genetoexons[gene_id]:
                genetoexons[gene_id][rec.transcript_id] = []
            genetoexons[gene_id][rec.transcript_id].append((rec.start, rec.end))

        elif rec.feature == 'transcript':
            end5 = rec.start if rec.strand == '+' else rec.end
            genetoinfo[gene_id][-1].append(end5)

    logging.info('continuing to parse annot')
    # FOR JUNCS TO GENE, DO BY CHROM AS WELL
    for chrom in chrom_to_gene_pos:
        chrom_to_gene_pos[chrom] = sorted(chrom_to_gene_pos[chrom])

    for gene in genetoexons:
        gene_to_all_exons[gene] = set()
        chrom = genetoinfo[gene][0]
        for t in genetoexons[gene]:
            exons = sorted(list(genetoexons[gene][t]))
            gene_to_all_exons[gene].update(set(exons))
            if len(exons) > 1:
                juncs = [(exons[x][1], exons[x + 1][0]) for x in range(len(exons) - 1)]
                for j in juncs:
                    if j not in juncs_to_gene[chrom]:
                        juncs_to_gene[chrom][j] = set()
                    juncs_to_gene[chrom][j].add(gene)
        allexons = sorted(list(gene_to_all_exons[gene]))
        newexons = []
        laststart, lastend = -1, -1
        for s, e in allexons:
            if s > lastend:
                if laststart != -1:
                    newexons.append((laststart, lastend))
                laststart, lastend = s, e
            else:
                lastend = max(lastend, e)
        if laststart != -1:
            newexons.append((laststart, lastend))
        gene_to_all_exons[gene] = newexons

    intronLocs, intronToGenome = {}, {}
    for g in genetoexons:
        chrom, start, end, strand, _ = genetoinfo[g]
        for t in genetoexons[g]:
            myexons = sorted(genetoexons[g][t])
            mylocs = [[0, myexons[0][0] - 500, myexons[0][0]]]  # add start of transcript
            runningtot = 0
            for i in range(len(myexons) - 1):
                runningtot += myexons[i][1] - myexons[i][0]  # add size of last exon
                mylocs.append([runningtot, myexons[i][1], myexons[i + 1][0]])  # add intron
            runningtot += myexons[-1][1] - myexons[-1][0]
            mylocs.append([runningtot, myexons[-1][1], myexons[-1][1] + 500])
            if strand == '-':
                mylocs = [[runningtot - mylocs[x][0], mylocs[x][1], mylocs[x][2]] for x in range(len(mylocs))]
            intronLocs[t] = sorted([x[0] for x in mylocs])
            intronToGenome[t] = {x[0]: (x[1], x[2]) for x in mylocs}

    gene_to_paralogs = paralog_groups_for_genes(
        get_paralog_ref(os.path.realpath(__file__).split('flair_fusion')[0] + 'dgd_Hsa_all_v71.tsv'), genetoinfo)

    tchim = {}
    if transcriptome_alignment:
        logging.info('loading transcriptomic chimeras')
        # reads chimeric aligned to transcripts but whole in the genome alignment.  A fusion
        # from a deletion between genes closer than the genome aligner's maximum intron
        # length aligns to the genome as one spliced read, so is left out too
        unsplit = genome_unsplit_reads(genome_aligned_bam, chimeric_read_names(transcriptchimbam), min_fragment_size)
        logging.info(f'transcriptomic chimeras: {len(unsplit)} reads aligned to the genome in one piece left out')
        tchim = id_chimeras('transcriptomic', transcriptchimbam, genetoinfo, chrom_to_gene_pos,
                            gene_to_all_exons, juncs_to_gene, gene_to_paralogs, genetoname,
                            min_support, maxloci=max_loci, reqdisttostart=max_dist_to_tss,
                            maxpromiscuity=max_partners, intronLocs=intronLocs, intronToGenome=intronToGenome,
                            skip_reads=unsplit)
    logging.info('loading genomic chimeras')

    combchim = id_chimeras('genomic', genomechimbam, genetoinfo, chrom_to_gene_pos, gene_to_all_exons, juncs_to_gene, gene_to_paralogs, genetoname, min_support, maxloci=max_loci, reqdisttostart=max_dist_to_tss, maxpromiscuity=max_partners,
                            min_mapq=quality)

    logging.info('combining genomic and transcriptomic')
    #
    for f in tchim:
        if f not in combchim:
            combchim[f] = tchim[f]
        else:
            combchim[f]['reads'] = combchim[f]['reads'] | tchim[f]['reads']
            combchim[f]['disttostart'].extend(tchim[f]['disttostart'])
            combchim[f]['qdist'].extend(tchim[f]['qdist'])
            genes = f.split('__')
            for g in genes:
                gends, tends = combchim[f][g][1:], tchim[f][g][1:]
                if gends[0] < gends[1]:
                    combchim[f][g][1] = min(gends[0], tends[0])
                    combchim[f][g][2] = max(gends[1], tends[1])
                else:
                    combchim[f][g][1] = max(gends[0], tends[0])
                    combchim[f][g][2] = min(gends[1], tends[1])

    bedout = open(output + '.prelimfusions.bed', 'w')
    fusionreads = set()
    for f in combchim:
        fusionreads.update(combchim[f]['reads'])
        genes = f.split('__')
        for g in genes:
            startdiststr = ','.join([str(x) for x in combchim[f]['disttostart']])
            bedname = g + '__' + '--'.join(genes)
            strand = '+' if combchim[f][g][1] < combchim[f][g][2] else '-'
            positions = sorted(combchim[f][g][1:])
            Bed(combchim[f][g][0], positions[0], positions[1],
                name=bedname, score=len(combchim[f]['reads']), strand=strand,
                extraCols=[startdiststr]).write(bedout)

    bedout.close()

    if os.path.getsize(output + '.prelimfusions.bed') == 0:
        report_nofusions(output)
        return

    logging.info('obtaining fusion reads')

    seenreads = set()
    freadsname = output + '.fusionreads.prelim.fa'
    out_fa = open(freadsname, 'w')
    bamfile = pysam.AlignmentFile(genomechimbam, 'rb')
    for a in bamfile:
        if not a.is_secondary and not a.is_supplementary and a.query_name in fusionreads:
            out_fa.write('>' + a.query_name + '\n' + a.get_forward_sequence() + '\n')
            seenreads.add(a.query_name)
    bamfile.close()
    if transcriptome_alignment:
        bamfile = pysam.AlignmentFile(transcriptchimbam, 'rb')
        for a in bamfile:
            if not a.is_secondary and not a.is_supplementary and a.query_name in fusionreads and a.query_name not in seenreads:
                out_fa.write('>' + a.query_name + '\n' + a.get_forward_sequence() + '\n')
        bamfile.close()
    out_fa.close()

    logging.info('generating synthetic reference')

    makesynthcommand = [sys.executable, path + 'make_synthetic_fusion_reference.py', '-a', gtf, '-g', genome,
                        '-o', output, '-c', output + '.prelimfusions.bed']
    pipettor.run([makesynthcommand])
    if os.path.getsize(output + '-syntheticFusionGenome.fa') == 0:
        report_nofusions(output)
        return

    # FIXME: pipettor by default captures stderr to include in an error message, this hides logging from
    # lower level.  Maybe don't capture when running flair subtools
    faidxcommand = ['samtools', 'faidx', output + '-syntheticFusionGenome.fa']
    pipettor.run([faidxcommand])

    logging.info('aligning to synthetic fusion genome')
    align_to_synth_genome(output + '-syntheticFusionGenome.fa', freadsname, output + '.syntheticAligned.nosplice.bam',
                          ['-s', str(min_fragment_size), '-t', str(threads), '-un', '--secondary=no', '-G', '1000k'])
    align_to_synth_genome(output + '-syntheticFusionGenome.fa', freadsname, output + '.syntheticAligned.withsplice.bam',
                          ['-s', str(min_fragment_size), '-t', str(threads), '--secondary=no', '-G', '1000k'])

    rname_to_read = {}
    with pysam.AlignmentFile(output + '.syntheticAligned.withsplice.bam', 'rb') as bamfile:
        for a in bamfile:
            if a.is_mapped and not a.is_secondary and not a.is_supplementary:
                rname_to_read[a.query_name] = a
    with pysam.AlignmentFile(output + '.syntheticAligned.nosplice.bam', 'rb') as bamfile:
        for a in bamfile:
            if a.is_mapped and not a.is_secondary and not a.is_supplementary:
                if a.query_name not in rname_to_read:
                    rname_to_read[a.query_name] = a
                else:
                    nosplice_match = a.get_cigar_stats()[0][0]
                    nosplice_insert = a.get_cigar_stats()[0][1]
                    withsplice_match = rname_to_read[a.query_name].get_cigar_stats()[0][0]
                    withsplice_insert = rname_to_read[a.query_name].get_cigar_stats()[0][1]
                    if nosplice_match > withsplice_match and nosplice_insert < withsplice_insert:
                        rname_to_read[a.query_name] = a
    template = pysam.AlignmentFile(output + '.syntheticAligned.withsplice.bam', 'rb')
    outbam = pysam.AlignmentFile(output + '.syntheticAligned.unsorted.bam', 'wb', template=template)
    for a in rname_to_read:
        outbam.write(rname_to_read[a])
    template.close()
    outbam.close()

    pipettor.run([('samtools', 'sort', '-o', output + '.syntheticAligned.bam', output + '.syntheticAligned.unsorted.bam')])
    pipettor.run([('samtools', 'index', output + '.syntheticAligned.bam')])
    pipettor.run([('rm', output + '.syntheticAligned.withsplice.bam', output + '.syntheticAligned.withsplice.bam.bai',
                   output + '.syntheticAligned.nosplice.bam', output + '.syntheticAligned.nosplice.bam.bai',
                   output + '.syntheticAligned.unsorted.bam')])

    logging.info('getting ss')
    ipcmd = ('intron-prospector', f'--genome-fasta={output}-syntheticFusionGenome.fa', f'--intron-bed6={output}.syntheticAligned.IPSJ.bed', '-C', '0.0', '--sj-filter=all', f'{output}.syntheticAligned.bam')
    pipettor.run([ipcmd])

    fusiontobp = {}
    for bed in BedReader(f'{output}-syntheticBreakpointLoc.bed', numStdCols=3):
        fusiontobp[bed.chrom] = bed.chromStart

    crossing = count_breakpoint_reads(f'{output}.syntheticAligned.bam', fusiontobp)
    good_sj = select_junctions(BedReader(f'{output}.syntheticAligned.IPSJ.bed', numStdCols=6), fusiontobp, crossing)

    out = open(f'{output}.syntheticAligned.SJ.bed', 'w')
    for line in good_sj:
        out.write('\t'.join(line) + '\n')
    out.close()

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

    logging.info('converting coordinates from synthetic to genomic')
    convert_synthetic_isos(output + '.syntheticAligned.flair.isoforms.bed',
                           output + '.syntheticAligned.flair.isoform.read.map.txt', freadsname,
                           output + '-syntheticBreakpointLoc.bed', output + '.fusions.isoforms.bed', min_dist_between_breakpoints)

    goodisos = set()
    for bed in BedReader(output + '.fusions.isoforms.bed', fixScores=True):
        goodisos.add(bed.name)

    with open(output + '.fusions.isoforms.fa', 'w') as out:
        good = False
        for line in open(output + '.syntheticAligned.flair.isoforms.fa'):
            if line[0] == '>':
                name = line[1:].rstrip()
                if name in goodisos:
                    good = True
                else:
                    good = False
            if good:
                out.write(line)

    with open(output + '.fusion.isoform.read.map.txt', 'w') as out:
        for line in open(output + '.syntheticAligned.flair.isoform.read.map.txt'):
            name = line.split('\t', 1)[0]
            if name in goodisos:
                out.write(line)

    bed_to_gtf(output + '.fusions.isoforms.bed', output + '.fusions.isoforms.gtf', is_flair_bed=True)

    with pysam.FastaFile(genome) as genome_fa:
        make_big_bed(genome_fa, output + '.syntheticAligned.flair.chrom.sizes', output + '.fusions.isoforms')

    # removing extra FLAIR files
    if not keep_intermediate:
        for filename in glob.glob(output + '.syntheticAligned.flair*'):
            os.remove(filename)
