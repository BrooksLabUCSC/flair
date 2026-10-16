"""
The table of final fusions, OUTPUT.fusions.tsv: per fusion, its preliminary and final
read counts, its isoforms, and its best supported breakpoint and isoform, sorted by
final reads.
"""
from collections import defaultdict

FUSION_TABLE_COLUMNS = ('fusion', 'gene_ids', 'prelim_reads', 'final_reads', 'isoforms',
                        'top_breakpoint', 'top_breakpoint_reads', 'top_isoform', 'top_isoform_reads')


def prelim_read_counts(prelim_bed):
    "fusion (its synthetic contig name) -> its reads in the preliminary fusions"
    counts = {}
    for line in open(prelim_bed):
        f = line.split('\t')
        counts[f[3].split('__', 1)[1]] = int(f[4])
    return counts


def isoform_contigs(synthetic_isoforms_bed):
    "isoform -> the synthetic fusion contig it is on"
    contigs = {}
    for line in open(synthetic_isoforms_bed):
        f = line.split('\t')
        contigs[f[3]] = f[0]
    return contigs


def isoform_read_counts(read_map):
    "isoform -> its number of reads"
    counts = {}
    for line in open(read_map):
        iso, reads = line.rstrip('\n').split('\t')
        counts[iso] = len(reads.split(','))
    return counts


def isoform_pieces(fusion_isoforms_bed):
    "isoform -> its pieces as (position in fusion, chrom, start, end, strand)"
    pieces = defaultdict(list)
    for line in open(fusion_isoforms_bed):
        f = line.rstrip('\n').split('\t')
        pieces[f[3]].append((int(f[20]), f[0], int(f[1]), int(f[2]), f[5]))
    return pieces


def breakpoint(piece5, piece3):
    """the junction of a 5' piece and the 3' piece after it, as chrom:pos(strand) of the
    5' piece's last base and of the 3' piece's first base, 1-based, joined by --"""
    _, chrom5, start5, end5, strand5 = piece5
    _, chrom3, start3, end3, strand3 = piece3
    donor = end5 if strand5 == '+' else start5 + 1
    acceptor = start3 + 1 if strand3 == '+' else end3
    return f'{chrom5}:{donor}({strand5})--{chrom3}:{acceptor}({strand3})'


def isoform_breakpoints(pieces):
    "an isoform's junctions between consecutive loci, in fusion order, comma-separated"
    pieces = sorted(pieces)
    return ','.join(breakpoint(pieces[k], pieces[k + 1]) for k in range(len(pieces) - 1))


def fusion_row(genes, isos, prelim_reads, iso_contig, iso_reads, iso_pieces, genetoname):
    """a fusion's row of the table, from its final isoforms, with the preliminary
    reads of the synthetic contigs they are on"""
    bp_reads = defaultdict(int)
    for iso in isos:
        bp_reads[isoform_breakpoints(iso_pieces[iso])] += iso_reads[iso]
    top_iso = min(isos, key=lambda iso: (-iso_reads[iso], iso))
    top_bp = min(bp_reads, key=lambda bp: (-bp_reads[bp], bp))
    contigs = {iso_contig[iso] for iso in isos}
    return ['--'.join(genetoname.get(g, g) for g in genes), '--'.join(genes),
            sum(prelim_reads.get(contig, 0) for contig in contigs),
            sum(iso_reads[iso] for iso in isos), len(isos), top_bp, bp_reads[top_bp], top_iso, iso_reads[top_iso]]


def fusion_rows(prelim_reads, iso_contig, iso_genes, iso_reads, iso_pieces, genetoname):
    "the table's rows, a fusion for each set of final genes, most final reads first"
    by_genes = defaultdict(list)
    for iso in iso_pieces:
        by_genes[iso_genes[iso]].append(iso)
    rows = [fusion_row(genes, isos, prelim_reads, iso_contig, iso_reads, iso_pieces, genetoname)
            for genes, isos in by_genes.items()]
    return sorted(rows, key=lambda row: (-row[3], row[0]))


def write_table(rows, tsv):
    with open(tsv, 'w') as out:
        out.write('\t'.join(FUSION_TABLE_COLUMNS) + '\n')
        for row in rows:
            out.write('\t'.join(str(x) for x in row) + '\n')


def write_fusion_table(output, iso_genes, genetoname):
    """OUTPUT.fusions.tsv, from the run's preliminary fusions and final isoforms, with
    iso_genes, isoform -> its final genes (convert_synthetic_isos)"""
    rows = fusion_rows(prelim_read_counts(output + '.prelimfusions.bed'),
                       isoform_contigs(output + '.syntheticAligned.flair.isoforms.bed'), iso_genes,
                       isoform_read_counts(output + '.fusions.isoform.read.map.txt'),
                       isoform_pieces(output + '.fusions.isoforms.bed'),
                       genetoname)
    write_table(rows, output + '.fusions.tsv')
