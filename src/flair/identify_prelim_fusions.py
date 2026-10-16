#!/usr/bin/env python3

import logging
import pysam
from collections import Counter, defaultdict
from functools import partial
from statistics import median
from flair.convert_synthetic_to_genome_bed import identify_fusion_problems, fusion_role


def binarySearch(arr, t):
    """Nearest value in a sorted list.  The loop used to fall off the end returning
    None when i and j converged without a branch firing, and the caller then indexed
    with it."""
    if t <= arr[0]:
        return arr[0]
    if t >= arr[-1]:
        return arr[-1]
    i, j, mid = 0, len(arr) - 1, 0
    while i < j:
        mid = int((i + j) / 2)
        if arr[mid] == t:
            return arr[mid]
        elif t < arr[mid]:
            if mid > 0 and t > arr[mid - 1]:
                if abs(arr[mid] - t) < abs(arr[mid - 1] - t):
                    return arr[mid]
                else:
                    return arr[mid - 1]
            j = mid
        else:
            if mid < len(arr) - 1 and t < arr[mid + 1]:
                if abs(arr[mid] - t) < abs(arr[mid + 1] - t):
                    return arr[mid]
                else:
                    return arr[mid + 1]
            i = mid + 1
    return arr[i]


def getGenomicPreciseLoc(tname, bpCoord, genedir, intronLocs, intronToGenome):
    closestSS = binarySearch(intronLocs[tname], bpCoord)
    bpIntronEnds = intronToGenome[tname][closestSS]
    genomeSS = None
    diffFromSS = bpCoord - closestSS
    if closestSS == 0:  # start of gene
        if genedir == '+':
            genomeSS = bpIntronEnds[1]
        else:
            genomeSS = bpIntronEnds[0]
    elif closestSS == max(intronLocs[tname]):  # end of gene
        if genedir == '+':
            genomeSS = bpIntronEnds[0]
        else:
            genomeSS = bpIntronEnds[1]
    elif genedir == '+':
        if bpCoord >= closestSS:
            genomeSS = bpIntronEnds[1]
        else:
            genomeSS = bpIntronEnds[0]
    else:
        if bpCoord >= closestSS:
            genomeSS = bpIntronEnds[0]
        else:
            genomeSS = bpIntronEnds[1]
    genomepos = None
    if genedir == '+':
        genomepos = genomeSS + diffFromSS
    else:
        genomepos = genomeSS - diffFromSS
    return genomepos


def get_exon_intron_blocks(read):
    align_start = read.reference_start
    align_end = read.reference_end
    ref_pos = align_start
    intron_blocks = []
    has_match = False
    for block in read.cigartuples:
        if block[0] == 3:  # intron
            if has_match:
                intron_blocks.append([ref_pos, ref_pos + block[1]])
            # this fixes weird bug if there's an intron, then an insertion, then another intron???
            elif len(intron_blocks) > 0:
                intron_blocks[-1][1] += block[1]
            has_match = False
            ref_pos += block[1]
        elif block[0] in {0, 7, 8, 2}:  # consumes reference
            ref_pos += block[1]
            if block[0] in {0, 7, 8}:
                has_match = True
    intron_blocks = [tuple(x) for x in intron_blocks]
    if len(intron_blocks) == 0:
        exon_blocks = [(align_start, align_end),]
    else:
        exon_blocks = [(align_start, intron_blocks[0][0])] + [(intron_blocks[x][1], intron_blocks[x + 1][0]) for x in range(len(intron_blocks) - 1)] + [(intron_blocks[-1][1], align_end)]
    return intron_blocks, exon_blocks


def getCorrectGene(chrom_to_gene_pos, gene_to_all_exons, juncs_to_gene, chrom, readblocks, thisdir):
    intron_blocks, exon_blocks = readblocks
    if chrom not in chrom_to_gene_pos:
        my_gene = (chrom, exon_blocks[0][0], exon_blocks[-1][1])  # chrom + ':' + str(round(exon_blocks[0][0], -4))
    else:
        found_genes = []
        for j in intron_blocks:
            if j in juncs_to_gene[chrom]:
                found_genes.extend(list(juncs_to_gene[chrom][j]))
        if len(found_genes) > 0:
            my_gene = Counter(found_genes).most_common()[0][0]
        else:
            # a gene on the piece's strand that its aligned blocks overlap the exons of
            # by the most bases.  A piece in a gene's introns, or spliced across its exons
            # without aligning to them, is not the gene's: it's a non-genic locus
            s, e = exon_blocks[0][0], exon_blocks[-1][1]
            gene_overlaps = []
            for start, end, strand, gene in chrom_to_gene_pos[chrom]:
                if start > e:  # assumes sorted
                    break
                if strand == thisdir and min(e, end) > max(start, s):  # require strand match
                    totoverlap = exonic_overlap(exon_blocks, gene_to_all_exons.get(gene, ()))
                    if totoverlap > 0:
                        gene_overlaps.append((totoverlap, gene))
            if len(gene_overlaps) > 0:
                gene_overlaps.sort(reverse=True)
                my_gene = gene_overlaps[0][1]
            else:
                my_gene = (chrom, exon_blocks[0][0], exon_blocks[-1][1])  # chrom + ':' + str(round(exon_blocks[0][0], -4))
    return my_gene


def name_nongenic_group(ognongenic_to_genes, group, chrom, start, end):
    "name every locus in an overlapping group after the group's merged span"
    gname = f'{chrom}:{start}-{end}'
    for n in group:
        ognongenic_to_genes[n] = gname


def group_nongenic_loci(nongenicloci):
    """Map each non-genic locus to the merged span of the overlapping loci it belongs
    to.  nongenicloci must be sorted."""
    ognongenic_to_genes = {}
    lastchrom, laststart, lastend, group = None, -1, -1, []
    for chrom, s, e in nongenicloci:
        if s > lastend or chrom != lastchrom:
            name_nongenic_group(ognongenic_to_genes, group, lastchrom, laststart, lastend)
            group = []
            lastchrom, laststart, lastend = chrom, s, e
        elif e > lastend:
            lastend = e
        group.append((chrom, s, e))
    name_nongenic_group(ognongenic_to_genes, group, lastchrom, laststart, lastend)
    return ognongenic_to_genes


def exonic_overlap(blocks, exons):
    """Bases of the sorted blocks that fall in the sorted, merged exons."""
    total, i = 0, 0
    for s, e in blocks:
        while i < len(exons) and exons[i][1] <= s:
            i += 1
        j = i
        while j < len(exons) and exons[j][0] < e:
            total += min(e, exons[j][1]) - max(s, exons[j][0])
            j += 1
    return total


def exon_sharing_genes(chrom_to_gene_pos, gene_to_all_exons):
    """gene -> the genes on its strand sharing exon sequence with it, as a gene and a
    readthrough annotation containing it.  chrom_to_gene_pos is chrom -> sorted
    (start, end, strand, gene), gene_to_all_exons gene -> sorted merged exons."""
    sharing = defaultdict(set)
    for genes in chrom_to_gene_pos.values():
        active = []
        for start, end, strand, gene in genes:
            active = [a for a in active if a[1] > start]
            for _, _, strand2, gene2 in active:
                if strand2 == strand and exonic_overlap(gene_to_all_exons.get(gene, ()), gene_to_all_exons.get(gene2, ())) > 0:
                    sharing[gene].add(gene2)
                    sharing[gene2].add(gene)
            active.append((start, end, strand, gene))
    return sharing


# a shorter gene is taken if its exons hold this much of what the assigned gene's do
SHORTER_GENE_MIN_COVER = 0.95


def prefer_shorter_gene(gene, blocks, sharing, genetoinfo, gene_to_all_exons):
    """The gene for a read piece aligned to blocks (sorted genomic intervals): of gene
    and the genes sharing exons with it, the shortest whose exons cover the blocks
    about as well as gene's, so a readthrough or other long annotation gives way to
    the gene it contains.  Otherwise a fusion's reads are split between the two,
    each needing min_support on its own, and both are reported."""
    if gene not in sharing:
        return gene
    cover = exonic_overlap(blocks, gene_to_all_exons[gene])
    if cover == 0:  # nothing to compare: a gene covering none of it is no better
        return gene
    best, best_len = gene, genetoinfo[gene][2] - genetoinfo[gene][1]
    for other in sorted(sharing[gene]):
        other_len = genetoinfo[other][2] - genetoinfo[other][1]
        if other_len < best_len and exonic_overlap(blocks, gene_to_all_exons[other]) >= SHORTER_GENE_MIN_COVER * cover:
            best, best_len = other, other_len
    return best


# chimeric pieces on opposite strands this close are a fold-back, not a fusion
FOLDBACK_MAX_DIST = 1000

# a preliminary fusion's loci, which make its synthetic contig, are padded by this
# beyond its reads
LOCUS_PADDING = 1000


def is_foldback(pieces):
    """Is a chimeric read a fold-back: two of its pieces on opposite strands of one
    chromosome overlapping or within FOLDBACK_MAX_DIST.  cDNA that folds back on
    itself gives a read of a region joined to its own reverse complement, and the
    piece on the other strand would be assigned to an antisense gene there, as
    SLC7A5--SLC7A5-AS1 in MCF7.  pieces are id_chimeras' pieces, [(qstart, refpos),
    (qend, refpos), gene, strand, chrom]."""
    spans = [(piece[4], piece[3], min(piece[0][1], piece[1][1]), max(piece[0][1], piece[1][1])) for piece in pieces]
    for i, (chrom, strand, start, end) in enumerate(spans):
        for chrom2, strand2, start2, end2 in spans[i + 1:]:
            if chrom2 == chrom and strand2 != strand and max(start, start2) - min(end, end2) <= FOLDBACK_MAX_DIST:
                return True
    return False


def is_low_quality(read, min_mapq):
    """Is this piece of a chimeric read below min_mapq: a piece that pieces about as well
    elsewhere places a partner the read may not have, as the multi-copy subtelomeric
    DDX11L genes, or TAOK3--RPS2 in MCF7 at MAPQ 1 to 15."""
    return min_mapq > 0 and read.is_mapped and not read.is_secondary and read.mapping_quality < min_mapq


# an annotated transcript start this far downstream of a read's 5' end still counts as
# its start, for reads running a little past the annotated TSS
TSS_DOWNSTREAM_BUFFER = 500


def dist_to_upstream_tss(end5, strand, tstarts, buffer=TSS_DOWNSTREAM_BUFFER):
    """Distance from a 5' gene's read 5' end to the nearest of the gene's transcript
    starts upstream of it on the gene's strand, or None if there is none.  A start
    downstream of the read's 5' end can't be its promoter: for reads beginning in ETV6
    exon 5 (SeraCare ETV6--NTRK3), the nearest start was an alternative one 16 kb
    downstream.  Starts up to buffer downstream count, at distance 0."""
    if strand == '+':
        dists = [end5 - t for t in tstarts if t <= end5 + buffer]
    else:
        dists = [t - end5 for t in tstarts if t >= end5 - buffer]
    return max(0, min(dists)) if dists else None


def fusion_tss_distance(ends5, strand, tstarts):
    """The smallest distance from any of a fusion's 5' ends (the most upstream one and
    the median one) to an upstream transcript start of its 5' gene, or None."""
    dists = [d for d in (dist_to_upstream_tss(e, strand, tstarts) for e in ends5) if d is not None]
    return min(dists) if dists else None


def needs_transcriptome(read, min_clip):
    """Can aligning this read to annotated transcripts find a fusion the genome alignment
    doesn't: it is unmapped, chimeric, or leaves min_clip or more unaligned at an end.  A
    read aligned to the genome in one piece, chimeric in the transcript alignment, is the
    transcript alignment putting consecutive sequence on two annotated genes, as
    overlapping or nested genes (MACC1--MACC1-OT1 in SeraCare) or segmental duplications
    (UPK3B in MCF7).  A clipped end can be a short partner piece minimap2 didn't align
    as supplementary (MYO6--SENP6 in MCF7)."""
    if read.is_unmapped or read.has_tag('SA'):
        return True
    cigar = read.cigartuples
    clip5 = cigar[0][1] if cigar[0][0] in (4, 5) else 0
    clip3 = cigar[-1][1] if cigar[-1][0] in (4, 5) else 0
    return clip5 >= min_clip or clip3 >= min_clip


def scan_genome_alignment(genome_bam, min_clip, chim_bam=None, reads_fa=None):
    """One pass over the genome alignment: writes the alignments of its chimeric reads to
    chim_bam, and the reads that need_transcriptome to reads_fa, each skipped when None,
    and returns the names of the latter.  Only those reads are aligned to annotated
    transcripts, a small share of them: 11% in MCF7, 1.5% in SeraCare."""
    needed = set()
    with pysam.AlignmentFile(genome_bam, 'rb') as pieces:
        chim_out = pysam.AlignmentFile(chim_bam, 'wb', template=pieces) if chim_bam else None
        fa_out = open(reads_fa, 'w') if reads_fa else None
        for read in pieces.fetch(until_eof=True):
            if read.is_secondary:
                continue
            if chim_out is not None and read.is_mapped and read.has_tag('SA'):
                chim_out.write(read)
            if read.is_supplementary or not needs_transcriptome(read, min_clip):
                continue
            needed.add(read.query_name)
            if fa_out is not None and read.query_sequence:
                fa_out.write('>' + read.query_name + '\n' + read.get_forward_sequence() + '\n')
        if chim_out is not None:
            chim_out.close()
        if fa_out is not None:
            fa_out.close()
    return needed


def query_span(read):
    "a piece's start and end in the read, counting hard-clipped bases"
    qstart, qend = read.query_alignment_start, read.query_alignment_end
    cigar = read.cigartuples
    if cigar[0][0] == 5:  # just hard clipping
        qstart += cigar[0][1]
        qend += cigar[0][1]
    return qstart, qend


def genomic_piece(read, chrom_to_gene_pos, gene_to_all_exons, juncs_to_gene, sharing, genetoinfo):
    """[(read 5' end, its position), (read 3' end, its position), gene, strand, chrom] of
    a piece of a read aligned to the genome, in read coordinates"""
    refstart, refend = read.reference_start, read.reference_end
    readdir = '-' if read.is_reverse else '+'
    qstart, qend = query_span(read)
    refchr = read.reference_name
    readblocks = get_exon_intron_blocks(read)
    genename = getCorrectGene(chrom_to_gene_pos, gene_to_all_exons, juncs_to_gene, refchr, readblocks, readdir)
    if type(genename) is not tuple:  # tuples are non-genic loci
        genename = prefer_shorter_gene(genename, readblocks[1], sharing, genetoinfo, gene_to_all_exons)
    if readdir == '+':
        return [(qstart, refstart), (qend, refend), genename, readdir, refchr]
    readlen = read.infer_read_length()
    return [(readlen - qend, refend), (readlen - qstart, refstart), genename, readdir, refchr]


def transcriptomic_piece(read, genetoinfo, intronLocs, intronToGenome, sharing, gene_to_all_exons):
    "a piece of a read aligned to an annotated transcript, as genomic_piece, in genome coordinates"
    readdir = '-' if read.is_reverse else '+'
    qstart, qend = query_span(read)
    tname, genename = read.reference_name.split('|')
    refchr, genedir = genetoinfo[genename][0], genetoinfo[genename][3]  # can do this because already required the read to be forward strand
    refstart = getGenomicPreciseLoc(tname, read.reference_start, genedir, intronLocs, intronToGenome)
    refend = getGenomicPreciseLoc(tname, read.reference_end, genedir, intronLocs, intronToGenome)
    refstart, refend = min(refstart, refend), max(refstart, refend)
    genename = prefer_shorter_gene(genename, [(refstart, refend)], sharing, genetoinfo, gene_to_all_exons)
    if genedir == '+':
        return [(qstart, refstart), (qend, refend), genename, readdir, refchr]
    return [(qstart, refend), (qend, refstart), genename, genedir, refchr]


def read_pieces(mode, bam, piece_of, min_mapq, only_reads):
    """read -> its chimeric pieces, made by piece_of.  Aligned to transcripts, a piece must
    be on the transcript's strand.  Mapping quality is only used for genome alignments:
    aligned to transcripts, a read of a gene with several overlapping isoforms has a low
    one regardless.  Reads with a low piece are left out after the pass, as the piece can
    come after the others."""
    low_quality = set()
    read_to_pieces = {}
    with pysam.AlignmentFile(bam, "rb") as bam_file:
        for read in bam_file:
            if mode == 'genomic' and is_low_quality(read, min_mapq):
                low_quality.add(read.query_name)
            if read.is_mapped and not read.is_secondary and read.has_tag('SA') \
                    and (mode == 'genomic' or not read.is_reverse) \
                    and (only_reads is None or read.query_name in only_reads):
                read_to_pieces.setdefault(read.query_name, []).append(piece_of(read))
    for read in low_quality:
        read_to_pieces.pop(read, None)
    if low_quality:
        logging.info(f'{mode} chimeras: {len(low_quality)} reads with a piece below MAPQ {min_mapq} left out')
    return read_to_pieces


def drop_foldbacks(mode, read_to_pieces):
    foldbacks = [read for read, pieces in read_to_pieces.items() if is_foldback(pieces)]
    for read in foldbacks:
        del read_to_pieces[read]
    logging.info(f'{mode} chimeras: {len(foldbacks)} fold-back reads of {len(read_to_pieces) + len(foldbacks)} left out')


def consensus_gene(genes):
    """the gene most of a locus's pieces were given, an annotated gene over a non-genic
    locus (chrom:start-end) on a tie, then the first by name"""
    counts = Counter(genes)
    return min(counts, key=lambda g: (-counts[g], ':' in g, g))


def _locus_extents(reads, read_to_pieces):
    "per fusion position, the (chrom, strand, start, end) its reads' pieces span"
    extents = []
    for i in range(len(read_to_pieces[reads[0]])):
        pieces = [sorted(read_to_pieces[read])[i] for read in reads]
        extents.append((pieces[0][4], pieces[0][3], min(min(p[0][1], p[1][1]) for p in pieces),
                        max(max(p[0][1], p[1][1]) for p in pieces)))
    return extents


def _same_loci(extents1, extents2):
    """do two fusions have their loci in the same places: each position's on one
    strand, overlapping once padded by LOCUS_PADDING, as their synthetic contigs would
    be, which would then hold the same sequence"""
    return len(extents1) == len(extents2) and all(
        (c1, s1) == (c2, s2) and b1 - LOCUS_PADDING < e2 + LOCUS_PADDING and b2 - LOCUS_PADDING < e1 + LOCUS_PADDING
        for (c1, s1, b1, e1), (c2, s2, b2, e2) in zip(extents1, extents2))


def _merge_same_loci(groups, read_to_pieces):
    "groups of reads, those whose loci are in the same places (_same_loci) merged, transitively"
    extents = [_locus_extents(reads, read_to_pieces) for reads in groups]
    parent = list(range(len(groups)))

    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    # only groups with their loci on the same chromosomes and strands can have the same loci
    by_place = defaultdict(list)
    for i, ext in enumerate(extents):
        by_place[tuple((c, s) for c, s, _, _ in ext)].append(i)
    for members in by_place.values():
        for x, i in enumerate(members):
            for j in members[x + 1:]:
                if _same_loci(extents[i], extents[j]):
                    parent[root(j)] = root(i)
    merged = defaultdict(list)
    for i, reads in enumerate(groups):
        merged[root(i)].extend(reads)
    return list(merged.values())


def reads_by_genes(read_to_pieces, gene_to_all_exons, max_loci):
    """The genes of a fusion's loci in read order -> its reads.  Reads are grouped by
    the genes of their pieces, then groups whose loci are in the same places at every
    position (_same_loci), on one strand and overlapping, are merged: reads of one
    fusion whose pieces were given different genes, by ending in a gene's exons or
    short of them, were otherwise split into several fusions, each a synthetic contig
    of the same sequence its reads were divided between at random.  A merged fusion's
    genes are its reads' consensus (consensus_gene).  Reads of 2 to max_loci genes, in
    fusions with an annotated 5' gene; overlapping non-genic loci are named as one."""
    nongenic_loci = sorted(piece[2] for pieces in read_to_pieces.values() for piece in pieces if type(piece[2]) is tuple)
    nongenic_to_locus = group_nongenic_loci(nongenic_loci)
    by_genes = defaultdict(list)
    for read in read_to_pieces:
        read_genes = tuple(nongenic_to_locus.get(piece[2], piece[2]) for piece in sorted(read_to_pieces[read]))
        if 2 <= len(set(read_genes)) <= max_loci:
            by_genes[read_genes].append(read)
    genes_to_reads = defaultdict(list)
    for reads in _merge_same_loci(list(by_genes.values()), read_to_pieces):
        pieces = [sorted(read_to_pieces[r]) for r in reads]
        genes = tuple(consensus_gene([nongenic_to_locus.get(p[i][2], p[i][2]) for p in pieces])
                      for i in range(len(pieces[0])))
        if genes[0] in gene_to_all_exons:  # 5' gene is annotated
            genes_to_reads[genes].extend(reads)
    return dict(genes_to_reads)


def fusion_partners(fusions, gene_to_paralogs):
    "(gene, role in the fusion) -> the sets of its partners, paralogs as their group"
    gene_to_partners = defaultdict(set)
    for fusion_genes in fusions:
        for pos, gene in enumerate(fusion_genes):
            others = set(fusion_genes) - {gene}
            other_groups = frozenset(gene_to_paralogs.get(g, g) for g in others)
            gene_to_partners[(gene, fusion_role(pos, len(fusion_genes)))].add(other_groups)
    return gene_to_partners


def fusion_read_positions(reads, read_to_pieces, num_loci):
    """of reads with a piece per locus: the reads, the read gaps between consecutive pieces
    per junction, and per locus the reads' positions at the pieces' 5' and 3' ends"""
    read_gaps = [[] for x in range(num_loci - 1)]
    piece_ends = [[[], []] for x in range(num_loci)]
    good_reads = []
    for read in reads:
        if len(read_to_pieces[read]) != num_loci:
            continue
        pieces = sorted(read_to_pieces[read])
        good_reads.append(read)
        for i in range(num_loci - 1):
            # positive where the pieces overlap in the read, negative where bases between them are unaligned
            read_gaps[i].append(pieces[i][1][0] - pieces[i + 1][0][0])
        for i in range(num_loci):
            piece_ends[i][0].append(pieces[i][0][1])
            piece_ends[i][1].append(pieces[i][1][1])
    return good_reads, read_gaps, piece_ends


def outer_bound(positions, low):
    """the outermost of positions held by 2 or more reads within 300 bases of each other,
    or the median if further out or none are, padded by LOCUS_PADDING: the lowest if low, else the
    highest"""
    poslist = sorted(positions)
    simplemed = median(poslist)
    groups, g = [], [-500]
    for p in poslist:
        if p - g[-1] > 300:
            if g[0] != -500:
                groups.append(g)
            g = [p]
        else:
            g.append(p)
    groups.append(g)
    goodpos = [p for g in groups if len(g) > 1 for p in g]
    if low:
        outpos = int(min(simplemed, min(goodpos))) - LOCUS_PADDING if len(goodpos) > 0 else int(simplemed) - LOCUS_PADDING
    else:
        outpos = int(max(simplemed, max(goodpos))) + LOCUS_PADDING if len(goodpos) > 0 else int(simplemed) + LOCUS_PADDING
    return max(0, outpos)


def locus_bounds(piece_ends):
    "per locus, its padded outer 5' and 3' bounds in the direction of its reads"
    bounds = []
    for starts, ends in piece_ends:
        strand = '+' if median(ends) > median(starts) else '-'
        bounds.append([outer_bound(starts, strand == '+'), outer_bound(ends, strand == '-')])
    return bounds


def tss_check(gene_info, bounds5, median5, max_dist_to_tss):
    """Does the 5' gene's locus run in its direction, implying a plausible promoter, with
    the distance to its TSS in max_dist_to_tss.  Returns that and the distance, -1 when no
    start is upstream, with the check off."""
    gene_strand, gene_tss = gene_info[3], gene_info[-1]
    # the most upstream 5' end, held by 2 or more reads, which with 5' truncated
    # reads is the one near the promoter, or the median one: either counts
    end5 = bounds5[0] - LOCUS_PADDING if bounds5[0] > bounds5[1] else bounds5[0] + LOCUS_PADDING
    tss_dist = fusion_tss_distance((end5, median5), gene_strand, gene_tss)
    forward = (gene_strand == '+' and bounds5[0] < bounds5[1]) or (gene_strand == '-' and bounds5[0] > bounds5[1])
    if not forward or (max_dist_to_tss is not None and (tss_dist is None or tss_dist > max_dist_to_tss)):
        return False, None
    return True, -1 if tss_dist is None else tss_dist


def breakpoints_consistent(read_gaps):
    "do the reads' pieces have few gaps between them, or very consistent ones"
    spacings = []
    for gaps in read_gaps:
        gaps = sorted(gaps)
        spacing = [gaps[i] - gaps[i - 1] for i in range(1, len(gaps))]
        # a single supporting read gives no consecutive differences,
        # and median([]) raises.  One read agrees with itself
        spacings.append(median(spacing) if len(spacing) > 0 else 0)
    return max([abs(median(x)) for x in read_gaps]) <= 10 \
        or (max([abs(min(x)) for x in read_gaps]) <= 10 and max(spacings) <= 3)


def call_fusion(fusion_genes, reads, read_to_pieces, genetoinfo, chroms, max_dist_to_tss):
    "the fusion's reads, TSS distance, read gaps and loci, or None if it fails a check"
    num_loci = len(fusion_genes)
    good_reads, read_gaps, piece_ends = fusion_read_positions(reads, read_to_pieces, num_loci)
    # the typical 5' end of the reads, for the TSS distance with the most upstream one.
    # A few reads starting further upstream than the rest set that one: in MCF7,
    # ESR1--CCDC170's reads mostly start at an annotated ESR1 TSS, but 3 start 14 kb
    # upstream, 19 kb from the nearest annotated start upstream of them
    median5 = median(piece_ends[0][0])
    bounds = locus_bounds(piece_ends)
    if fusion_genes[0] not in genetoinfo:
        return None
    passes, tss_dist = tss_check(genetoinfo[fusion_genes[0]], bounds[0], median5, max_dist_to_tss)
    if not passes or not breakpoints_consistent(read_gaps):
        return None
    info = {'reads': set(good_reads), 'disttostart': [tss_dist], 'read_gaps': read_gaps}
    for i in range(num_loci):
        info[fusion_genes[i]] = [chroms[i], bounds[i][0], bounds[i][1]]
    return info


def _num_genes(pieces):
    "the number of distinct genes or non-genic loci of a read's pieces"
    return len({piece[2] for piece in pieces})


def _query_span(piece):
    return min(piece[0][0], piece[1][0]), max(piece[0][0], piece[1][0])


def _singly_aligned(pieces):
    "the number of read bases in exactly one of a read's pieces"
    edges = sorted(edge for piece in pieces for edge in ((_query_span(piece)[0], 1), (_query_span(piece)[1], -1)))
    total, depth, prev = 0, 0, 0
    for pos, step in edges:
        if depth == 1:
            total += pos - prev
        depth += step
        prev = pos
    return total


def _use_transcript_pieces(transcript_pieces, genome_pieces):
    """a read's transcript pieces over its genome ones if they give it more genes and
    explain as much of the read, by bases aligned once, not left unaligned or aligned
    twice: a 3-gene fusion's short middle exon, which the genome alignment missed or
    aligned as part of the piece before it, rather than a short alignment to another
    gene of sequence the genome alignment put in one piece"""
    return _num_genes(transcript_pieces) > _num_genes(genome_pieces) and \
        _singly_aligned(transcript_pieces) >= _singly_aligned(genome_pieces)


def id_chimeras(genome_bam, annot, min_support, max_loci=10, max_dist_to_tss=None, max_partners=4,
                transcriptome_bam=None, min_mapq=0, transcriptome_reads=None):
    """The preliminary fusions of the chimeric reads of the genome alignment and, if
    given, of the alignment to annotated transcripts (of transcriptome_reads), as one
    set of reads: each read's genome pieces, or its transcript ones if they give it
    more genes and explain as much of the read, as a 3-gene fusion's short middle exon
    the genome alignment missed (_use_transcript_pieces).
    The transcript alignment can split a read into more pieces of one gene, aligned to
    different transcripts, so pieces alone don't decide.  Reads are grouped into
    fusions by reads_by_genes, and min_support and the filters apply to each fusion's
    reads from both alignments together.  annot is flair fusion's Annotation.
    Returns the fusions' genes joined by __ -> their reads, TSS distance, read gaps and
    loci (call_fusion)."""
    sharing = exon_sharing_genes(annot.chrom_to_gene_pos, annot.gene_to_all_exons)

    # each read's pieces, from the alignment to transcripts, then the genome alignment
    read_to_pieces = {}
    if transcriptome_bam is not None:
        piece_of = partial(transcriptomic_piece, genetoinfo=annot.genetoinfo, intronLocs=annot.intronLocs,
                           intronToGenome=annot.intronToGenome, sharing=sharing,
                           gene_to_all_exons=annot.gene_to_all_exons)
        read_to_pieces = read_pieces('transcriptomic', transcriptome_bam, piece_of, min_mapq, transcriptome_reads)
        drop_foldbacks('transcriptomic', read_to_pieces)
    piece_of = partial(genomic_piece, chrom_to_gene_pos=annot.chrom_to_gene_pos,
                       gene_to_all_exons=annot.gene_to_all_exons, juncs_to_gene=annot.juncs_to_gene,
                       sharing=sharing, genetoinfo=annot.genetoinfo)
    read_to_genome_pieces = read_pieces('genomic', genome_bam, piece_of, min_mapq, None)
    drop_foldbacks('genomic', read_to_genome_pieces)

    # the genome pieces win unless the transcript ones find a gene they missed
    for read, genome_pieces in read_to_genome_pieces.items():
        if read not in read_to_pieces or not _use_transcript_pieces(read_to_pieces[read], genome_pieces):
            read_to_pieces[read] = genome_pieces
    num_from_genome = sum(read_to_pieces[read] is genome_pieces for read, genome_pieces in read_to_genome_pieces.items())
    logging.info(f'chimeric reads: {num_from_genome} with pieces from the genome alignment, '
                 f'{len(read_to_pieces) - num_from_genome} from the alignment to transcripts')

    # reads grouped into fusions by their genes, each fusion then filtered as a whole
    genes_to_reads = reads_by_genes(read_to_pieces, annot.gene_to_all_exons, max_loci)
    supported = [genes for genes, reads in genes_to_reads.items() if len(reads) >= min_support]
    gene_to_partners = fusion_partners(supported, annot.gene_to_paralogs)
    fusion_to_info = {}
    for fusion_genes in supported:
        chroms = [annot.genetoinfo[g][0] if g in annot.genetoinfo else g.split(':')[0] for g in fusion_genes]
        if not identify_fusion_problems(fusion_genes, gene_to_partners, max_partners, annot.genetoname,
                                        annot.gene_to_paralogs, chroms, len(genes_to_reads[fusion_genes])):
            continue
        info = call_fusion(fusion_genes, genes_to_reads[fusion_genes], read_to_pieces, annot.genetoinfo, chroms,
                           max_dist_to_tss)
        if info is not None:
            fusion_to_info['__'.join(fusion_genes)] = info
    return fusion_to_info
