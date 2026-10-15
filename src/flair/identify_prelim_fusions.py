#!/usr/bin/env python3

import logging
import pysam
from collections import Counter, defaultdict
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


def getCorrectGene(chrom_to_gene_pos, gene_to_all_exons, juncs_to_gene, chrom, readblocks, thisdir):  # noqa: C901 - FIXME: reduce complexity
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


def is_foldback(aligns):
    """Is a chimeric read a fold-back: two of its pieces on opposite strands of one
    chromosome overlapping or within FOLDBACK_MAX_DIST.  cDNA that folds back on
    itself gives a read of a region joined to its own reverse complement, and the
    piece on the other strand would be assigned to an antisense gene there, as
    SLC7A5--SLC7A5-AS1 in MCF7.  aligns are id_chimeras' pieces, [(qstart, refpos),
    (qend, refpos), gene, strand, chrom]."""
    spans = [(a[4], a[3], min(a[0][1], a[1][1]), max(a[0][1], a[1][1])) for a in aligns]
    for i, (chrom, strand, start, end) in enumerate(spans):
        for chrom2, strand2, start2, end2 in spans[i + 1:]:
            if chrom2 == chrom and strand2 != strand and max(start, start2) - min(end, end2) <= FOLDBACK_MAX_DIST:
                return True
    return False


def is_low_quality(read, min_mapq):
    """Is this piece of a chimeric read below min_mapq: a piece that aligns about as well
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
    with pysam.AlignmentFile(genome_bam, 'rb') as aligns:
        chim_out = pysam.AlignmentFile(chim_bam, 'wb', template=aligns) if chim_bam else None
        fa_out = open(reads_fa, 'w') if reads_fa else None
        for read in aligns.fetch(until_eof=True):
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


def id_chimeras(mode, bam, genetoinfo, chrom_to_gene_pos, gene_to_all_exons, juncs_to_gene, gene_to_paralogs,  # noqa: C901 - FIXME: reduce complexity
                genetoname, minsup, maxloci=10, reqdisttostart=None, maxpromiscuity=4, intronLocs=None, intronToGenome=None,
                min_mapq=0, only_reads=None):
    isrevtosign = {True: '-', False: '+'}
    # mapping quality is only used for genome alignments: aligned to transcripts, a
    # read of a gene with several overlapping isoforms has a low one regardless.  Reads
    # with a low piece are left out after the pass, as the piece can come after the others
    low_quality = set()
    withsup = pysam.AlignmentFile(bam, "rb")
    readToAligns = {}
    sharing = exon_sharing_genes(chrom_to_gene_pos, gene_to_all_exons)
    for read in withsup:
        if mode == 'genomic' and is_low_quality(read, min_mapq):
            low_quality.add(read.query_name)
        if read.is_mapped and not read.is_secondary and read.has_tag('SA') and (mode == 'genomic' or not read.is_reverse) \
                and (only_reads is None or read.query_name in only_reads):  # if aligned to transcriptome, must match strand of transcript
            rname = read.query_name
            if rname not in readToAligns:
                readToAligns[rname] = []

            refstart, refend, readdir = read.reference_start, read.reference_end, isrevtosign[read.is_reverse]
            qstart, qend = read.query_alignment_start, read.query_alignment_end
            readlen = read.infer_read_length()
            cigar = read.cigartuples
            if cigar[0][0] == 5:  # just hard clipping
                qstart += cigar[0][1]
                qend += cigar[0][1]

            if mode == 'genomic':
                refchr = read.reference_name
                readblocks = get_exon_intron_blocks(read)
                genename = getCorrectGene(chrom_to_gene_pos, gene_to_all_exons, juncs_to_gene, refchr, readblocks, readdir)
                if type(genename) is not tuple:  # tuples are non-genic loci
                    genename = prefer_shorter_gene(genename, readblocks[1], sharing, genetoinfo, gene_to_all_exons)
                if readdir == '+':
                    readToAligns[rname].append([(qstart, refstart), (qend, refend), genename, readdir, refchr])
                else:
                    readToAligns[rname].append([(readlen - qend, refend), (readlen - qstart, refstart), genename, readdir, refchr])
            else:
                tname, genename = read.reference_name.split('|')
                refchr, genedir = genetoinfo[genename][0], genetoinfo[genename][3]  # can do this because already required the read to be forward strand
                refstart = getGenomicPreciseLoc(tname, refstart, genedir, intronLocs, intronToGenome)
                refend = getGenomicPreciseLoc(tname, refend, genedir, intronLocs, intronToGenome)
                refstart, refend = min(refstart, refend), max(refstart, refend)
                genename = prefer_shorter_gene(genename, [(refstart, refend)], sharing, genetoinfo, gene_to_all_exons)
                if genedir == '+':
                    readToAligns[rname].append([(qstart, refstart), (qend, refend), genename, readdir, refchr])
                else:
                    readToAligns[rname].append([(qstart, refend), (qend, refstart), genename, genedir, refchr])
    withsup.close()
    for r in low_quality:
        readToAligns.pop(r, None)
    if low_quality:
        logging.info(f'{mode} chimeras: {len(low_quality)} reads with a piece below MAPQ {min_mapq} left out')

    foldbacks = [r for r, aligns in readToAligns.items() if is_foldback(aligns)]
    for r in foldbacks:
        del readToAligns[r]
    logging.info(f'{mode} chimeras: {len(foldbacks)} fold-back reads of {len(readToAligns) + len(foldbacks)} left out')
    nongenicloci = [a[2] for aligns in readToAligns.values() for a in aligns if type(a[2]) is tuple]
    nongenicloci.sort()
    ognongenic_to_genes = group_nongenic_loci(nongenicloci)

    interestingloci = {}
    for read in readToAligns:
        alignedloci = sorted(readToAligns[read])

        readgenes = [x[2] for x in alignedloci]
        if 2 <= len(set(readgenes)) <= maxloci:
            if readgenes[0] in gene_to_all_exons:  # 5' gene is annotated
                readgenes = [ognongenic_to_genes[x] if x in ognongenic_to_genes else x for x in readgenes]
                info = tuple(readgenes)
                if info not in interestingloci:
                    interestingloci[info] = []
                interestingloci[info].append(read)

    goodcov = []
    for l in interestingloci:
        if len(interestingloci[l]) >= minsup:
            goodcov.append(l)

    locustopartners = defaultdict(set)  # (gene, role in the fusion) -> partner sets
    for fgenes in goodcov:
        for pos, i in enumerate(fgenes):
            other = set(fgenes) - {i, }
            newother = frozenset([gene_to_paralogs[g] if g in gene_to_paralogs else g for g in other])
            locustopartners[(i, fusion_role(pos, len(fgenes)))].add(newother)

    fusiontoinfo = {}
    for fgenes in goodcov:
        genomic_chroms = [genetoinfo[g][0] if g in genetoinfo else g.split(':')[0] for g in fgenes]
        is_good_fusion = identify_fusion_problems(fgenes, locustopartners, maxpromiscuity, genetoname, gene_to_paralogs, genomic_chroms, len(interestingloci[fgenes]))
        if is_good_fusion:
            numloci = len(fgenes)
            qdist = [[] for x in range(numloci - 1)]
            alignblocks = [[[], []] for x in range(numloci)]
            goodreads = []
            for r in interestingloci[fgenes]:
                if len(readToAligns[r]) == numloci:
                    alignedloci = sorted(readToAligns[r])
                    goodreads.append(r)
                    for i in range(numloci - 1):
                        qdist[i].append(alignedloci[i][1][0] - alignedloci[i + 1][0][0])
                    for i in range(numloci):
                        alignblocks[i][0].append(alignedloci[i][0][1])
                        alignblocks[i][1].append(alignedloci[i][1][1])

            # the typical 5' end of the reads, for the TSS distance with the most upstream one
            # below.  A few reads starting further upstream than the rest set that one: in
            # MCF7, ESR1--CCDC170's reads mostly start at an annotated ESR1 TSS, but 3 start
            # 14 kb upstream, 19 kb from the nearest annotated start upstream of them
            median5 = median(alignblocks[0][0])

            # check that gene orders for all reads are consistent
            # to start, no clustering, take simple min/max
            for i in range(numloci):
                strand = '+' if median(alignblocks[i][1]) > median(alignblocks[i][0]) else '-'
                for j in range(2):
                    poslist = sorted(alignblocks[i][j])
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
                    goodpos = []
                    for g in groups:
                        if len(g) > 1:
                            goodpos.extend(g)
                    if (strand == '+' and j == 0) or (strand == '-' and j == 1):
                        outpos = int(min(simplemed, min(goodpos))) - 1000 if len(goodpos) > 0 else int(simplemed) - 1000
                    else:
                        outpos = int(max(simplemed, max(goodpos))) + 1000 if len(goodpos) > 0 else int(simplemed) + 1000
                    outpos = max(0, outpos)
                    alignblocks[i][j] = outpos

            # check that 5' gene is in the forward direction, implies plausible promoter
            if fgenes[0] in genetoinfo:
                firstgenedir = genetoinfo[fgenes[0]][3]
                firstgenetstarts = genetoinfo[fgenes[0]][-1]
                # the most upstream 5' end, held by 2 or more reads, which with 5' truncated
                # reads is the one near the promoter, or the median one: either counts
                end5 = alignblocks[0][0]
                if alignblocks[0][0] > alignblocks[0][1]:
                    end5 -= 1000
                else:
                    end5 += 1000
                mindisttostart = fusion_tss_distance((end5, median5), firstgenedir, firstgenetstarts)

                if ((firstgenedir == '+' and alignblocks[0][0] < alignblocks[0][1])
                        or (firstgenedir == '-' and alignblocks[0][0] > alignblocks[0][1]))\
                        and (reqdisttostart is None or (mindisttostart is not None and mindisttostart <= reqdisttostart)):
                    if mindisttostart is None:  # no upstream start, with the check off
                        mindisttostart = -1
                    simscores = []
                    for qdistlist in qdist:
                        simscore = []
                        qdistlist = sorted(qdistlist)
                        for i in range(1, len(qdistlist)):
                            simscore.append(qdistlist[i] - qdistlist[i - 1])
                        # a single supporting read gives no consecutive differences,
                        # and median([]) raises.  One read agrees with itself
                        simscores.append(median(simscore) if len(simscore) > 0 else 0)
                    if max([abs(median(x)) for x in qdist]) <= 10 \
                            or (max([abs(min(x)) for x in qdist]) <= 10 and max(simscores) <= 3):  # alignments have to either have few gaps or be very consistent
                        fname = '__'.join(fgenes)
                        fusiontoinfo[fname] = {'reads': set(goodreads), 'disttostart': [mindisttostart], 'qdist': qdist}
                        for i in range(numloci):
                            fusiontoinfo[fname][fgenes[i]] = [genomic_chroms[i], alignblocks[i][0], alignblocks[i][1]]

    return fusiontoinfo
