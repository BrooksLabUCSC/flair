import logging
import re
from flair.io_utils import open_sorted_bed
from flair.pycbio.hgdata.bed import BedBlock, BedReader
from flair.pycbio.sys import fileOps
from flair.flair_bed import FlairBed
from flair.isoform_data import Exon, Junc
from flair.gene_identification import (get_genes_with_shared_juncs, genes_sharing_splice_sites,
                                       get_spliced_exon_overlaps, get_unspliced_exon_overlaps)

def get_synth_info(breakpointfile):
    synthchrtoinfo = {}
    for bed in BedReader(breakpointfile, numStdCols=4):
        synthchrtoinfo[bed.chrom] = '--'.join(bed.name.split('--')[1:])
    return synthchrtoinfo

def get_paralog_ref(paralogfile):
    grouptogenes = {}
    genetoparalogs = {}
    # 5  244     2       69321074        69338940        1       ENSG00000205572 SERF1B  small EDRK-rich factor 1B (centromeric) [Source:HGNC Symbol;Acc:10756]
    for line in fileOps.iterLines(paralogfile):
        line = line.split('\t')
        group = line[1]
        gname = line[6]
        if group not in grouptogenes:
            grouptogenes[group] = set()
        grouptogenes[group].add(gname)

    for group in grouptogenes:
        gset = grouptogenes[group]
        for g in gset:
            genetoparalogs[g] = group
    return genetoparalogs

def paralog_groups_for_genes(genetoparalogs, gene_ids):
    """Key the paralog groups by the annotation's gene ids.  The paralog table
    has unversioned Ensembl ids, while GENCODE ids carry a version
    (ENSG00000215915.10), so looking them up as written never matched and no
    fusion was rejected as between paralogs."""
    return {g: genetoparalogs[g.split('.')[0]] for g in gene_ids
            if g.split('.')[0] in genetoparalogs}


# mitochondrial sequence names in UCSC, GENCODE, Ensembl and RefSeq genomes
MITOCHONDRIA = {'chrM', 'chrMT', 'M', 'MT', 'NC_012920.1'}


def is_reference_chromosome(chrom):
    """Is chrom a nuclear chromosome, named as UCSC and GENCODE (chr12), RefSeq
    (NC_000012.12) or Ensembl (12, X) do: not mitochondria, nor an unplaced scaffold
    (GL000009.2, KI270706.1, RefSeq's NT_ and NW_)"""
    if chrom in MITOCHONDRIA:
        return False
    return chrom.startswith('chr') or chrom.startswith('NC_') or chrom.isdigit() or chrom in ('X', 'Y')


NONGENIC_LOCUS = re.compile(r'.+:\d+-\d+$')


def is_nongenic_locus(name):
    "is name a non-genic locus, chrom:start-end, rather than a gene id"
    return NONGENIC_LOCUS.match(name) is not None


def fusion_role(i, numgenes):
    "a gene's role in a fusion by its position: 5' end, 3' end or between them"
    return '5p' if i == 0 else '3p' if i == numgenes - 1 else 'mid'


def identify_fusion_problems(fgenes, locustopartners, maxpromiscuity, genetoname, genetoparalogs, genomic_chroms, isosup):
    """locustopartners is (gene, fusion_role) -> its partner sets: a gene is promiscuous
    with more than maxpromiscuity partners in the same role.  Counted over all roles, a
    gene of an amplicon was rejected as the 3' partner of a well supported fusion for
    partners it had as a 5' one: SULF2 in MCF7, 3' of ARFGEF2 (492 reads) and 5' of
    PRICKLE2, FTH1, GNAS, KRT8 and an intergenic locus."""
    ispromiscuous, areparalogs, areig, allnames, gcount = False, [], [], set(), []
    for i, g in enumerate(fgenes):
        gcount.append(list(fgenes).count(g))
        if len(locustopartners[(g, fusion_role(i, len(fgenes)))]) > maxpromiscuity:
            ispromiscuous = True
        if g in genetoname:
            # the whole name: genes of one name, as a gene's X and Y copies, are rejected,
            # paralogs by the paralog table.  The first 4 letters used to be compared,
            # which rejected unrelated genes sharing a stem, as BCAS4--BCAS3 in MCF7
            allnames.add(genetoname[g])
            if genetoname[g][:3] in {'IGK', 'IGH', 'IGL', 'IGF'}:
                areig.append(True)
            else:
                areig.append(False)
        else:
            allnames.add(g)
            # if g[:3] == 'chr':
            #     areig.append(True)
            # else:
            areig.append(False)
        if g in genetoparalogs:
            other = set(fgenes) - {g, }
            other = {genetoparalogs[g2] if g2 in genetoparalogs else g2 for g2 in other}
            if genetoparalogs[g] in other or len(other) == 0:
                areparalogs.append(True)
            else:
                areparalogs.append(False)
        else:
            areparalogs.append(False)
    overall = (all([x == 1 for x in gcount]) and isosup >= 1 and not ispromiscuous
               and any(x is False for x in areparalogs) and all(x is False for x in areig)
               and len(allnames) > 1 and all(is_reference_chromosome(x) for x in genomic_chroms))
    return overall

def get_locus_bounds(synthinfo):
    locuslen = [abs(x[3] - x[2]) for x in synthinfo]
    locusbounds = []
    laststart = 0
    for i in range(len(locuslen)):
        locusbounds.append((laststart, locuslen[i] + laststart))
        laststart += locuslen[i]
    return locusbounds

def separate_exons_by_locus(esizes, estarts, numloci, locusbounds, start, thickStart, thickEnd):
    """The first exon start, relative to start, and the exon indexes of each locus, and
    each locus's part of the coding region, thickStart to thickEnd in synthetic
    coordinates: None for a locus the coding region doesn't reach, or when there is
    none (thickStart == thickEnd)"""
    starts = [None for x in range(numloci)]
    exonindexes = [[] for x in range(numloci)]
    thickEdges = [None for x in range(numloci)]
    for order, (locusstart, locusend) in enumerate(locusbounds):
        if thickStart < thickEnd and thickStart < locusend and locusstart < thickEnd:
            thickEdges[order] = (max(thickStart, locusstart), min(thickEnd, locusend))
    for i in range(len(esizes)):
        thisstart, thisend = start + estarts[i], start + estarts[i] + esizes[i]
        for order in range(numloci):
            # <= on both edges: with a strict test on the left an exon starting
            # exactly on a boundary matched no locus, failing its own locus's left
            # test and the previous locus's right test, and locus 0 starts at 0
            if locusbounds[order][0] <= thisstart and thisend <= locusbounds[order][1]:
                if starts[order] is None:
                    starts[order] = estarts[i]  # thisstart #- locusbounds[order][0]
                exonindexes[order].append(i)
    return starts, exonindexes, thickEdges

def convert_to_genomic_coords(numloci, synthinfo, exonindexes, starts, locusbounds, esizes, estarts, start, iso, thickEdges):
    """Each locus's piece of a synthetic isoform in genome coordinates.  A piece's
    coding region is its part of the isoform's (separate_exons_by_locus), within its
    exons' span: a piece the coding region doesn't reach has none (thickStart and
    thickEnd None, written as chromEnd), rather than being all coding"""
    genomicbounds = []
    bedinfo = []
    for order in range(numloci):
        genename, genomicchr, leftbound, rightbound = synthinfo[order]
        locusesizes = [esizes[i] for i in exonindexes[order]]
        locusestarts = [estarts[i] for i in exonindexes[order]]
        locusestarts = [x - starts[order] for x in locusestarts]
        totlen = locusestarts[-1] + locusesizes[-1]
        locusdir = '+' if leftbound < rightbound else '-'
        # name = 'gene' + str(order + 1) + '_' + iso
        thickStart, thickEnd = None, None
        if leftbound > rightbound:  # reverse direction
            temp = []
            for i in range(len(locusestarts) - 1, -1, -1):
                temp.append(totlen - (locusestarts[i] + locusesizes[i]))
            locusestarts = temp
            locusesizes = locusesizes[::-1]
            chromStart = leftbound - (((start + starts[order]) - locusbounds[order][0]) + totlen)
            if thickEdges[order] is not None:
                thickStart = leftbound - (thickEdges[order][1] - locusbounds[order][0])
                thickEnd = leftbound - (thickEdges[order][0] - locusbounds[order][0])
        else:
            chromStart = leftbound + ((start + starts[order]) - locusbounds[order][0])
            if thickEdges[order] is not None:
                thickStart = leftbound + (thickEdges[order][0] - locusbounds[order][0])
                thickEnd = leftbound + (thickEdges[order][1] - locusbounds[order][0])
        chromEnd = chromStart + totlen
        if thickStart is not None:
            thickStart, thickEnd = max(thickStart, chromStart), min(thickEnd, chromEnd)
            if thickStart >= thickEnd:
                thickStart, thickEnd = None, None

        genomicbounds.append((genomicchr, chromStart, chromEnd, locusdir))
        blocks = [BedBlock(chromStart + locusestarts[i], chromStart + locusestarts[i] + locusesizes[i])
                  for i in range(len(locusesizes))]
        bedinfo.append((genomicchr, chromStart, chromEnd, locusdir, thickStart, thickEnd, blocks))
    return genomicbounds, bedinfo

def check_too_close(numloci, genomicbounds, min_dist):
    for i in range(numloci - 1):
        if genomicbounds[i][0] == genomicbounds[i + 1][0] and genomicbounds[i][3] == genomicbounds[i + 1][3]:  # chrom and strand are same
            if genomicbounds[i][3] == '+' and 0 < genomicbounds[i + 1][1] - genomicbounds[i][2] < min_dist:
                return True
            # bounds are start < end on either strand: on -, the gap runs from the
            # downstream piece's end up to the upstream piece's start
            elif genomicbounds[i][3] == '-' and 0 < genomicbounds[i][1] - genomicbounds[i + 1][2] < min_dist:
                return True
    return False

def get_gene_id(gene_og_to_new, gene_og, curr_gene_count):
    if gene_og not in gene_og_to_new:
        gene_id = f'FLG{curr_gene_count:08d}'
        gene_info = tuple(gene_og.split('--'))
        gene_og_to_new[gene_og] = (gene_id, gene_info)
        curr_gene_count += 1
    gene_id, gene_info = gene_og_to_new[gene_og]
    return curr_gene_count, gene_id, gene_info

def breakpoint_sites(blocks, strand, order, numloci):
    """A fusion piece's splice sites at its fusion junctions: its 3' edge, unless it
    is the last piece, and its 5' edge, unless it is the first.  They are set by
    splicing, as its junctions are, while its other ends are where reads ended."""
    sites = []
    if order < numloci - 1:
        sites.append(blocks[-1].end if strand == '+' else blocks[0].start)
    if order > 0:
        sites.append(blocks[0].start if strand == '+' else blocks[-1].end)
    return sites


def piece_gene(chrom, strand, blocks, order, numloci, annots_by_chrom):
    """The gene of a fusion piece in genome coordinates, as flair transcriptome gives
    an isoform its gene (get_genes_with_shared_juncs), by its junctions and its
    splice sites at the fusion junctions (breakpoint_sites), not its ends: a piece
    of only a first exon has its donor.  Of several genes the piece reads through,
    the one with a splice site at a fusion junction, the gene fused, then the one
    with more annotated transcripts.  A piece with no annotated junction or splice
    site gets the gene whose exons cover most of its own, more than half of them (as
    get_spliced_exon_overlaps, then get_unspliced_exon_overlaps), the one with more
    annotated transcripts on a tie.  A piece in no gene's exons is None, a non-genic
    locus name_nongenic_loci names.  The gene depends only on the piece and the
    annotation, so flair combine gives a combined fusion isoform its genes the same
    way, whatever the samples' preliminary fusions."""
    annots = annots_by_chrom.get(chrom)
    if annots is None:  # no genes on the chromosome
        return None
    juncs = tuple(Junc(blocks[i].end, blocks[i + 1].start) for i in range(len(blocks) - 1))
    sites = breakpoint_sites(blocks, strand, order, numloci)
    genes = get_genes_with_shared_juncs(juncs, strand, annots, extra_sites=sites)
    if genes:
        fused = genes_sharing_splice_sites(sites, strand, annots)
        return min(genes, key=lambda g: (g not in fused, -annots.gene_to_num_transcripts.get(g, 0), g))
    exons = [Exon(blk.start, blk.end) for blk in blocks]
    hits = get_spliced_exon_overlaps(strand, exons, annots) or get_unspliced_exon_overlaps(strand, exons, annots)
    if hits:
        return min(hits, key=lambda hit: (-hit[0], -annots.gene_to_num_transcripts.get(hit[1], 0), hit[1]))[1]
    return None


def name_nongenic_loci(converted):
    """Name the pieces piece_gene left in no gene (None) by their synthetic contig's
    locus: chrom:start-end spanning all such pieces of it, so isoforms of one fusion
    with different ends there have one non-genic locus.  converted is a list of
    (bed, bedinfos, genes), and genes are changed in place."""
    spans = {}
    for bed, bedinfos, genes in converted:
        for order, gene in enumerate(genes):
            if gene is None:
                chrom, start, end = bedinfos[order][:3]
                key = (bed.chrom, order)
                prev = spans.get(key, (chrom, start, end))
                spans[key] = (chrom, min(start, prev[1]), max(end, prev[2]))
    for bed, bedinfos, genes in converted:
        for order, gene in enumerate(genes):
            if gene is None:
                genes[order] = '{}:{}-{}'.format(*spans[(bed.chrom, order)])


def convert_synthetic_isos(isoformsbed, breakpointfile, outname, min_dist, annots_by_chrom):
    """The fusion isoforms of the synthetic genome's flair transcriptome, in genome
    coordinates, a bed line for each piece.  Each piece's gene is from its
    junctions and fusion-junction splice sites (piece_gene), with annots_by_chrom,
    an AnnotsByChrom, or its non-genic locus (name_nongenic_loci).  An isoform is left out
    if it doesn't cross a breakpoint, its pieces are too close (check_too_close) or
    two pieces are of one gene; every other output of the fusion isoforms is
    filtered to those of this bed.  Returns isoform -> its pieces' genes, in fusion
    order, non-genic loci as chrom:start-end (name_nongenic_loci)"""
    synthchrtoinfo = get_synth_info(breakpointfile)

    converted = []
    for bed in BedReader(isoformsbed, fixScores=True, bedClass=FlairBed):
        iso, start = bed.name, bed.chromStart
        esizes = [len(blk) for blk in bed.blocks]
        estarts = [blk.start - start for blk in bed.blocks]
        synthinfo = [x.split('..') for x in synthchrtoinfo[bed.chrom].split('--')]
        synthinfo = [[y[0], y[1], int(y[2]), int(y[3])] for y in synthinfo]

        locusbounds = get_locus_bounds(synthinfo)
        if start < locusbounds[0][1] and locusbounds[-1][0] < bed.chromEnd:
            numloci = len(synthinfo)

            starts, exonindexes, thickEdges = separate_exons_by_locus(esizes, estarts, numloci, locusbounds, start, bed.thickStart, bed.thickEnd)

            assigned_exons = sorted(i for indexes in exonindexes for i in indexes)
            if assigned_exons != list(range(len(esizes))):
                logging.warning('skipping %s: not all synthetic exons could be assigned to exactly one locus', iso)
                continue

            if None not in starts:
                genomicbounds, bedinfos = convert_to_genomic_coords(numloci, synthinfo, exonindexes, starts, locusbounds, esizes, estarts, start, iso, thickEdges)
                if not check_too_close(numloci, genomicbounds, min_dist):
                    genes = [piece_gene(genomicchr, locusdir, blocks, order, numloci, annots_by_chrom)
                             for order, (genomicchr, _, _, locusdir, _, _, blocks) in enumerate(bedinfos)]
                    converted.append((bed, bedinfos, genes))

    name_nongenic_loci(converted)
    out = open_sorted_bed(outname)
    gene_og_to_new = {}
    curr_gene_count = 1
    iso_genes = {}
    for bed, bedinfos, genes in converted:
        if len(set(genes)) < len(genes):
            # pieces of one gene: a rearrangement within a gene, as of a
            # segmental duplication, which flair fusion doesn't call
            logging.info('skipping %s: pieces of one gene, %s', bed.name, '--'.join(genes))
            continue
        iso_genes[bed.name] = tuple(genes)
        for fusion_order, (genomicchr, chromStart, chromEnd, locusdir, thickStart, thickEnd, blocks) in enumerate(bedinfos):
            curr_gene_count, fusion_id, ref_gene_info = get_gene_id(gene_og_to_new, '--'.join(genes), curr_gene_count)
            flair_genes = []
            for gene in ref_gene_info:
                curr_gene_count, gene_id, ref_gene_ids = get_gene_id(gene_og_to_new, gene, curr_gene_count)
                flair_genes.append(gene_id)
            # convert reference genes that are chromosome regions to none in the ref_gene_mappings field
            ref_gene_info = tuple([None if is_nongenic_locus(x) else x for x in ref_gene_info])
            outbed = FlairBed(genomicchr, chromStart, chromEnd, name=bed.name, score=bed.score, strand=locusdir, blocks=blocks,
                              gene_id=fusion_id, ref_gene_mappings=ref_gene_info, thickStart=thickStart, thickEnd=thickEnd,
                              read_support=bed.read_support, frac_support=bed.frac_support, productivity=bed.productivity, transcript_class='fusion',
                              fused_genes=tuple(flair_genes), pos_in_fusion=fusion_order, samples=bed.samples,
                              aaseq_id=bed.aaseq_id)
            outbed.write(out)
    out.close()
    return iso_genes
