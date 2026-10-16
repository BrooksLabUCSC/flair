"""Gene identification: the annotated genes of an isoform, by the junctions and splice
sites it shares with them, or else by exon overlap, the genes of flair transcriptome's
first-pass isoforms, and the grouping of isoforms of no annotated gene into novel
genes by shared splicing.  Used by flair transcriptome, and by flair combine and flair
fusion to name genes as transcriptome does."""
import hashlib
import logging
from flair.isoform_data import Exon, Gene, binary_search
from flair.junction_correct import UNKNOWN_STRAND

# overlap fraction thresholds for gene assignment
MIN_ISOFORM_OVERLAP_FRAC = 0.5
MIN_ANNOT_OVERLAP_FRAC = 0

# search window for binary search of single-exon annotations
ANNOT_SE_SEARCH_WINDOW = 2
# in gene assignment, a gene's matched splice site this close to the isoform's
# 5' splice site counts as matching the isoform's 5' end
FIVE_PRIME_SS_WINDOW = 20
# in the gene assignment fallback with fixed ends, a spliced isoform's terminal exons
# are taken to be this long, instead of reaching its read ends
FALLBACK_TERMINAL_EXON_LEN = 100


def _gene_span(gene_id, annots):
    """span of the gene's annotated junctions, so a long first or last exon does
    not make a gene look longer or overlap more genes"""
    juncs = annots.gene_to_annot_juncs[gene_id]
    return min(j.start for j in juncs), max(j.end for j in juncs)

def _drop_subset_genes(gene_to_matches):
    """keep the genes whose set of matches to the isoform (junctions or splice
    sites) is not a proper subset of another gene's set"""
    return [g for g, matches in gene_to_matches.items()
            if not any(matches < other for other in gene_to_matches.values())]

def _group_overlapping_genes(genes, annots):
    """group genes whose genomic spans overlap, transitively"""
    groups, group_end = [], None
    for (start, end), gene_id in sorted((_gene_span(g, annots), g) for g in genes):
        if groups and start < group_end:
            groups[-1].append(gene_id)
            group_end = max(group_end, end)
        else:
            groups.append([gene_id])
            group_end = end
    return groups

def _pick_group_gene(group, gene_to_juncs, gene_to_sites, five_prime_ss, strand, annots):
    """Pick one gene from a group of overlapping genes, none of whose matches to
    the isoform contain another's, by the highest score: 2 for each of the
    isoform's junctions the gene has annotated, 1 for each other splice site of
    the isoform it has annotated, and 1 if its matched splice site nearest the
    isoform's 5' end is within FIVE_PRIME_SS_WINDOW of the isoform's 5' splice
    site.  Ties go to the gene with the most annotated transcripts, so a gene and
    a readthrough of it, usually annotated with few transcripts, give the gene,
    then to the gene with the shortest junction span, then the gene id."""
    def score(gene_id):
        juncs = gene_to_juncs.get(gene_id, ())
        junc_sites = {site for j in juncs for site in (j.start, j.end)}
        sites = gene_to_sites.get(gene_id, set())
        matched = junc_sites | sites
        five_prime_match = min(matched) if strand == '+' else max(matched)
        return (2 * len(juncs) + len(sites - junc_sites)
                + (1 if abs(five_prime_match - five_prime_ss) < FIVE_PRIME_SS_WINDOW else 0))

    def length(gene_id):
        start, end = _gene_span(gene_id, annots)
        return end - start
    return min(group, key=lambda g: (-score(g), -annots.gene_to_num_transcripts.get(g, 0), length(g), g))

def _genes_sharing_splice_sites(sites, strand, annots):
    """map of same-strand gene -> set of the isoform's splice sites it has
    annotated"""
    gene_to_sites = {}
    for site in sites:
        for gene_id in annots.splice_site_to_genes.get(site, ()):
            if annots.gene_to_strand[gene_id] == strand:
                gene_to_sites.setdefault(gene_id, set()).add(site)
    return gene_to_sites

def get_genes_with_shared_juncs(juncs, strand, annots, extra_sites=()):
    """Assign a spliced isoform to annotated genes by the junctions and individual
    splice sites it shares with them, with extra_sites, splice sites it has outside
    its junctions, as a fusion piece's at its fusion junctions.  A gene is dropped
    if its matches, junctions and splice sites together, are a proper subset of
    another gene's.  The remaining genes are grouped by genomic overlap of their
    junction spans: each group of overlapping genes, such as a gene cluster sharing
    exons or a readthrough gene and its parts, gives one gene, chosen by
    _pick_group_gene, while genes that do not overlap each other, as for a
    readthrough isoform without a readthrough annotation, are all returned.
    Returns a tuple of gene ids, empty if no gene shares a junction or splice site.
    When some gene shares a junction, genes sharing only splice sites compete only
    within a group with such a gene; they don't make groups of their own, which
    would report a neighboring gene that shares a splice site as if the isoform
    read through it."""
    gene_to_juncs = {}
    for j in juncs:
        if j in annots.junc_to_gene:
            for transcript_id, gene_id in annots.junc_to_gene[j]:
                if gene_id not in gene_to_juncs:
                    gene_to_juncs[gene_id] = set()
                gene_to_juncs[gene_id].add(j)
    sites = {site for j in juncs for site in (j.start, j.end)} | set(extra_sites)
    if not sites:
        return ()
    gene_to_sites = _genes_sharing_splice_sites(sites, strand, annots)
    # a gene with a matching junction has its splice sites, strand or not
    for gene_id, gene_juncs in gene_to_juncs.items():
        gene_to_sites.setdefault(gene_id, set()).update(site for j in gene_juncs for site in (j.start, j.end))
    gene_to_matches = {g: {('junc', j) for j in gene_to_juncs.get(g, ())} | {('site', s) for s in gene_sites}
                       for g, gene_sites in gene_to_sites.items()}
    genes = _drop_subset_genes(gene_to_matches)
    groups = _group_overlapping_genes(genes, annots)
    if gene_to_juncs:
        groups = [group for group in groups if any(g in gene_to_juncs for g in group)]
    # the isoform's first splice site, at the inner edge of its 5' exon for an
    # isoform with junctions
    five_prime_ss = max(sites) if strand == '-' else min(sites)
    return tuple(sorted(_pick_group_gene(group, gene_to_juncs, gene_to_sites, five_prime_ss, strand, annots)
                        for group in groups))


def get_single_exon_gene_overlaps(strand, iso_readrec, annots):
    gene_hits = {}
    exon = iso_readrec.exons[0]
    index = binary_search(exon, annots.all_annot_SE[strand])
    # FIXME: how does this ever work? all_annot_SE is [(start, end, strand, gene_id), ...]
    # max(0, ...): a negative slice start would be read as an offset from the end of
    # the list, which for a long list gives an empty window
    for annot_exon_info in annots.all_annot_SE[strand][max(0, index - ANNOT_SE_SEARCH_WINDOW):index + ANNOT_SE_SEARCH_WINDOW]:
        # FIXME: make overlap a function
        overlap = min(exon.end, annot_exon_info.end) - max(exon.start, annot_exon_info.start)
        if overlap > 0:
            # base coverage of long-read isoform by the annotated isoform
            frac_of_iso = float(overlap) / (exon.end - exon.start)
            # base coverage of the annotated isoform by the long-read isoform
            frac_of_annot = float(overlap) / (annot_exon_info.end - annot_exon_info.start)
            if frac_of_iso > MIN_ISOFORM_OVERLAP_FRAC and frac_of_annot > MIN_ANNOT_OVERLAP_FRAC:
                if annot_exon_info.name not in gene_hits or frac_of_iso > gene_hits[annot_exon_info.name][0]:
                    gene_hits[annot_exon_info.name] = [frac_of_iso, frac_of_annot]
    return gene_hits

def get_spliced_exon_overlaps(strand, exons, annots):
    gene_hits = []
    for annot_gene in annots.spliced_exons[strand]:
        annot_exons = sorted(list(annots.spliced_exons[strand][annot_gene]))
        # check if there is overlap in the genes
        # FIXME: not clear how this checks for overlap
        if (min((annot_exons[-1].end, exons[-1].end)) > max((annot_exons[0].start, exons[0].start))):
            covered_pos = set()
            for ex in exons:
                for aex in annot_exons:
                    for p in range(max((aex.start, ex.start)), min((aex.end, ex.end))):
                        covered_pos.add(p)
            if len(covered_pos) > sum([x.end - x.start for x in exons]) * 0.5:
                gene_hits.append([len(covered_pos), annot_gene, strand])
    return gene_hits

def _merge_intervals(intervals):
    "merge overlapping (start, end) intervals, sorted by start"
    merged = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return merged

def get_unspliced_exon_overlaps(strand, exons, annots):
    """Genes on the strand whose single-exon (unspliced) annotated transcripts cover
    more than half of the exons' length, as [bases covered, gene_id, strand]"""
    gene_to_intervals = {}
    for annot_exon in annots.all_annot_SE[strand]:
        gene_to_intervals.setdefault(annot_exon.name, []).append((annot_exon.start, annot_exon.end))
    exon_len = sum(e.end - e.start for e in exons)
    gene_hits = []
    for gene_id, intervals in gene_to_intervals.items():
        covered = sum(max(0, min(e.end, a_end) - max(e.start, a_start))
                      for e in exons for a_start, a_end in _merge_intervals(intervals))
        if covered > exon_len * 0.5:
            gene_hits.append([covered, gene_id, strand])
    return gene_hits

def fixed_end_exons(juncs):
    """exons of a spliced isoform from its junctions, with terminal exons of
    FALLBACK_TERMINAL_EXON_LEN rather than reaching its read ends"""
    exons = [Exon(max(0, juncs[0].start - FALLBACK_TERMINAL_EXON_LEN), juncs[0].start)]
    exons.extend(Exon(juncs[i].end, juncs[i + 1].start) for i in range(len(juncs) - 1))
    exons.append(Exon(juncs[-1].end, juncs[-1].end + FALLBACK_TERMINAL_EXON_LEN))
    return exons

def _get_transcript_gene_from_annot(iso_readrec, annots):
    """Return (transcript_id, gene_id) if iso matches an annotated junction chain, else (None, None).
    Each junction chain is named once, before end variants are split off; the variants
    keep this transcript_id."""
    if iso_readrec.juncs != () and iso_readrec.juncs in annots.juncchain_to_transcript:
        transcript_id, gene_id = annots.juncchain_to_transcript[iso_readrec.juncs]
        return transcript_id, (gene_id, )
    else:
        return None, None


def _find_gene_id_by_splicing(iso_readrec, annots):
    """Genes sharing junctions or splice sites with a spliced isoform, or, for a
    single-exon isoform, the best overlapping single-exon gene, as a tuple of gene
    ids, or None"""
    if iso_readrec.juncs != ():
        gene_hits = get_genes_with_shared_juncs(iso_readrec.juncs, iso_readrec.strand, annots)
        if gene_hits:
            return gene_hits
    else:
        gene_hits = get_single_exon_gene_overlaps(iso_readrec.strand, iso_readrec, annots)
        if gene_hits:
            return (sorted(gene_hits.items(), key=lambda x: x[1], reverse=True)[0][0], )
    return None


def _find_gene_id_by_exon_overlap(iso_readrec, annots):
    """The gene whose exons best overlap the isoform's, as a one-gene tuple, or None"""
    # look for exon overlap.  A spliced isoform's terminal exons
    # are given a fixed length from its terminal junctions (fixed_end_exons), so this
    # doesn't depend on its read ends.  There was an 'ambig' strand branch here; nothing
    # assigns that strand, the two branches were identical, and annots.spliced_exons is
    # keyed by '+' and '-' only, so it would have raised
    if iso_readrec.juncs != ():
        # spliced genes first, then single-exon genes, such as an unspliced lncRNA
        # the isoform is a spliced form of
        exons = fixed_end_exons(iso_readrec.juncs)
        gene_hits = (get_spliced_exon_overlaps(iso_readrec.strand, exons, annots)
                     or get_unspliced_exon_overlaps(iso_readrec.strand, exons, annots))
    else:
        gene_hits = get_spliced_exon_overlaps(iso_readrec.strand, iso_readrec.exons, annots)
    if gene_hits:
        gene_hits.sort(reverse=True)
        return (gene_hits[0][1], )
    else:
        return None


def _find_gene_id_by_overlap(iso_readrec, annots):
    """Find gene_id for an isoform without a matching junction chain, using junction or exon overlap."""
    # this all requires that we already trust the strand of the transcript
    # returns tuple of matching genes, will go into ref_gene_id field
    return (_find_gene_id_by_splicing(iso_readrec, annots)
            or _find_gene_id_by_exon_overlap(iso_readrec, annots))


def get_gene_name_unknown_strand(isoform, annots):
    """Gene identification for a spliced isoform of UNKNOWN_STRAND, whose junctions
    were only weakly stranded, which also gives its strand.  A matching annotated
    junction chain gives the gene.  Otherwise the isoform is tried on each strand,
    first by shared splicing and then by exon overlap; the first of these to find
    genes, all on one strand, gives the gene.  Returns (gene ids, transcript id,
    strand), with gene ids and strand None if no gene was found or genes were found
    on both strands."""
    transcript_id, gene_id = _get_transcript_gene_from_annot(isoform, annots)
    if transcript_id is not None:
        return gene_id, transcript_id, annots.gene_to_strand[gene_id[0]]
    for find_gene_id in (_find_gene_id_by_splicing, _find_gene_id_by_exon_overlap):
        strand_to_genes = {}
        for strand in ('+', '-'):
            isoform.strand = strand
            gene_id = find_gene_id(isoform, annots)
            if gene_id:
                # gene strand, as genes sharing a junction are found on either strand
                strand_to_genes.setdefault(annots.gene_to_strand[gene_id[0]], gene_id)
        isoform.strand = UNKNOWN_STRAND
        if len(strand_to_genes) == 1:
            strand, gene_id = strand_to_genes.popitem()
            return gene_id, None, strand
        elif len(strand_to_genes) > 1:
            break
    return None, None, None


def get_gene_name_firstpass(isoform, annots):
    transcript_id, gene_id = _get_transcript_gene_from_annot(isoform, annots)
    if transcript_id is None:
        gene_id = _find_gene_id_by_overlap(isoform, annots)
    return gene_id, transcript_id

def add_gene_isoform(genes, gene_id, isoform, strand, is_novel):
    # gene_id is a tuple of annotated gene ids, or one novel locus name.  A bare
    # ','.join put a comma between every character of the novel name
    key = gene_id if isinstance(gene_id, str) else ','.join(gene_id)
    hashed_id = int(hashlib.md5(key.encode('utf-8')).hexdigest(), 16)
    if is_novel:
        gene_id = ()
    if hashed_id not in genes:
        genes[hashed_id] = Gene(hashed_id, gene_id, isoform.chrom, strand)
    # this command also sets the gene_id in the isoform object
    genes[hashed_id].add_isoform(isoform)


def build_genes(firstpass, annots, region_chrom, sjc_with_overlap_groups):
    """Assign gene names to firstpass isoforms, building Gene objects.

    Returns (genes, novel_gene_isos_to_group) where:
    - genes: dict of gene_id -> Gene for isoforms matched to known genes
    - novel_gene_isos_to_group: isoforms needing novel gene assignment
    """
    genes = {}
    novel_gene_isos_to_group = {'+': [], '-': []}
    for iso_key in list(firstpass):
        isoform = firstpass[iso_key]
        if isoform.strand == UNKNOWN_STRAND:
            # only weakly stranded by its junctions; the strand has to come from a
            # gene, as there can't be isoforms of unknown strand in the output
            gene_id, isoform_id, strand = get_gene_name_unknown_strand(isoform, annots)
            if gene_id is None:
                logging.debug(f"isoform dropped: strand unknown and not resolved by gene identification: {iso_key}")
                firstpass.remove(iso_key)
                continue
            isoform.strand = strand
        else:
            gene_id, isoform_id = get_gene_name_firstpass(isoform, annots)
        isoform.ref_transcript_id = isoform_id
        if gene_id is not None:
            # removing this strand correction breaks the unusual junction (due to underlying variant?) test
            # currently just using first gene in list, ideally would use the most 5' annotated gene for this strand correction, if we want to do strand correction at all
            isoform.strand = annots.gene_to_strand[gene_id[0]]
            add_gene_isoform(genes, gene_id, isoform, annots.gene_to_strand[gene_id[0]], is_novel=False)
        else:
            novel_gene_isos_to_group[isoform.strand].append((isoform.start, isoform.end, iso_key))

    for strand in novel_gene_isos_to_group:
        generate_non_gene_iso_groups_strand(genes, novel_gene_isos_to_group, strand, region_chrom, sjc_with_overlap_groups)

    return genes


def _sweep_overlap_groups(spans):
    """group (start, end, key) spans whose coordinates overlap, transitively;
    returns lists of keys"""
    groups, group_end = [], None
    for start, end, key in sorted(spans):
        if groups and start < group_end:
            groups[-1].append(key)
            group_end = max(group_end, end)
        else:
            groups.append([key])
            group_end = end
    return groups

def _junction_span(isoform):
    return isoform.juncs[0].start, isoform.juncs[-1].end

def _group_by_read_span(isos, firstpass):
    "group isoforms by overlap of the spans between their read ends"
    return _sweep_overlap_groups([(firstpass[k].start, firstpass[k].end, k) for k in isos])

def _splice_site_components(spliced, firstpass):
    """connected components of spliced isoforms sharing a splice site, a shared
    junction sharing both of its sites"""
    parent = {k: k for k in spliced}

    def find(k):
        while parent[k] != k:
            parent[k] = parent[parent[k]]
            k = parent[k]
        return k
    site_owner = {}
    for k in spliced:
        for site in {site for j in firstpass[k].juncs for site in (j.start, j.end)}:
            if site in site_owner:
                parent[find(k)] = find(site_owner[site])
            else:
                site_owner[site] = k
    components = {}
    for k in spliced:
        components.setdefault(find(k), []).append(k)
    return list(components.values())

def _attach_single_exon_isos(single, groups, firstpass):
    """add each single-exon isoform to the group whose exons it overlaps most;
    returns those overlapping no group's exons"""
    group_exons = [[(e.start, e.end) for k in group for e in firstpass[k].exons] for group in groups]
    unattached = []
    for k in single:
        start, end = firstpass[k].start, firstpass[k].end
        overlaps = [sum(max(0, min(end, e_end) - max(start, e_start)) for e_start, e_end in exons)
                    for exons in group_exons]
        best = max(range(len(groups)), key=lambda i: overlaps[i], default=None)
        if best is not None and overlaps[best] > 0:
            groups[best].append(k)
        else:
            unattached.append(k)
    return unattached

def _group_novel_by_splicing(novel, firstpass):
    """Spliced isoforms sharing a junction or splice site are grouped, and groups
    whose junction spans overlap are merged, all without read ends.  A single-exon
    isoform joins the spliced group whose exons it overlaps most; single-exon
    isoforms overlapping no spliced group are grouped by their read-end spans."""
    spliced = [k for k in novel if firstpass[k].juncs]
    single = [k for k in novel if not firstpass[k].juncs]
    comp_spans = []
    for members in _splice_site_components(spliced, firstpass):
        spans = [_junction_span(firstpass[k]) for k in members]
        comp_spans.append((min(s for s, e in spans), max(e for s, e in spans), tuple(members)))
    groups = [[k for members in merged for k in members] for merged in _sweep_overlap_groups(comp_spans)]
    unattached = _attach_single_exon_isos(single, groups, firstpass)
    return groups + _group_by_read_span(unattached, firstpass)

def generate_non_gene_iso_groups_strand(genes, novel_gene_isos_to_group, strand, chrom, firstpass):
    """Group the strand's isoforms that have no annotated gene into novel genes, by
    _group_novel_by_splicing, and create their Gene objects.  A novel gene is named
    for the span of its isoforms."""
    novel = [iso_key for start, end, iso_key in novel_gene_isos_to_group[strand]]
    used_ids = set()
    for group in _group_novel_by_splicing(novel, firstpass):
        group_start = min(firstpass[k].start for k in group)
        group_end = max(firstpass[k].end for k in group)
        gene_id = f'{chrom}:{group_start}-{group_end}:{strand}'
        # groups can now share a span; a shared name would merge them
        suffix = 1
        while gene_id in used_ids:
            suffix += 1
            gene_id = f'{chrom}:{group_start}-{group_end}:{strand}.{suffix}'
        used_ids.add(gene_id)
        for k in group:
            add_gene_isoform(genes, gene_id, firstpass[k], strand, is_novel=True)
