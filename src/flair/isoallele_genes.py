"""One gene's isoform alleles, for flair isoalleles: its reads' groups, by haplotype
and somatic set, and the sequence and protein of each group of each isoform.

Each sample's reads over the gene are read once (_scan_sample): the reads flair
quantify assigned to one of the gene's isoforms (the sample's read map) give their
haplotype, from the HP tag, and their alleles at the isoform's candidate somatic
sites and unphased germline hets, and in a pair the normal's reads give the depth at
each candidate somatic site.

Somatic sites, per isoform (site_class): flair quantify assigns a read to a
multi-exon isoform only if it is aligned SPLICE_SITE_FLANK bases into the first and
last exons, so every assigned read covers the internal exons and the terminal
exons next to their splice sites ('guaranteed'), but only some reach the rest of the
terminal exons, such as UTRs ('terminal').  The guaranteed sites define the
somatic sets; if an isoform has none, its terminal sites do.  An isoform's other
terminal sites are added to a group's sequence if the group's reads that reach the
site all carry it, and at least one does (terminal_sites.tsv).

A read's group is its haplotype (H1, H2 or untagged) and the set of defining
somatic sites it carries, numbered per gene by reads, most first: H2-S1.  A read
that doesn't inform every defining site of its isoform may carry the one it misses,
so its set isn't known: its group is H2-S? and gets no sequence.
"""
import collections
import copy
import io
from flair.count_sam_transcripts import SPLICE_SITE_FLANK
from flair.isoallele_sequences import (build_transcript, call_cds, exon_index, frameshift_residue, protein_change,
                                       protein_id, variant_label)
from flair.isoallele_variants import (HOMOZYGOUS, PHASED_HET, UNTAGGED, ReadAlleles, called_variants, fetch,
                                      germline_variants, haplotype)

GUARANTEED, TERMINAL = 'guaranteed', 'terminal'
HAPLOTYPES = ('H1', 'H2')
UNKNOWN_SET = 'S?'


def exonic(iso, key):
    "whether a variant lies wholly inside one of the isoform's exons"
    return exon_index(sorted((b.start, b.end) for b in iso.blocks), key.pos, key.ref) is not None


def site_class(iso, key):
    """GUARANTEED if every read flair quantify assigns to the isoform covers a variant,
    TERMINAL if only those reaching far enough into a terminal exon do, None if it
    isn't in an exon.  A single-exon isoform needs only half its length covered, so all
    of its sites are terminal."""
    start, end = key.pos - 1, key.pos - 1 + len(key.ref)
    blocks = sorted((b.start, b.end) for b in iso.blocks)
    i = next((j for j, (s, e) in enumerate(blocks) if s <= start and end <= e), None)
    if i is None:
        return None
    if len(blocks) == 1:
        return TERMINAL
    if 0 < i < len(blocks) - 1:
        return GUARANTEED
    if i == 0:  # leftmost exon: guaranteed only next to its splice site, on its right
        return GUARANTEED if start >= blocks[0][1] - SPLICE_SITE_FLANK else TERMINAL
    return GUARANTEED if end <= blocks[-1][0] + SPLICE_SITE_FLANK else TERMINAL


class GeneOutput:
    "the rows of every output, of one gene or a batch of genes"

    TABLES = ('counts', 'sequences', 'somatic_sets', 'variants', 'excluded', 'warnings', 'terminal_sites', 'genes')

    def __init__(self):
        for table in self.TABLES:
            setattr(self, table, [])
        self.fasta = []
        self.bed = []
        self.proteins = {}  # protein id -> (aaseq, [isoform|group])
        self.labels = collections.defaultdict(list)  # sample -> [(read, isoform|group)]


ReadInfo = collections.namedtuple('ReadInfo', ('isoform', 'haplotype', 'alleles'))


def _scan_sample(ctx, sample, region, isoforms, candidates, hets, depth_sites):
    """One pass over a sample's reads in the gene's region: each read assigned to one of
    the isoforms as ReadInfo, its alleles at the isoform's candidate somatic sites and
    unphased hets; the reads informing each of depth_sites; and the gene's spliced
    reads and those assigned, for genes.tsv.  Also returns the assigned reads not in
    the bam."""
    chrom, start, end = region
    assigned = ctx.map_readers[sample].assigned(isoforms)
    strand = collections.Counter(i.strand for i in isoforms.values()).most_common(1)[0][0]
    exons = sorted({(b.start, b.end) for i in isoforms.values() for b in i.blocks})
    sites_of = {name: [k for k in candidates if exonic(iso, k)] + [h for h in hets if exonic(iso, h)]
                for name, iso in isoforms.items()}
    reads, depth, seen = {}, collections.Counter(), set()
    spliced = assigned_spliced = 0
    for read in ctx.bams[sample].fetch(chrom, start, end):
        if read.is_secondary or read.is_supplementary or read.is_unmapped:
            continue
        alleles = ReadAlleles(read)
        for key in depth_sites:
            if read.reference_start < key.pos + len(key.ref) and key.pos - 1 < read.reference_end \
                    and alleles.allele(key) is not None:
                depth[key] += 1
        name = read.query_name
        if read.is_reverse == (strand == '-') and any(op == 3 for op, _ in read.cigartuples) \
                and any(s < e2 and s2 < e for s, e in read.get_blocks() for s2, e2 in exons):
            spliced += 1
            assigned_spliced += name in assigned
        iso = assigned.get(name)
        if iso is None or name in seen:
            continue
        seen.add(name)
        reads[name] = ReadInfo(iso, haplotype(read), {k: alleles.allele(k) for k in sites_of[iso]})
    coverage = (spliced, assigned_spliced, len(assigned))
    return reads, depth, coverage, len(set(assigned) - seen)


def _somatic_candidates(ctx, region, germline):
    "the gene's somatic variants before the normal depth check, by ctx.somatic_mode"
    chrom, start, end = region
    if ctx.somatic_mode == 'variant':
        return [k for k in ctx.variant_keys if k.chrom == chrom and start < k.pos <= end]
    if ctx.somatic_mode == 'somatic_vcf':
        return sorted(called_variants(fetch(ctx.somatic_vcf, ctx.somatic_vcf_path, chrom, start, end)))
    if ctx.somatic_mode == 'tumor_vcf':
        normal_calls = {g.key for g in germline}
        return sorted(called_variants(fetch(ctx.tumor_vcf, ctx.tumor_vcf_path, chrom, start, end)) - normal_calls)
    return []


def _place_unphased(ctx, key, iso_name, reads_by_sample):
    """An unphased het's haplotype on one isoform, from the isoform's tagged reads, in
    every sample, that carry it: the haplotype of at least ctx.carrier_min of them, of
    at least ctx.min_tagged in all.  None if the reads can't place it."""
    alt_reads = collections.Counter(info.haplotype for reads in reads_by_sample.values() for info in reads.values()
                                    if info.isoform == iso_name and info.haplotype != UNTAGGED
                                    and info.alleles.get(key) == 'alt')
    tagged = sum(alt_reads.values())
    if tagged >= ctx.min_tagged:
        top, n = min(alt_reads.items(), key=lambda x: (-x[1], x[0]))
        if n / tagged >= ctx.carrier_min:
            return top
    return None


def process_gene(ctx, gene_id, isos, out):  # noqa: C901 - one gene's whole analysis, in its order
    """One gene's tables, sequences and warnings, added to out (GeneOutput).  ctx is the
    worker's open files and settings (flair_isoalleles.WorkerContext)."""
    region = (isos[0].chrom, min(i.chromStart for i in isos), max(i.chromEnd for i in isos))
    isoforms = {i.name: i for i in isos}
    phased = fetch(ctx.phased_vcf, ctx.phased_vcf_path, *region)
    if ctx.germline_vcf is None:
        germline = germline_variants(phased)
    else:
        germline = germline_variants(fetch(ctx.germline_vcf, ctx.germline_vcf_path, *region), phased)
    candidates = _somatic_candidates(ctx, region, germline)
    candidate_set = set(candidates)
    # a somatic variant is not also germline, as a --variant the normal calls
    germline = [g for g in germline if g.key not in candidate_set]
    hets = [g.key for g in germline if g.kind not in (HOMOZYGOUS, PHASED_HET)]

    reads_by_sample, coverage, missing = {}, {}, {}
    depth = collections.Counter()
    for sample in ctx.samples:  # the normal first, in a pair: its depth decides the somatic sites
        depth_sites = candidates if sample == ctx.normal else []
        reads, sample_depth, coverage[sample], missing[sample] = _scan_sample(
            ctx, sample, region, isoforms, candidates, hets, depth_sites)
        reads_by_sample[sample] = reads
        depth.update(sample_depth)

    if ctx.somatic_mode == 'tumor_vcf':
        somatic = [k for k in candidates if depth[k] >= ctx.min_normal_depth]
        for key in candidates:
            if depth[key] < ctx.min_normal_depth:
                out.excluded.append([gene_id, variant_label(key), depth[key],
                                     ','.join(n for n, i in isoforms.items() if exonic(i, key)),
                                     f'tumor-only call with {depth[key]} normal reads (< {ctx.min_normal_depth}): '
                                     'germline can\'t be ruled out'])
    else:
        somatic = candidates

    for sample in ctx.output_samples:
        spliced, assigned_spliced, n_assigned = coverage[sample]
        share = round(assigned_spliced / spliced, 3) if spliced else ''
        out.genes.append([gene_id, sample, spliced, assigned_spliced, share, n_assigned])
        if missing[sample]:
            out.warnings.append([gene_id, '', 'reads not in bam',
                                 f'{sample}: {missing[sample]} assigned reads not found in the bam'])

    # each isoform's somatic sites that define its sets, and its other terminal ones
    defining, terminal = {}, {}
    for name, iso in isoforms.items():
        classes = {k: site_class(iso, k) for k in somatic}
        guaranteed = [k for k, c in classes.items() if c == GUARANTEED]
        term = [k for k, c in classes.items() if c == TERMINAL]
        defining[name], terminal[name] = (guaranteed, term) if guaranteed else (term, [])
        if not guaranteed:
            for key in term:
                out.terminal_sites.append([gene_id, name, variant_label(key), '', '', '',
                                           'defines sets (the isoform has no guaranteed somatic site)'])

    def carried(info):
        "the defining sites a read carries, sorted, and whether it misses one"
        sites = defining[info.isoform]
        return (tuple(sorted(k for k in sites if info.alleles.get(k) == 'alt')),
                any(info.alleles.get(k) is None for k in sites))

    set_reads = collections.Counter(s for reads in reads_by_sample.values() for info in reads.values()
                                    for s, uncovered in (carried(info),) if s and not uncovered)
    set_ids = {s: f'S{i}' for i, (s, _) in enumerate(sorted(set_reads.items(), key=lambda x: (-x[1], x[0])), 1)}
    sets_by_id = {sid: s for s, sid in set_ids.items()}
    for s, sid in set_ids.items():
        out.somatic_sets.append([gene_id, sid, ','.join(variant_label(k) for k in s), set_reads[s]])

    def group(info):
        sites, uncovered = carried(info)
        if uncovered:
            return f'{info.haplotype}-{UNKNOWN_SET}'
        return f'{info.haplotype}-{set_ids[sites]}' if sites else info.haplotype

    group_reads = collections.defaultdict(collections.Counter)  # (isoform, group) -> {sample: reads}
    groups_seen = collections.defaultdict(set)
    for sample, reads in reads_by_sample.items():
        partial = collections.Counter()
        for read, info in reads.items():
            g = group(info)
            group_reads[(info.isoform, g)][sample] += 1
            groups_seen[info.isoform].add(g)
            partial[info.isoform] += g.endswith(UNKNOWN_SET)
            if ctx.keep_labels:
                out.labels[sample].append((read, f'{info.isoform}|{g}'))
        for iso, n in sorted(partial.items()):
            if n:
                out.warnings.append([gene_id, iso, 'reads missing somatic sites',
                                     f'{sample}: {n} reads assigned to the isoform don\'t cover every somatic site '
                                     'that defines its sets, so their set isn\'t known: they are in -S? groups, '
                                     'which get no sequence'])
    for (iso, g), counts in sorted(group_reads.items()):
        out.counts.append([iso, f'{iso}|{g}', [counts[s] for s in ctx.output_samples]])

    low_depth = [k for k in somatic if ctx.normal is not None and depth[k] < ctx.min_normal_depth]
    reference_ids = _reference_protein_ids(ctx, isoforms)
    for name, iso in isoforms.items():
        if not groups_seen[name]:
            continue  # no reads: no sequences
        _isoform_sequences(ctx, gene_id, name, iso, germline, reads_by_sample, terminal, sets_by_id,
                           groups_seen[name], group, group_reads, low_depth, reference_ids, out)

    for g in germline:
        if g.kind in (HOMOZYGOUS, PHASED_HET):
            in_isos = [n for n, i in isoforms.items() if exonic(i, g.key)]
            if in_isos:
                out.variants.append([gene_id, ','.join(in_isos), variant_label(g.key), f'germline {g.kind}',
                                     g.haplotype or 'both', g.phase_set if g.phase_set is not None else ''])
    for key in somatic:
        in_isos = [n for n, i in isoforms.items() if exonic(i, key)]
        sids = [sid for s, sid in set_ids.items() if key in s]
        out.variants.append([gene_id, ','.join(in_isos), variant_label(key), 'somatic', ','.join(sids) or 'no reads',
                             f'normal reads {depth[key]}' if ctx.normal is not None else ''])


def _haplotype_alleles(ctx, gene_id, name, iso, germline, reads_by_sample, out):
    """The isoform's germline variants on each haplotype, untagged reads' sequences
    having only the homozygous ones, and the phase blocks of its phased hets"""
    alleles = {'H1': [], 'H2': [], UNTAGGED: []}
    blocks = set()
    for g in germline:
        if not exonic(iso, g.key):
            continue
        if g.kind == HOMOZYGOUS:
            for hap in alleles:
                alleles[hap].append(g.key)
        elif g.kind == PHASED_HET:
            if g.haplotype in alleles:
                alleles[g.haplotype].append(g.key)
            if g.phase_set is not None:
                blocks.add(g.phase_set)
        else:
            placed = _place_unphased(ctx, g.key, name, reads_by_sample)
            if placed in alleles:
                alleles[placed].append(g.key)
                out.variants.append([gene_id, name, variant_label(g.key), 'germline unphased het, placed by reads',
                                     placed, ''])
            else:
                out.warnings.append([gene_id, name, 'unplaced het',
                                     f'{variant_label(g.key)}: the reads can\'t place it, the reference allele is used'])
                out.variants.append([gene_id, name, variant_label(g.key), 'germline unphased het, not placed', '', ''])
    if len(blocks) > 1:
        out.warnings.append([gene_id, name, 'several phase blocks',
                             f'phase blocks {",".join(map(str, sorted(blocks)))}: H1 and H2 are only consistent '
                             'within each'])
    return alleles, blocks


def _terminal_added(gene_id, name, terminal, wanted, reads_by_sample, group, out):
    """the terminal-exon somatic sites added to each group's sequence: those the
    group's reads that reach them all carry, at least one"""
    added = collections.defaultdict(list)
    for key in terminal[name]:
        for g in sorted(wanted):
            seen = [info.alleles.get(key) for reads in reads_by_sample.values() for info in reads.values()
                    if info.isoform == name and group(info) == g and info.alleles.get(key) is not None]
            n_alt = seen.count('alt')
            decision = ('added' if seen and n_alt == len(seen) else
                        'dropped (no read reaches it)' if not seen else 'dropped (reads without it)')
            if decision == 'added':
                added[g].append(key)
            out.terminal_sites.append([gene_id, name, variant_label(key), g, len(seen), n_alt, decision])
    return added


def _reference_protein_ids(ctx, isoforms):
    """the protein ids of the gene's isoforms, by protein: a variant sequence encoding
    the protein of another of the gene's isoforms gets its id"""
    ids = {}
    for iso in isoforms.values():
        if iso.aaseq_id:
            seq, origins, sizes, _ = build_transcript(ctx.genome, iso, [])
            protein = call_cds(iso, seq, origins, sizes, False)[1]
            if protein:
                ids.setdefault(protein, iso.aaseq_id)
    return ids


def _isoform_sequences(ctx, gene_id, name, iso, germline, reads_by_sample, terminal, sets_by_id,
                       groups, group, group_reads, low_depth, reference_ids, out):
    "one isoform's sequences: H1, H2, and every somatic group its reads show"
    alleles, blocks = _haplotype_alleles(ctx, gene_id, name, iso, germline, reads_by_sample, out)
    for key in low_depth:
        if exonic(iso, key):
            out.warnings.append([gene_id, name, 'LOW NORMAL COVERAGE AT SOMATIC SITE',
                                 f'{variant_label(key)}: fewer than {ctx.min_normal_depth} normal reads, germline '
                                 'can\'t be ruled out; kept as somatic, as --variant or --somatic_vcf gives it'])

    wanted = set(HAPLOTYPES) | {g for g in groups if '-S' in g and not g.endswith(UNKNOWN_SET)}
    added = _terminal_added(gene_id, name, terminal, wanted, reads_by_sample, group, out)

    # the reference CDS and protein, to compare each group's with.  An isoform the annotation
    # calls productive although its stop is premature stays exempt from the premature stop rule
    ref_seq, ref_origins, ref_sizes, _ = build_transcript(ctx.genome, iso, [])
    ref_prod, ref_protein, _, _ = call_cds(iso, ref_seq, ref_origins, ref_sizes, False)
    exempt = iso.productivity == 'PRO' and ref_prod == 'PTC'
    germline_protein = {}
    for g in sorted(wanted):
        hap, _, sid = g.partition('-')
        variants = alleles.get(hap, []) + list(sets_by_id.get(sid, ())) + added[g]
        seq, origins, sizes, skipped = build_transcript(ctx.genome, iso, variants)
        for key in skipped:
            out.warnings.append([gene_id, name, 'variant not applied',
                                 f'{g}: {variant_label(key)} isn\'t wholly in one exon, overlaps another variant, '
                                 'or its ref allele doesn\'t match the genome'])
        applied = [k for k in variants if k not in skipped]
        productivity, aaseq, thick, note = call_cds(iso, seq, origins, sizes, exempt)
        change = protein_change(ref_protein, aaseq, frameshift_residue(iso, applied))
        somatic_applied = [k for k in applied if k not in alleles.get(hap, [])]
        if hap not in germline_protein:
            germline_protein[hap] = call_cds(iso, *build_transcript(ctx.genome, iso, alleles.get(hap, []))[:3],
                                             exempt)[1]
        somatic_change = (protein_change(germline_protein[hap], aaseq, frameshift_residue(iso, somatic_applied))
                          if somatic_applied else 'none')
        pid = protein_id(aaseq, reference_ids) if aaseq else ''
        if aaseq:
            out.proteins.setdefault(pid, (aaseq, []))[1].append(f'{name}|{g}')
        seq_name = f'{name}|{g}'
        main = [variant_label(k) for k in applied if k not in added[g]]
        term_added = [variant_label(k) for k in added[g] if k not in skipped]
        out.fasta.append(f'>{seq_name}\n{seq}\n')
        out.sequences.append([gene_id, name, g, ','.join(main), ','.join(term_added), productivity, note, pid,
                              change, somatic_change,
                              f'{iso.chrom}:{thick[0]}-{thick[1]}' if thick[1] > thick[0] else '',
                              ','.join(map(str, sorted(blocks)))])
        bed = copy.deepcopy(iso)
        bed.name = seq_name
        bed.source_isoform, bed.allele_group = name, g
        bed.thickStart, bed.thickEnd = thick
        bed.productivity, bed.aaseq_id = productivity, pid or None
        bed.transcript_class = 'with_vars' if applied else iso.transcript_class
        reads = group_reads[(name, g)]
        bed.read_support, bed.frac_support = sum(reads.values()), None
        bed.samples = [s for s in ctx.output_samples if reads[s]]
        bed.itemRgb = bed._get_rgb()
        fh = io.StringIO()
        bed.write(fh)
        out.bed.append(fh.getvalue())
