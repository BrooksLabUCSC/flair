"""Clusters of junction chains sharing splice sites, the genes of their chains and
events, and the events of a region, found per cluster and written to its tables.
Used by flair spliceevents"""

from statistics import median
from collections import Counter
from flair.isoform_data import Isoform, Junc
from flair.gene_identification import novel_genes, group_by_splicing
from flair.interval_index import IntervalIndex
from flair.spliceevents_events import (END_GROUP_WINDOW, EventCaller, add_reads, alternative_site_events,
                                       complex_skipping_events, end_events, exon_skipping_events, extract_splicing_info,
                                       extract_terminal_info, group_ends, intron_retention_events, shared_site_events,
                                       terminal_exon_events, total)
from flair.spliceevents_regions import event_in_regions, gene_junction_span
from flair.spliceevents_tables import EventTables, is_expressed


class ReadCoverage:
    """The reads of a strand of a region, per sample, by their spans rounded to 10
    bases, for the reads covering an event's locus, whatever their gene"""
    def __init__(self):
        self.spans = {}
        self.index = None

    def add(self, start, end, sample):
        add_reads(self.spans, (round(start, -1), round(end, -1)), sample, 1)

    def covering(self, ref):
        "reads per sample whose spans contain ref, (start, end)"
        if self.index is None:
            self.index = IntervalIndex()
            for read_span, counts in self.spans.items():
                self.index.add(read_span[0], read_span[1], (read_span, counts))
            self.index.build()
        # a span containing ref overlaps the range between its ends, which can be in either
        # order, as a terminal exon event's (its introns' last start, first end)
        low, high = min(ref), max(ref)
        return total([counts for read_span, counts in self.index.overlap(low, max(high, low + 1))
                      if read_span[0] < ref[0] and ref[1] < read_span[1]])


def chain_isoforms(chrom, chains):
    """ReadChains as Isoforms, by (strand, junctions), with all their reads, spanning
    the median of their starts and of their ends"""
    isoforms = {}
    for chain in chains:
        key = (chain.strand, chain.juncs)
        if key not in isoforms:
            isoforms[key] = Isoform(chrom, chain.strand, tuple(Junc(*j) for j in chain.juncs))
        isoforms[key].reads.extend(chain.ends)
    for isoform in isoforms.values():
        isoform.start, isoform.end = int(median(isoform.starts)), int(median(isoform.ends))
    return isoforms


def _splicing_clusters(chain_isoforms):
    """The clusters of junction chains, (strand, junctions) keys of chain_isoforms,
    linked by shared splice sites (group_by_splicing), each named for its strand and
    order: each chain's cluster, and each cluster's strand.  Clusters whose junctions
    only overlap, as a gene nested in another's intron, are not merged: they share
    no event"""
    cluster_of, key_strand = {}, {}
    for strand in ('+', '-'):
        keys = sorted(k for k in chain_isoforms if k[0] == strand)
        for group in sorted(group_by_splicing(keys, chain_isoforms, merge_overlapping=False), key=lambda g: min(k[1][0][0] for k in g)):
            cluster = f'{strand}{len(key_strand)}'
            key_strand[cluster] = strand
            cluster_of.update((k, cluster) for k in group)
    return cluster_of, key_strand


def group_chains(chrom, samples, sample_chains, novel_gene):
    """Each sample's reads grouped three ways, each a list per sample of
    key -> junctions -> [ReadEnds]: by cluster of junction chains linked by shared
    splice sites, on one strand, for the events of splicing, which don't depend on
    genes; by each chain's 5'-most gene, for the first exons and starts; and by its
    3'-most gene, for the last exons and ends.  A gene's 5' fragments are with its
    starts and its 3' fragments with its ends, under no junctions, and nowhere else.
    Also each key's strand, each chain's genes, by (strand, junctions), a novel gene's
    if of none annotated, and each strand's ReadCoverage, of its spliced reads"""
    chain_genes = {}
    coverage = {'+': ReadCoverage(), '-': ReadCoverage()}
    for sample, chains in zip(samples, sample_chains):
        for chain in chains:
            if not chain.fragment:
                key = (chain.strand, chain.juncs)
                chain_genes[key] = chain.genes or (novel_gene[key],)
                for read in chain.ends:
                    coverage[chain.strand].add(read.start, read.end, sample)
    isoforms = chain_isoforms(chrom, (chain for chains in sample_chains for chain in chains if not chain.fragment))
    cluster_of, key_strand = _splicing_clusters(isoforms)
    return _chain_views(sample_chains, chain_genes, cluster_of, key_strand), key_strand, chain_genes, coverage


def _chain_views(sample_chains, chain_genes, cluster_of, key_strand):
    "group_chains' three groupings of each sample's chains, adding the genes' strands to key_strand"
    views = {'cluster': [], 'five': [], 'three': []}
    for chains in sample_chains:
        by = {'cluster': {}, 'five': {}, 'three': {}}
        for chain in chains:
            if chain.fragment:
                gene = chain.genes[0]
                by['five' if chain.fragment == 'first' else 'three'].setdefault(gene, {})[()] = chain.ends
                key_strand[gene] = chain.strand
                continue
            key = (chain.strand, chain.juncs)
            genes = chain_genes[key]
            for view, group in (('cluster', cluster_of[key]), ('five', genes[0]), ('three', genes[-1])):
                by[view].setdefault(group, {})[chain.juncs] = chain.ends
                key_strand[group] = chain.strand
        for view in views:
            views[view].append(by[view])
    return views


class AnnotatedSplicing:
    """The annotated introns and splice sites of each strand, and the genes they are
    annotated in, to give events their genes and tell their novel junctions"""
    def __init__(self, annots=None):
        self.intron_genes = {'+': {}, '-': {}}
        self.site_genes = {'+': {}, '-': {}}
        self.num_transcripts = {}
        self.names = {}
        if annots is not None:
            self.num_transcripts = annots.gene_to_num_transcripts
            self.names = annots.gene_to_name
            for gene, juncs in annots.gene_to_annot_juncs.items():
                strand = annots.gene_to_strand.get(gene)
                if strand not in self.intron_genes:
                    continue
                for j in juncs:
                    self.intron_genes[strand].setdefault((j.start, j.end), set()).add(gene)
                    for site in (j.start, j.end):
                        self.site_genes[strand].setdefault(site, set()).add(gene)

    def is_annotated(self, junc, strand):
        return junc in self.intron_genes[strand]

    def gene_names(self, genes):
        "the annotation's names of genes, a list, comma-separated, empty if it names none of them"
        names = [self.names.get(gene, '') for gene in genes]
        return ','.join(names) if any(names) else ''


def _top(scores):
    "the keys of the highest score, several on a tie"
    if not scores:
        return []
    best = max(scores.values())
    return sorted(key for key, score in scores.items() if score == best)


def event_gene_scores(introns, strand, annotated, junc_reads):
    """the genes sharing an intron or splice site with an event, on its strand, and their
    scores: 2 for an intron the gene has annotated, or failing that 1 for each of its
    splice sites the gene has, each intron weighted by its share of the reads of all
    the introns, junc_reads, so a rare alternative counts little"""
    total_reads = sum(junc_reads.get(junc, 0) for junc in introns)
    scores = Counter()
    for junc in introns:
        weight = junc_reads.get(junc, 0) / total_reads if total_reads > 0 else 1
        intron_genes = annotated.intron_genes[strand].get(junc, set())
        for gene in intron_genes:
            scores[gene] += 2 * weight
        for site in junc:
            for gene in annotated.site_genes[strand].get(site, set()) - intron_genes:
                scores[gene] += weight
    # genes of the same introns tie, whatever order their weights were added in
    return {gene: round(score, 9) for gene, score in scores.items()}


def annotated_event_genes(scores, annotated):
    """the genes of the highest event_gene_scores.  A tie goes to the gene of most
    annotated transcripts, so a gene and a readthrough of it give the gene; genes still
    tied are all given"""
    best = _top(scores)
    return _top({gene: annotated.num_transcripts[gene] for gene in best})


def read_event_genes(introns, cluster_chains, chain_genes, strand):
    """the genes of most reads of an event's introns, by the genes the reads were given
    in correction, each gene of a readthrough counted: for an event of no annotated
    intron or splice site"""
    reads = Counter()
    for chains in cluster_chains:
        for juncs, ends in chains.items():
            if introns.intersection(juncs):
                for gene in chain_genes[(strand, juncs)]:
                    reads[gene] += len(ends)
    return _top(reads)


def _keys(view):
    return sorted(set().union(*[sample_groups.keys() for sample_groups in view]))


def gene_terminal_events(caller, view, gene, terminal, end_type, annot_sites):
    """A gene's first (afe) or last (ale) exon event, and start (tss) or end (tts)
    event, from its chains in view, and the first or last exon event's junctions"""
    juncs, terminal_exons, ends = extract_terminal_info(caller.samples, view, gene, caller.strand, terminal)
    events, event_juncs = terminal_exon_events(caller, terminal, terminal_exons, juncs, annot_sites)
    events += end_events(caller, end_type, group_ends(ends, END_GROUP_WINDOW))
    for event in events:
        event.gene = gene
    return events, event_juncs


def cluster_events(caller, cluster_chains, terminal_juncs, annotated, chain_genes):
    """The splicing events of a cluster of junction chains, each sample's, each event
    given the gene annotated with most of its introns, those of all its alternatives,
    and their splice sites, weighted by their reads (annotated_event_genes), or with
    none, that of most of its reads (read_event_genes); several only on a tie"""
    info = extract_splicing_info(caller.samples, cluster_chains, caller.strand)
    es, es_juncs, seen = exon_skipping_events(caller, info)
    ces = complex_skipping_events(caller, info, es_juncs, seen)
    alt3, alt3_groups = alternative_site_events(caller, 'alt3', info.ss5to3, es_juncs, terminal_juncs)
    alt5, alt5_groups = alternative_site_events(caller, 'alt5', info.ss3to5, es_juncs, terminal_juncs)
    events = (es + ces + alt3 + shared_site_events(caller, 'alt5ss', alt3_groups)
              + alt5 + shared_site_events(caller, 'alt3ss', alt5_groups) + intron_retention_events(caller, info))
    junc_reads = {junc: sum(counts.values()) for junc, counts in info.juncs.items()}
    for event in events:
        scores = event_gene_scores(event.introns, caller.strand, annotated, junc_reads)
        genes = annotated_event_genes(scores, annotated) or read_event_genes(event.introns, cluster_chains, chain_genes, caller.strand)
        event.gene = ','.join(genes)
        # the other genes sharing a splice site, best first
        event.additional_genes = sorted(set(scores) - set(genes), key=lambda gene: (-scores[gene], gene))
    return events


def process_region_events(temp_prefix, chrom, samples, sample_chains, novel_gene, thresholds, annot_afe_ss, annot_ale_ss,
                          report_regions, gene_juncs, annotated=None, owned=None):
    """A region's events, from each sample's ReadChains, written to its tables.  The
    events of splicing are found per cluster of junction chains linked by shared
    splice sites, without genes (cluster_events); the alternative first exons and
    starts per each chain's 5'-most gene, and the last exons and ends per its 3'-most
    gene, so a readthrough's reads count with the gene each end lies in.  An event's
    locus reads are all the strand's reads covering it.  Only the events inside
    report_regions, and given owned, starting in it, are written (event_in_regions)"""
    annotated = annotated or AnnotatedSplicing()
    views, key_strand, chain_genes, coverage = group_chains(chrom, samples, sample_chains, novel_gene)
    callers = {strand: EventCaller(chrom, strand, samples, thresholds, coverage[strand]) for strand in ('+', '-')}
    events = []  # (event, of a GENE_WIDE_EVENTS type the span of its gene's junctions, for the regions)
    terminal_juncs = {'+': set(), '-': set()}
    for view, terminal, annot_ss, end_type in (('five', 'afe', annot_afe_ss, 'tss'), ('three', 'ale', annot_ale_ss, 'tts')):
        for gene in _keys(views[view]):
            caller = callers[key_strand[gene]]
            gene_events, event_juncs = gene_terminal_events(caller, views[view], gene, terminal, end_type,
                                                            annot_ss.get(caller.strand, []))
            terminal_juncs[caller.strand] |= event_juncs
            gene_span = gene_junction_span(gene_juncs[gene])
            events.extend((event, gene_span) for event in gene_events)
    for cluster in sorted(_keys(views['cluster']), key=lambda c: int(c[1:])):
        caller = callers[key_strand[cluster]]
        cluster_chains = [chains.get(cluster, {}) for chains in views['cluster']]
        events.extend((event, None) for event in cluster_events(caller, cluster_chains, terminal_juncs[caller.strand],
                                                                annotated, chain_genes))
    reported = {}
    for event, gene_span in events:
        if event_in_regions(event, report_regions, gene_span, owned):
            if event.name in reported:
                raise AssertionError(f'two events named {event.name}, of {reported[event.name].gene} and {event.gene}')
            reported[event.name] = event
    with EventTables(temp_prefix, samples, thresholds.event_support) as tables:
        for event in reported.values():
            if is_expressed(event, samples, thresholds):
                tables.write(event, annotated)


def novel_genes_of_chains(chrom, sample_chains):
    """(strand, junctions) -> novel gene of the ReadChains of no annotated gene, of all
    samples together, so each sample's reads get the same novel genes, grouped and
    named as flair transcriptome does (novel_genes)"""
    novel = chain_isoforms(chrom, (chain for chains in sample_chains for chain in chains if not chain.genes))
    novel_gene = {}
    for strand in ('+', '-'):
        for gene, keys in novel_genes(sorted(k for k in novel if k[0] == strand), novel, chrom, strand):
            novel_gene.update((k, gene) for k in keys)
    return novel_gene


def gene_junctions(sample_chains, novel_gene, annots):
    """Each gene's junctions, for --region_bed: an annotated gene's annotated ones, or
    those of its reads for a novel gene, or an annotated one with none, of one exon"""
    read_juncs = {}
    for chains in sample_chains:
        for chain in chains:
            for gene in chain.genes or (novel_gene[(chain.strand, chain.juncs)],):
                read_juncs.setdefault(gene, set()).update(chain.juncs)
    return {gene: {(j.start, j.end) for j in annots.gene_to_annot_juncs.get(gene, ())} or juncs
            for gene, juncs in read_juncs.items()}
