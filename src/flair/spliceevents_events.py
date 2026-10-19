"""Splicing events of one cluster of junction chains, by type, from the reads per
sample of its junctions, splice sites, exons and ends.  Used by flair spliceevents"""

from bisect import bisect_left
from statistics import median
from collections import namedtuple, Counter


# the thresholds of event calling, from the command line
Thresholds = namedtuple('Thresholds', ('feature_support', 'event_frac_of_tot', 'junc_frac_of_event', 'event_support'))
# an alternative whose PSI doesn't vary is a row of its event, not part of its other
# row, when it has this part of the event's reads in some sample, as the main
# alternative often does
MIN_STEADY_FRACTION = 0.1


def thresholds_of(args):
    return Thresholds(args.feature_support, args.event_frac_of_tot, args.junc_frac_of_event, args.event_support)


# Reads are counted per sample in Counters, sample -> reads, without the samples of none


def add_reads(counts, key, sample, numreads):
    "add a sample's reads to counts[key]"
    if key not in counts:
        counts[key] = Counter()
    key_counts = counts[key]
    key_counts[sample] = key_counts.get(sample, 0) + numreads


def add_counts(counts, key, more):
    "add more, reads per sample, to counts[key]"
    if key not in counts:
        counts[key] = Counter()
    counts[key].update(more)


def total(all_counts):
    "the sum of reads per sample"
    result = Counter()
    for counts in all_counts:
        result.update(counts)
    return result


def psi_varies(psis, junc_frac_of_event):
    """Is a junction, splice site or end of an event alternative, from its PSIs in the
    samples with enough reads of the event: differing across them by
    junc_frac_of_event, or, in the only sample with a PSI, at least junc_frac_of_event
    from both 0 and 1, so an event in one sample, as in a small dataset, is reported"""
    if len(psis) == 0:
        return False
    if len(psis) == 1:
        return junc_frac_of_event <= psis[0] <= 1 - junc_frac_of_event
    return max(psis) - min(psis) >= junc_frac_of_event


def sample_psis(counts, totals, samples, min_total):
    "the PSIs, counts over totals, of the samples with at least min_total reads in totals"
    return [counts[s] / totals[s] for s in samples if totals[s] >= min_total]


def junc_between(ss1, ss2):
    "the junction between two splice sites"
    return (min(ss1, ss2), max(ss1, ss2))


def lies_left(upstream, strand):
    """whether a feature lies left of another, at lower coordinates, on strand: if
    upstream (5') of it on the + strand, or downstream (3') of it on the - strand"""
    return strand == '+' if upstream else strand == '-'


def internal_exons(juncs):
    "the exons between a chain's junctions"
    return [(juncs[i][1], juncs[i + 1][0]) for i in range(len(juncs) - 1)]


def event_name(etype, chrom, start, end, strand):
    """an event's name: its type and outer bounds, start 0-based, written 1-based as a
    genome browser's locus.  Not its gene, which may be several, or decided only once
    the event is found"""
    return f'{etype}|{chrom}:{start + 1}-{end}({strand})'


def row_name(ename, alternative):
    "an event row's name, in BED records: the event's and the row's alternative"
    return f'{ename}|{alternative}'


def coord_ranges(ranges):
    "introns or exons, (start, end) with 0-based starts, written 1-based"
    return ','.join(f'{start + 1}-{end}' for start, end in sorted(ranges))


def span(ranges):
    "the outer bounds of introns or exons"
    return min(r[0] for r in ranges), max(r[1] for r in ranges)


def junction_blocks(junc, exons=()):
    "the BED blocks of a junction: 10 bases each side of it, and the exons included between its ends"
    return [(junc[0] - 10, junc[0]), *exons, (junc[1], junc[1] + 10)]


# the BED colors of each event type: first exons as starts, last exons as ends, and
# the sites sharing alternative ones as the alternative ones
EVENT_COLORS = {'es': '66,105,245', 'ces': '66,105,245', 'ir': '144,66,245',
                'alt3': '245,215,66', 'alt3ss': '245,215,66', 'alt5': '43,184,39', 'alt5ss': '43,184,39',
                'tss': '255,0,0', 'afe': '255,0,0', 'tts': '0,0,255', 'ale': '0,0,255'}


# an event's row: its reads per sample, its included and excluded junctions (or ends),
# the skipping junctions of an exon skip, its exons, and its BED blocks, a list per
# BED record
EventRow = namedtuple('EventRow', ('counts', 'inc_juncs', 'exc_juncs', 'outer_juncs', 'exons', 'blocks'))


class SplicingEvent:
    """An event: its rows, EventRows by alternative, its reads per sample (totjunc),
    the reads covering its locus (totoverlap), and the reads of its other row, those
    of the alternatives too rare, or invariant and too small a part of it, to be rows"""
    def __init__(self, name, eventtype, chrom, strand, totjunc, gene=''):
        self.name = name
        self.eventtype = eventtype
        self.chrom = chrom
        self.strand = strand
        self.gene = gene
        self.totjunc = totjunc
        self.totoverlap = None
        self.rows = {}
        self.other = Counter()
        # all the introns of its alternatives, those of other too, for its gene
        self.introns = set()
        # the other annotated genes sharing an intron or splice site with it
        self.additional_genes = []

    def add(self, alternative, counts, inc_juncs, exc_juncs, outer_juncs=(), exons=(), blocks=()):
        self.rows[alternative] = EventRow(counts, set(inc_juncs), set(exc_juncs), set(outer_juncs), set(exons), list(blocks))

    def keep_other(self, other, min_reads):
        "count the reads of the event's other junctions or sites, if some sample has min_reads of them"
        if any(n >= min_reads for n in other.values()):
            self.other = other


def row_introns(event, alternative):
    """the introns spliced by the reads of an event row: those it includes, or the
    skipping junctions of an exon skip's exclusion"""
    row = event.rows[alternative]
    if event.eventtype in ('es', 'ces') and alternative == 'exc':
        return set(row.outer_juncs)
    return {j for j in row.inc_juncs if isinstance(j, tuple)}


class EventCaller:
    """What finding the events of a strand of a region needs: its chromosome, strand,
    samples, Thresholds and ReadCoverage, with the tests of an event's alternatives"""
    def __init__(self, chrom, strand, samples, thresholds, coverage):
        self.chrom = chrom
        self.strand = strand
        self.samples = samples
        self.thresholds = thresholds
        self.coverage = coverage

    def supported(self, counts, min_reads=None):
        "do some sample's reads reach min_reads, by default feature_support"
        if min_reads is None:
            min_reads = self.thresholds.feature_support
        return any(counts[s] >= min_reads for s in self.samples)

    def varies(self, counts, totals):
        "does the PSI of counts, of totals, vary over the samples with event_support reads (psi_varies)"
        return psi_varies(sample_psis(counts, totals, self.samples, self.thresholds.event_support),
                          self.thresholds.junc_frac_of_event)

    def some_psi_between(self, counts, totals):
        "is counts, of totals, between junc_frac_of_event and 1 minus it in some sample"
        frac = self.thresholds.junc_frac_of_event
        return any(frac <= counts[s] / totals[s] <= 1 - frac for s in self.samples if totals[s] > 0)

    def split(self, alternatives, totals):
        """the supported alternatives, key -> reads per sample, of an event: those that
        vary, which make it an event, and those that don't but are a large part of it
        in some sample (MIN_STEADY_FRACTION), rows too, and the reads of the rest, its
        other row"""
        varying, steady, other = [], [], Counter()
        for key, counts in alternatives.items():
            if self.supported(counts) and self.varies(counts, totals):
                varying.append(key)
            elif self.supported(counts) and any(counts[s] >= MIN_STEADY_FRACTION * totals[s] > 0 for s in self.samples):
                steady.append(key)
            else:
                other += counts
        return varying, steady, other

    def event(self, eventtype, start, end, totals):
        return SplicingEvent(event_name(eventtype, self.chrom, start, end, self.strand), eventtype, self.chrom,
                             self.strand, totals)

    def keep_other(self, event, other):
        event.keep_other(other, self.thresholds.feature_support)

    def covering(self, ref):
        "the reads per sample covering ref"
        return self.coverage.covering(ref)


TERMINAL_JUNCTION_DIST_FROM_ANNOT = 20


# the most of a first or last exon's junction reads, pooled over the samples, that
# have an exon before (or after) it: truncated reads can't add an exon, so a real
# first or last exon has none, while an internal exon whose reads in one sample are
# all truncated has them in the others
MAX_NONTERMINAL_FRACTION = 0.05


# for a first or last exon of --annot_basic, which may also be used as an internal one
MAX_NONTERMINAL_FRACTION_ANNOTATED = 0.25


MIN_TERMINAL_JUNCTION_SEPARATION = 100


MIN_TERMINAL_SS_FRACTION = 0.1


# read ends within this of the next are one start or end
END_GROUP_WINDOW = 100


# a cluster's splicing, reads per sample: of each junction, of each 3' splice site by
# its 5' site and each 5' site by its 3' site, of each pair of consecutive junctions,
# of each block, an internal exon or a read's first or last, and of each skipping
# junction by the junctions of the reads not skipping between its ends
ClusterSplicing = namedtuple('ClusterSplicing', ('juncs', 'ss5to3', 'ss3to5', 'junc_pairs', 'blocks', 'skipped'))


def extract_splicing_info(samples, cluster_chains, strand):
    "the ClusterSplicing of a cluster's chains, junctions -> [ReadEnds] per sample"
    info = ClusterSplicing({}, {}, {}, {}, {}, {})
    for sample, chains in zip(samples, cluster_chains):
        for juncs, ends in chains.items():
            numreads = len(ends)
            for junc in juncs:
                ss5, ss3 = junc if strand == '+' else junc[::-1]
                add_reads(info.ss5to3.setdefault(ss5, {}), ss3, sample, numreads)
                add_reads(info.ss3to5.setdefault(ss3, {}), ss5, sample, numreads)
                add_reads(info.juncs, junc, sample, numreads)
            for prev, nxt in zip(juncs, juncs[1:]):
                add_reads(info.blocks, (prev[1], nxt[0]), sample, numreads)
                add_reads(info.junc_pairs, (prev, nxt), sample, numreads)
            for read in ends:
                add_reads(info.blocks, (read.start, juncs[0][0]), sample, 1)
                add_reads(info.blocks, (juncs[-1][1], read.end), sample, 1)
    # once all the junctions are known: a chain counted before the skipping junction
    # was seen, as in a sample listed before the samples with skipping, was left out
    for sample, chains in zip(samples, cluster_chains):
        for juncs, ends in chains.items():
            add_skipped(info, juncs, sample, len(ends))
    return info


def add_skipped(info, juncs, sample, numreads):
    """count a chain's reads for each skipping junction from one of its junctions' start
    to another's end, by its junctions between them"""
    starts, ends = {j[0] for j in juncs}, {j[1] for j in juncs}
    for outer in info.juncs:
        if outer not in juncs and outer[0] in starts and outer[1] in ends:
            inner = tuple(j for j in juncs if outer[0] <= j[0] and j[1] <= outer[1])
            add_reads(info.skipped.setdefault(outer, {}), inner, sample, numreads)


def _flanking(outer, exon):
    "the junctions each side of exon inside the skipping junction outer"
    return (outer[0], exon[0]), (exon[1], outer[1])


def exon_skipping_events(caller, info):
    """Skipped exons (es): each internal exon, included by its flanking junctions and
    excluded by the junctions joining their outer ends.  Also the reads of each
    junction flanking an included exon by its skipping junction, which alt5 and alt3
    events count apart, and each event's junctions, which ces events don't repeat"""
    outers, inc, exc = {}, {}, {}
    for (prev, nxt), counts in info.junc_pairs.items():
        outer, exon = (prev[0], nxt[1]), (prev[1], nxt[0])
        outers.setdefault(exon, set()).add(outer)
        add_counts(inc, exon, counts)
        add_counts(exc, exon, info.juncs.get(outer, {}))
    events, es_juncs, seen = [], {}, set()
    for exon, exon_outers in outers.items():
        tot = inc[exon] + exc[exon]
        if not (caller.supported(inc[exon]) and caller.supported(exc[exon]) and caller.varies(inc[exon], tot)):
            continue
        # the skipping junctions whose inclusion is a large enough part of the event
        good = {outer for outer in exon_outers if caller.supported(info.junc_pairs[_flanking(outer, exon)])
                and caller.some_psi_between(info.junc_pairs[_flanking(outer, exon)], tot)}
        flanking = {j for outer in good for j in _flanking(outer, exon)}
        event = caller.event('es', *exon, tot)
        events.append(event)
        event.introns = {j for outer in exon_outers for j in _flanking(outer, exon)} | (exon_outers & info.juncs.keys())
        event.add('inc', inc[exon], flanking, (), good, [exon], [junction_blocks(outer, [exon]) for outer in good])
        event.add('exc', exc[exon], (), flanking, good, [exon])
        for outer in good:
            for j in _flanking(outer, exon):
                add_counts(es_juncs.setdefault(j, {}), outer, info.junc_pairs[_flanking(outer, exon)])
        ref = good or exon_outers
        event.totoverlap = caller.covering((max(j[0] for j in ref), min(j[1] for j in ref)))
        seen.add(frozenset(good | flanking))
    return events, es_juncs, seen


def complex_skipping_events(caller, info, es_juncs, seen):
    """Complex exon skipping (ces): each skipping junction against every way the
    reads not skipping include the exons between its ends, by their junctions: several
    exons, or one exon with alternative splice sites, which no other event puts
    together.  Not when an es event has the same junctions.  Adds the reads of the
    junctions flanking the included exons to es_juncs"""
    events = []
    for outer, inner_counts in info.skipped.items():
        exc = info.juncs[outer]
        tot = total([exc, *inner_counts.values()])
        if not (caller.supported(exc) and caller.varies(exc, tot)):
            continue
        varying, steady, other = caller.split(inner_counts, tot)
        if not varying or frozenset(set().union(*varying) | {outer}) in seen:
            continue
        rows = varying + steady
        inner_juncs = set().union(*rows)
        event = caller.event('ces', *outer, tot)
        events.append(event)
        event.introns = {outer}.union(*inner_counts)
        exons = {exon for inner in rows for exon in internal_exons(inner)}
        event.add('exc', exc, (), inner_juncs, [outer], exons, [junction_blocks(outer)])
        for inner in rows:
            inner_exons = internal_exons(inner)
            event.add(coord_ranges(inner_exons), inner_counts[inner], inner, inner_juncs - set(inner), [outer], inner_exons,
                      [junction_blocks(outer, inner_exons)])
            for j in inner:
                add_counts(es_juncs.setdefault(j, {}), outer, inner_counts[inner])
        caller.keep_other(event, other)
        event.totoverlap = caller.covering(outer)
    return events


def _apart_from_skipping(caller, site, alt, alternatives, es_juncs):
    """is the junction of site and alt its own alternative: not flanking a skipped exon,
    or with reads apart from the skipping, over 10% of them and feature_support, in some
    sample"""
    junc = junc_between(site, alt)
    if junc not in es_juncs:
        return True
    skipping = total([es_juncs[junc][outer] for outer in (junc_between(site, other) for other in alternatives if other != alt)
                      if outer in es_juncs[junc]])
    counts = alternatives[alt]
    return any(counts[s] > 0 and (counts[s] - skipping[s]) / counts[s] > 0.1
               and counts[s] - skipping[s] > caller.thresholds.feature_support for s in caller.samples)


def alternative_site_events(caller, eventtype, site_alternatives, es_juncs, terminal_juncs):
    """Splice sites alternative to each other, spliced to one site (alt5 or alt3, by
    eventtype), from site_alternatives, each site's reads by the site it is spliced
    to.  A junction flanking a skipped exon counts only with reads apart from the
    skipping (es_juncs), and the junctions of a first or last exon event,
    terminal_juncs, aren't repeated.  Also each group of alternative sites, by the
    sites spliced to all of them, with their reads, for shared_site_events"""
    events, groups = [], {}
    for site, alternatives in site_alternatives.items():
        if len(alternatives) <= 1:
            continue
        tot = total(alternatives.values())
        varying, steady, other = caller.split(alternatives, tot)
        if frozenset(junc_between(site, alt) for alt in varying) in terminal_juncs:
            continue
        if len(varying) > 1:
            groups.setdefault(frozenset(varying), {})[site] = tot
        # only if alternative splicing, not counting the reads of junctions skipping exons
        own = [alt for alt in varying if _apart_from_skipping(caller, site, alt, alternatives, es_juncs)]
        if len(own) <= 1 or not caller.supported(tot, caller.thresholds.event_support):
            continue
        rows = varying + steady
        juncs = {junc_between(site, alt) for alt in rows}
        # from the shared site to the furthest alternative one, so unique to the shared site
        start, end = span(juncs)
        event = caller.event(eventtype, start, end, tot)
        events.append(event)
        event.introns = {junc_between(site, alt) for alt in alternatives}
        for alt in rows:
            junc = junc_between(site, alt)
            event.add(coord_ranges([junc]), alternatives[alt], [junc], juncs - {junc}, blocks=[[junc]])
        caller.keep_other(event, other)
        # all the alternative sites lie on the same side of the shared one
        event.totoverlap = caller.covering((site, min(rows)) if site < min(rows) else (max(rows), site))
    return events, groups


def shared_site_events(caller, eventtype, groups):
    """Splice sites spliced alike to the same alternative sites (alt5ss or alt3ss, by
    eventtype), from the groups of alternative_site_events"""
    events = []
    for alts, site_totals in groups.items():
        tot = total(site_totals.values())
        if len(site_totals) <= 1 or not caller.supported(tot, caller.thresholds.event_support):
            continue
        varying, steady, other = caller.split(site_totals, tot)
        if len(varying) <= 1:
            continue
        rows = varying + steady
        juncs = {junc_between(site, alt) for site in rows for alt in alts}
        # unique: each site is in one group, and the event ends at its outermost site
        start, end = span(juncs)
        event = caller.event(eventtype, start, end, tot)
        events.append(event)
        event.introns = {junc_between(site, alt) for site in site_totals for alt in alts}
        for site in rows:
            inc = {junc_between(site, alt) for alt in alts}
            # 10 bases of the exon, on its side of the site: upstream of a 5' site, downstream of a 3' one
            exon_left = lies_left(eventtype == 'alt5ss', caller.strand)
            block = (site - 10, site) if exon_left else (site, site + 10)
            event.add(coord_ranges(inc), site_totals[site], inc, juncs - inc, blocks=[[block]])
        caller.keep_other(event, other)
        event.totoverlap = caller.covering((max(rows), min(alts)) if median(rows) < median(alts) else (max(alts), min(rows)))
    return events


def intron_retention_events(caller, info):
    "Retained introns (ir): each junction's reads, and those of the blocks containing it"
    events = []
    for junc, spliced in info.juncs.items():
        retained = total([counts for block, counts in info.blocks.items() if block[0] < junc[0] and junc[1] < block[1]])
        tot = spliced + retained
        if caller.supported(spliced) and caller.supported(retained) and caller.varies(retained, tot):
            event = caller.event('ir', *junc, tot)
            events.append(event)
            event.introns = {junc}
            event.add('spliced', spliced, [junc], (), blocks=[junction_blocks(junc)])
            event.add('retained', retained, (), [junc], blocks=[[junc]])
            event.totoverlap = caller.covering(junc)
    return events


def extract_terminal_info(samples, chains_by_sample, key, strand, terminal):
    """the junctions, first (afe) or last (ale) exons and read starts (or ends) of the
    chains of key, a 5' (or 3') gene: reads per sample of each junction, of each first
    (or last) exon's junction by its splice site, and (position, sample) per read.  Its
    fragments, under no junctions, have only their starts (or ends)"""
    juncs, terminal_exons, ends = {}, {}, []
    # a first exon's splice site, and the read's start, are upstream, a last exon's downstream
    take_start = lies_left(terminal == 'afe', strand)
    for sample, chains in zip(samples, chains_by_sample):
        for chain_juncs, reads in chains.get(key, {}).items():
            for junc in chain_juncs:
                add_reads(juncs, junc, sample, len(reads))
            if chain_juncs:
                junc, site = (chain_juncs[0], chain_juncs[0][0]) if take_start else (chain_juncs[-1], chain_juncs[-1][1])
                add_reads(terminal_exons.setdefault(site, {}), junc, sample, len(reads))
            ends.extend((read.start if take_start else read.end, sample) for read in reads)
    return juncs, terminal_exons, ends


def group_ends(ends, window):
    """read ends, (position, sample), clustered where they lie within window of the
    next, each cluster's reads per sample counted at its median"""
    grouped, cluster = {}, []
    for end in sorted(ends):
        if cluster and end[0] - cluster[-1][0] > window:
            _add_end_cluster(grouped, cluster)
            cluster = []
        cluster.append(end)
    if cluster:
        _add_end_cluster(grouped, cluster)
    return grouped


def _add_end_cluster(grouped, cluster):
    grouped.setdefault(int(median([pos for pos, _ in cluster])), Counter()).update([sample for _, sample in cluster])


def near_annotated_site(site, annot_sites):
    "is site within TERMINAL_JUNCTION_DIST_FROM_ANNOT of one of annot_sites, a sorted list"
    i = bisect_left(annot_sites, site - TERMINAL_JUNCTION_DIST_FROM_ANNOT)
    return i < len(annot_sites) and annot_sites[i] <= site + TERMINAL_JUNCTION_DIST_FROM_ANNOT


def is_terminal_site(junc_counts, juncs, min_read_support, annotated):
    """Is a splice site of first (or last) junctions really a first (or last) exon's:
    pooled over the samples, at least min_read_support reads have its junctions as
    their first (or last), and few of its junctions' reads have an exon before (or
    after) it, more being allowed of an annotated one"""
    terminal = sum(sum(counts.values()) for counts in junc_counts.values())
    all_reads = sum(sum(juncs[junc].values()) for junc in junc_counts)
    max_fraction = MAX_NONTERMINAL_FRACTION_ANNOTATED if annotated else MAX_NONTERMINAL_FRACTION
    return terminal >= min_read_support and all_reads - terminal <= max_fraction * all_reads


def terminal_exon_events(caller, eventtype, terminal_exons, juncs, annot_sites):
    """Alternative first (afe) or last (ale) exons of a gene, from terminal_exons, its
    reads' first or last junctions by their splice site.  Only the sites really of a
    first or last exon (is_terminal_site) are its alternatives, and only their reads
    are the event's: a truncated read's first junction, after an internal exon, is
    neither.  Also the event's junctions, which alt5 and alt3 events don't repeat"""
    site_counts = {site: total(junc_counts.values()) for site, junc_counts in terminal_exons.items()
                   if is_terminal_site(junc_counts, juncs, caller.thresholds.feature_support,
                                       near_annotated_site(site, annot_sites))}
    tot = total(site_counts.values())
    if len(site_counts) <= 1 or not caller.supported(tot, caller.thresholds.event_support):
        return [], set()
    varying, steady, other = caller.split(site_counts, tot)
    # the sites of a large enough part of the reads, far enough apart
    major = [site for site in varying
             if any(site_counts[site][s] / tot[s] >= MIN_TERMINAL_SS_FRACTION for s in caller.samples if tot[s] > 0)]
    if len(major) <= 1 or max(major) - min(major) < MIN_TERMINAL_JUNCTION_SEPARATION:
        return [], set()

    def alternative_juncs(sites):
        return {junc for site in sites for junc, counts in terminal_exons[site].items()
                if caller.supported(counts) and caller.some_psi_between(counts, tot)}

    rows = varying + steady
    event_juncs = alternative_juncs(rows)
    # each first or last exon by its introns: those of the event's, or failing that, all its own
    site_juncs = {site: ({j for j in event_juncs if site in j} or set(terminal_exons[site])) for site in rows}
    event = caller.event(eventtype, *span([j for js in site_juncs.values() for j in js]), tot)
    for site in rows:
        # 100 bases of the exon, on its side of the splice site: upstream of a first exon's, downstream of a last's
        exon_left = lies_left(eventtype == 'afe', caller.strand)
        block = (site - 100, site) if exon_left else (site, site + 100)
        inc = {j for j in event_juncs if site in j}
        event.add(coord_ranges(site_juncs[site]), site_counts[site], inc, event_juncs - inc, blocks=[[block]])
    caller.keep_other(event, other)
    event.totoverlap = caller.covering((max(j[0] for j in event_juncs), min(j[1] for j in event_juncs)))
    # the varying sites' junctions, as an alt5 or alt3 event's are compared
    return [event], {frozenset(alternative_juncs(varying))}


def end_events(caller, eventtype, ends):
    "Alternative starts (tss) or ends (tts) of a gene, from its reads' ends grouped by group_ends"
    tot = total(ends.values())
    if not caller.supported(tot, caller.thresholds.event_support):
        return []
    varying, steady, other = caller.split(ends, tot)
    if len(varying) <= 1:
        return []
    # the 1-based position of each end's base: the reads' starts are 0-based, on +
    # their tss and on - their tts
    from_starts = lies_left(eventtype == 'tss', caller.strand)
    positions = {end: end + 1 if from_starts else end for end in varying + steady}
    event = caller.event(eventtype, min(positions.values()) - 1, max(positions.values()), tot)
    # all the reads of a start or end are the gene's, so its locus is the event
    event.totoverlap = tot
    for end, pos in positions.items():
        event.add(str(pos), ends[end], {pos}, set(positions.values()) - {pos}, blocks=[[(pos - 1, pos)]])
    caller.keep_other(event, other)
    return [event]
