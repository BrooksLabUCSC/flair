"""The regions reads are corrected and events reported in: chromosomes cut in
intergenic gaps, regions extended over their reads, and the events each region
reports.  Used by flair spliceevents"""

from bisect import bisect_right
import pysam
from flair import SeqRange


def region_key(region):
    return (region.name, region.start, region.end)


def gene_clusters(region, transcripts):
    """The spans of region's annotated genes, both strands, merged where they overlap,
    and the number of exons in each, sorted"""
    spans = sorted((max(t.start, region.start), min(t.end, region.end), len(t.exons)) for t in transcripts
                   if t.chrom == region.name and t.start < region.end and region.start < t.end)
    clusters = []
    for start, end, exons in spans:
        if clusters and start < clusters[-1][1]:
            clusters[-1][1] = max(clusters[-1][1], end)
            clusters[-1][2] += exons
        else:
            clusters.append([start, end, exons])
    return clusters


def intergenic_cuts(region, transcripts, min_exons=0):
    """Where to cut region into partitions: the middle of the gaps between its
    gene clusters, each partition with at least min_exons annotated exons, but
    the last, so 0 cuts at every gap"""
    clusters = gene_clusters(region, transcripts)
    cuts, exons = [], 0
    for (_, end, n), (start, _, _) in zip(clusters, clusters[1:]):
        exons += n
        if exons >= min_exons:
            cuts.append((end + start) // 2)
            exons = 0
    return cuts


def _read_spans(bamfiles, region):
    "the (start, end) of each mapped, not secondary, alignment overlapping region, in all bam files"
    for bamfile in bamfiles:
        with pysam.AlignmentFile(bamfile, 'rb') as bam:
            for read in bam.fetch(region.name, region.start, region.end):
                if read.is_mapped and not read.is_secondary:
                    yield read.reference_start, read.reference_end


def correction_regions(report_regions, bamfiles):
    """The regions to correct reads in, for the regions events are reported in: each
    report region extended over the reads overlapping it in any sample, so a read
    running past its end is corrected with the annotation and junctions it reaches,
    as anywhere else.  Extended regions that overlap are merged.  Returns the
    regions, and a map of each one's key to the report regions in it."""
    extended = []
    for region in report_regions:
        start, end = region.start, region.end
        for read_start, read_end in _read_spans(bamfiles, region):
            start, end = min(start, read_start), max(end, read_end)
        extended.append((SeqRange(region.name, start, end), region))
    merged = combine_regions([ext for ext, _ in extended])
    reports = {region_key(region): [] for region in merged}
    for ext, region in extended:
        container = next(m for m in merged if m.name == ext.name and m.start <= ext.start and ext.end <= m.end)
        reports[region_key(container)].append(region)
    return merged, reports


def partition_regions(report_regions, bamfiles, cuts):
    """The partitions of report regions cut at cuts, a dict of region key to sorted
    positions: each piece owns the reads aligned to start in it, the first also
    those starting before the region, and is corrected over their whole length.
    Reads crossing a cut are only the piece's they start in, so the events past a
    cut lose them.  Returns, for each piece with reads, its region to correct reads
    in, the range its reads start in, its report region, and its number of reads"""
    pieces = []
    for region in report_regions:
        bounds = [region.start] + cuts.get(region_key(region), []) + [region.end]
        low, high = bounds[:-1], bounds[1:]
        reads = [0] * len(low)
        for read_start, read_end in _read_spans(bamfiles, region):
            i = max(bisect_right(bounds, read_start) - 1, 0)
            low[i], high[i] = min(low[i], read_start), max(high[i], read_end)
            reads[i] += 1
        for i in range(len(low)):
            if reads[i] > 0:
                owned = SeqRange(region.name, bounds[i] if i > 0 else 0, bounds[i + 1])
                pieces.append((SeqRange(region.name, low[i], high[i]), owned, region, reads[i]))
    return pieces


# events whose reads are all of a gene's: its starts and ends, first and last exons
GENE_WIDE_EVENTS = {'tss', 'tts', 'afe', 'ale'}


def _event_positions(event):
    "the splice sites and ends of an event's junctions, exons and ends"
    positions = []
    for row in event.rows.values():
        for feature in (*row.inc_juncs, *row.exc_juncs, *row.outer_juncs, *row.exons):
            positions.extend(feature if isinstance(feature, tuple) else (feature,))
    return positions


def gene_junction_span(juncs):
    """The span of a gene's junctions, its annotated ones, or a novel gene's of its
    reads, from the first one's start to the last one's end"""
    return (min(j[0] for j in juncs), max(j[1] for j in juncs))


def event_in_regions(event, regions, gene_span, owned=None):
    """Is an event inside one of the regions events are reported in: all its
    junctions and exons, so every read reaching it overlaps the region and its
    counts are complete.  An event of all its gene's reads (GENE_WIDE_EVENTS) only
    when all the gene's annotated junctions are, gene_span (gene_junction_span),
    as a user choosing a region around a gene would think of it: its reads are
    spliced, so nearly all have a junction in the region and were fetched.  The
    gene's first and last exons may run past the region.  Given owned, the range
    a partition's reads start in, only an event starting in it, so each event is
    reported by one partition"""
    if event.eventtype in GENE_WIDE_EVENTS:
        low, high = gene_span
    else:
        positions = _event_positions(event)
        low, high = min(positions), max(positions)
    if owned is not None and not owned.start <= low < owned.end:
        return False
    return any(region.name == event.chrom and region.start <= low and high <= region.end for region in regions)


def combine_regions(regions):
    "regions sorted, with those overlapping or touching merged"
    merged = []
    for region in sorted(regions):
        if merged and merged[-1].name == region.name and region.start <= merged[-1].end:
            merged[-1] = SeqRange(region.name, merged[-1].start, max(merged[-1].end, region.end))
        else:
            merged.append(SeqRange(region.name, region.start, region.end))
    return merged
