"""The regions reads are corrected and events reported in: which reads of a
partition to fetch, and which events it reports.  Used by flair spliceevents"""

from flair.read_processing import ReadWindows


def region_key(region):
    return (region.name, region.start, region.end)


def read_windows(reports):
    """The ReadWindows of the reads of a partition to correct: those overlapping one
    of its report regions.  A partition can be far wider than the report regions in
    it, as it is a cluster of overlapping reads and annotated transcripts; fetching
    the windows alone never reads the rest"""
    return ReadWindows(sorted(reports), False)


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


def event_in_regions(event, regions, gene_span):
    """Is an event inside one of the regions events are reported in: all its
    junctions and exons, so every read reaching it overlaps the region and its
    counts are complete.  An event of all its gene's reads (GENE_WIDE_EVENTS) only
    when all the gene's annotated junctions are, gene_span (gene_junction_span),
    as a user choosing a region around a gene would think of it: its reads are
    spliced, so nearly all have a junction in the region and were fetched.  The
    gene's first and last exons may run past the region.  No read crosses a
    partition, so an event is found by the one partition holding its reads"""
    if event.eventtype in GENE_WIDE_EVENTS:
        low, high = gene_span
    else:
        positions = _event_positions(event)
        low, high = min(positions), max(positions)
    return any(region.name == event.chrom and region.start <= low and high <= region.end for region in regions)
