"""flair spliceevents with --region_bed: which events are reported"""
from flair import SeqRange
from flair.flair_spliceevents import SplicingEvent, SplicingEventJunction, event_in_regions, gene_junction_span
from flair.isoform_data import Junc

REGION = [SeqRange('chr1', 1000, 5000)]


def _event(eventtype, *juncs, chrom='chr1'):
    event = SplicingEvent('e', eventtype, 'G1', chrom, '+', {'s': 0}, ['s'])
    event.events['inc'] = SplicingEventJunction('inc', [], {'s': 0}, set(juncs), set())
    return event


def test_an_event_inside_the_region_is_reported():
    assert event_in_regions(_event('es', (1200, 1500), (1600, 2000)), REGION, (0, 9000))


def test_an_event_reaching_past_the_region_is_not():
    # its reads past 5000 weren't all fetched
    assert not event_in_regions(_event('es', (1200, 1500), (1600, 5200)), REGION, (1100, 5500))


def test_an_event_on_another_chromosome_is_not():
    assert not event_in_regions(_event('es', (1200, 1500), chrom='chr2'), REGION, (0, 9000))


def test_a_gene_wide_event_needs_the_whole_gene_in_the_region():
    # a TSS event's PSIs are of all the gene's reads
    tss = _event('tss', 1200, 1300)
    assert event_in_regions(tss, REGION, (1100, 4900))
    assert not event_in_regions(tss, REGION, (1100, 6000))


def test_a_gene_wide_event_is_reported_when_the_genes_first_and_last_exons_run_past_the_region():
    # junctions 1500-2000 and 3000-4000 in the region, the gene's exons 900-1500 and 4000-5200 past it:
    # every spliced read of the gene has a junction in the region, so was fetched
    span = gene_junction_span({Junc(1500, 2000), Junc(3000, 4000)})
    assert span == (1500, 4000)
    assert event_in_regions(_event('tss', 1200, 1300), REGION, span)


def test_a_gene_wide_event_is_not_reported_with_an_annotated_junction_past_the_region():
    span = gene_junction_span({Junc(1500, 2000), Junc(4800, 5300)})
    assert span == (1500, 5300)
    assert not event_in_regions(_event('tss', 1200, 1300), REGION, span)
