"""The tables of splicing events: their columns, the rows each region writes, and
their combination into a run's.  Used by flair spliceevents"""

from flair.partition_runner import combine_temp_files_by_suffix
from flair.io_utils import sort_bed_file
from statistics import median
from flair.pycbio.hgdata.bed import Bed, BedBlock
from flair.spliceevents_events import EVENT_COLORS, row_introns, row_name
from flair.splicing_events_tsv import SplicingEventsWriter
from flair.splicing_totals_tsv import SplicingTotalsWriter


def junc_list(chrom, features):
    "junctions or exons, 1-based, or end positions, each with its chromosome"
    return [f'{chrom}:{f[0] + 1}-{f[1]}' if isinstance(f, tuple) else f'{chrom}:{f}' for f in sorted(features)]


def _names(field):
    "a comma separated field of names as a list"
    return field.split(',') if field else []


def is_expressed(event, samples, thresholds):
    """Is an event expressed enough to write: in some sample, its reads at least
    event_frac_of_tot of its locus reads"""
    locus = event.totoverlap
    return any(locus[s] > 0 and thresholds.event_frac_of_tot <= event.totjunc[s] / locus[s] for s in samples)


class EventTables:
    """A region's tables of events, open for writing: each row's counts and two PSIs,
    over the event's reads (PSIjunc) and over the reads covering its locus (PSItot),
    with the row's introns that aren't annotated, its BED lines scored by its median
    PSIjunc, and the event's two totals, which spliceevents_outliers divides by"""
    def __init__(self, prefix, samples, event_support):
        self.samples = samples
        self.event_support = event_support
        self.bed = open(prefix + '.spliceevents.bed', 'w')
        self.tables = {suffix: SplicingEventsWriter(prefix + suffix, samples) for suffix in EVENT_TABLE_SUFFIXES}
        self.totals = SplicingTotalsWriter(prefix + TOTALS_SUFFIX, samples)

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.bed.close()
        for writer in [*self.tables.values(), self.totals]:
            writer.close()

    def _write_row(self, columns, counts, event):
        "a row's counts and PSIs; returns its PSIjuncs, of the samples with event_support reads"
        self.tables['.spliceevents.counts.tsv'].writeRow(columns + [str(counts[s]) for s in self.samples])
        for suffix, totals in (('.spliceevents.PSIjunc.tsv', event.totjunc), ('.spliceevents.PSItot.tsv', event.totoverlap)):
            self.tables[suffix].writeRow(columns + [str(round(counts[s] / totals[s], 4)) if totals[s] >= self.event_support
                                                    else '' for s in self.samples])
        return [counts[s] / event.totjunc[s] for s in self.samples if event.totjunc[s] >= self.event_support]

    def _write_beds(self, event, alternative, blocks_list, psis):
        "a row's BED records, scored by its median PSIjunc"
        for blocks in blocks_list:
            start, end = blocks[0][0], blocks[-1][1]
            Bed(event.chrom, start, end, row_name(event.name, alternative), score=round(median(psis) * 100),
                strand=event.strand, thickStart=start, thickEnd=end, itemRgb=EVENT_COLORS[event.eventtype],
                blocks=[BedBlock(s, e) for s, e in blocks]).write(self.bed)

    def write(self, event, annotated):
        genes = [_names(event.gene), _names(annotated.gene_names(event.gene.split(','))),
                 list(event.additional_genes), _names(annotated.gene_names(event.additional_genes))]
        for alternative, row in event.rows.items():
            novel = {j for j in row_introns(event, alternative) if not annotated.is_annotated(j, event.strand)}
            columns = [event.name, alternative, event.eventtype, *genes,
                       *(junc_list(event.chrom, features)
                         for features in (row.inc_juncs, row.exc_juncs, row.outer_juncs, row.exons, novel))]
            self._write_beds(event, alternative, row.blocks, self._write_row(columns, row.counts, event))
        if any(event.other[s] > 0 for s in self.samples):
            self._write_row([event.name, 'other', event.eventtype, *genes, [], [], [], [], []], event.other, event)
        for label, counts in (('event_reads', event.totjunc), ('locus_reads', event.totoverlap)):
            self.totals.writeRow([event.name, event.eventtype, *genes, label] + [counts[s] for s in self.samples])


# the tables each region writes, combined into the run's: a file of events
# (splicing_events_tsv) per kind of sample column, their totals (splicing_totals_tsv)
# and their BED
EVENT_TABLE_SUFFIXES = ('.spliceevents.counts.tsv', '.spliceevents.PSIjunc.tsv', '.spliceevents.PSItot.tsv')
TOTALS_SUFFIX = '.spliceevents.totals.tsv'
REGION_OUTPUT_SUFFIXES = ('.spliceevents.bed', *EVENT_TABLE_SUFFIXES, TOTALS_SUFFIX)


def combine_region_tables(output, prefixes):
    """the regions' tables as the run's, OUTPUT.spliceevents.*, the BED sorted by position
    as every final BED is.  Each region's table parts carry the header, kept once"""
    combine_temp_files_by_suffix(output, prefixes, ['.spliceevents.bed'])
    sort_bed_file(output + '.spliceevents.bed')
    tables = [*EVENT_TABLE_SUFFIXES, TOTALS_SUFFIX]
    combine_temp_files_by_suffix(output, prefixes, tables, headers=tables)
