"""The tables of splicing events: their columns, the rows each region writes, and
their combination into a run's.  Used by flair spliceevents"""

import shutil
from flair.partition_runner import combine_temp_files_by_suffix
from statistics import median
from flair.pycbio.hgdata.bed import Bed, BedBlock
from flair.spliceevents_events import EVENT_COLORS, psi_varies, row_introns, row_name, sample_psis


def junc_string(chrom, features):
    "junctions or exons, 1-based, or end positions, each with its chromosome"
    return ','.join(f'{chrom}:{f[0] + 1}-{f[1]}' if isinstance(f, tuple) else f'{chrom}:{f}' for f in sorted(features))


def is_expressed(event, samples, thresholds):
    """Is an event expressed enough to write: in some sample, its reads at least
    event_frac_of_tot of its locus reads, and some row, but an exon skip's exclusion,
    varying as a part of the locus reads"""
    locus = event.totoverlap
    if not any(locus[s] > 0 and thresholds.event_frac_of_tot <= event.totjunc[s] / locus[s] for s in samples):
        return False
    return any(psi_varies(sample_psis(row.counts, locus, samples, 1), thresholds.junc_frac_of_event)
               for alternative, row in event.rows.items() if event.eventtype not in ('es', 'ces') or alternative != 'exc')


class EventTables:
    """A region's tables of events, open for writing: each row's counts and two PSIs,
    over the event's reads (PSIjunc) and over the reads covering its locus (PSItot),
    with the row's introns that aren't annotated, its BED lines scored by its median
    PSIjunc, and the event's two totals, which spliceevents_outliers divides by"""
    def __init__(self, prefix, samples, event_support):
        self.samples = samples
        self.event_support = event_support
        self.files = {suffix: open(prefix + suffix, 'w') for suffix in REGION_OUTPUT_SUFFIXES}

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        for fh in self.files.values():
            fh.close()

    def _line(self, suffix, columns):
        self.files[suffix].write('\t'.join(columns) + '\n')

    def _write_row(self, columns, counts, event):
        "a row's counts and PSIs; returns its PSIjuncs, of the samples with event_support reads"
        self._line('.spliceevents.counts.tsv', columns + [str(counts[s]) for s in self.samples])
        for suffix, totals in (('.spliceevents.PSIjunc.tsv', event.totjunc), ('.spliceevents.PSItot.tsv', event.totoverlap)):
            self._line(suffix, columns + [str(round(counts[s] / totals[s], 4)) if totals[s] >= self.event_support else ''
                                          for s in self.samples])
        return [counts[s] / event.totjunc[s] for s in self.samples if event.totjunc[s] >= self.event_support]

    def _write_beds(self, event, alternative, blocks_list, psis):
        "a row's BED records, scored by its median PSIjunc, none if no sample has a PSI"
        if not psis:
            return
        for blocks in blocks_list:
            start, end = blocks[0][0], blocks[-1][1]
            Bed(event.chrom, start, end, row_name(event.name, alternative), score=round(median(psis) * 100),
                strand=event.strand, thickStart=start, thickEnd=end, itemRgb=EVENT_COLORS[event.eventtype],
                blocks=[BedBlock(s, e) for s, e in blocks]).write(self.files['.spliceevents.bed'])

    def write(self, event, annotated):
        genes = [event.gene, annotated.gene_names(event.gene.split(',')),
                 ','.join(event.additional_genes), annotated.gene_names(event.additional_genes)]
        for alternative, row in event.rows.items():
            novel = {j for j in row_introns(event, alternative) if not annotated.is_annotated(j, event.strand)}
            columns = [event.name, alternative, event.eventtype, *genes,
                       *(junc_string(event.chrom, features)
                         for features in (row.inc_juncs, row.exc_juncs, row.outer_juncs, row.exons, novel))]
            self._write_beds(event, alternative, row.blocks, self._write_row(columns, row.counts, event))
        if any(event.other[s] > 0 for s in self.samples):
            self._write_row([event.name, 'other', event.eventtype, *genes, '', '', '', '', ''], event.other, event)
        for label, counts in (('event_reads', event.totjunc), ('locus_reads', event.totoverlap)):
            self._line('.spliceevents.totals.tsv', [event.name, event.eventtype, *genes, label] + [str(counts[s]) for s in self.samples])


# the tables each region writes, combined into the run's
REGION_OUTPUT_SUFFIXES = ('.spliceevents.bed', '.spliceevents.counts.tsv', '.spliceevents.PSIjunc.tsv',
                          '.spliceevents.PSItot.tsv', '.spliceevents.totals.tsv')


# the columns of the counts and PSI tables before the samples': a row of an event,
# its event, and its junctions, and of the totals table
EVENT_GENE_COLUMNS = ['gene', 'gene_name', 'additional_genes', 'additional_gene_names']


EVENT_ROW_COLUMNS = ['event', 'alternative', 'event_type', *EVENT_GENE_COLUMNS, 'junctions_included', 'junctions_excluded',
                     'outer_junctions', 'exons', 'novel_junctions']


EVENT_TOTALS_COLUMNS = ['event', 'event_type', *EVENT_GENE_COLUMNS, 'total']


def write_with_header(output, prefixes, suffix, header):
    "the partitions' tables of one suffix, as OUTPUT + suffix under header"
    with open(output + suffix, 'w') as out:
        out.write('\t'.join(header) + '\n')
        for prefix in prefixes:
            with open(prefix + suffix) as fh:
                shutil.copyfileobj(fh, out)


def combine_region_tables(output, prefixes, samples):
    "the regions' tables as the run's, OUTPUT.spliceevents.*"
    combine_temp_files_by_suffix(output, prefixes, ['.spliceevents.bed'])
    for suffix, columns in (('.spliceevents.counts.tsv', EVENT_ROW_COLUMNS), ('.spliceevents.PSIjunc.tsv', EVENT_ROW_COLUMNS),
                            ('.spliceevents.PSItot.tsv', EVENT_ROW_COLUMNS), ('.spliceevents.totals.tsv', EVENT_TOTALS_COLUMNS)):
        write_with_header(output, prefixes, suffix, columns + samples)
