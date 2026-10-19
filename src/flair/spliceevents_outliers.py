#!/usr/bin/env python3
"""
Outlier samples of flair spliceevents' events: samples whose PSI of an event row is
far from the median of all samples', from OUTPUT.spliceevents.counts.tsv and
OUTPUT.spliceevents.totals.tsv.

NOTE: the outliers of starts and ends (tss and tts) are left out of the filtered
table, OUTPUT.spliceevents.outliers.filtered.tsv, on purpose: their events, from
read ends, give far too many false positive outliers.  They are still in the
table of every outlier.  When start and end events improve, add them to
SIMPLE_EVENTS and take them out of LEFT_OUT_EVENTS.
"""
import argparse
import csv
from collections import defaultdict
from itertools import product
from statistics import median
from flair import FlairInputDataError
from flair.pycbio.sys import cli

# a sample is an outlier when its PSI is at least MIN_DELTA_PSI from the median and
# more than MIN_DEVIATIONS deviations (half the interquartile range) from it, of at
# least MIN_SAMPLES samples with a PSI
MIN_DELTA_PSI = 0.1
MIN_DEVIATIONS = 3
MIN_SAMPLES = 5
# filtered outliers are this many deviations from the median
FILTERED_MIN_DEVIATIONS = 6
# an outlier of a complex exon skip or alternative splice site that a simpler event
# also explains is kept only if its deviation is this much better than that event's
BETTER_DEVIATION_FRACTION = 0.8
# no deviation, the samples' PSIs equal but this one's
NO_DEVIATION = 100000

OUTLIER_COLUMNS = ['event', 'alternative', 'event_type', 'gene', 'gene_name', 'additional_genes', 'additional_gene_names',
                   'sample', 'medianPSI', 'dev(IQR/2)', 'sample_val', 'tot_not_NA_samples', 'event_reads;total_locus_reads', 'delta_PSI_to_med', 'dev_from_med']
# the denominator of each PSI, a row of the totals table
PSI_TOTALS = {'junc': 'event_reads', 'tot': 'locus_reads'}


def build_parser():
    desc = ('outlier samples of flair spliceevents events, whose PSI of an event is far from the '
            "median of all samples'")
    parser = argparse.ArgumentParser(prog='spliceevents_outliers', description=desc)
    parser.add_argument('-i', '--input', required=True,
                        help='output prefix of flair spliceevents: reads its OUTPUT.spliceevents.counts.tsv '
                             'and OUTPUT.spliceevents.totals.tsv')
    parser.add_argument('-o', '--output', required=True,
                        help='output prefix: writes OUTPUT.spliceevents.outliers.tsv, every outlier, and '
                             f'OUTPUT.spliceevents.outliers.filtered.tsv, those more than {FILTERED_MIN_DEVIATIONS} '
                             'deviations from the median and not explained by a simpler event, without starts and '
                             'ends (tss, tts), whose outliers are mostly false positives')
    parser.add_argument('--psi', choices=sorted(PSI_TOTALS), default='junc',
                        help="the PSI to call outliers of: 'junc', of the event's reads (PSIjunc), more sensitive, "
                             "or 'tot', of the reads covering its locus (PSItot), for strong overall events "
                             '(default: %(default)s)')
    parser.add_argument('--event_support', type=int, default=20,
                        help='minimum reads of a sample for its PSI of an event to be used, as flair '
                             'spliceevents --event_support (default: %(default)s)')
    return parser


def read_tsv(path):
    "the header and rows, each a dict, of a tab-separated table"
    with open(path, newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        return reader.fieldnames, list(reader)


def read_events(input_prefix, psi):
    """the event rows of the counts table, in table order, each with its per-sample
    counts and the per-sample totals of its PSI, and the samples"""
    header, rows = read_tsv(input_prefix + '.spliceevents.counts.tsv')
    if 'novel_junctions' not in header:
        raise FlairInputDataError(f'{input_prefix}.spliceevents.counts.tsv has no novel_junctions column: '
                                  'it is from an earlier flair spliceevents, rerun it')
    samples = header[header.index('novel_junctions') + 1:]
    _, total_rows = read_tsv(input_prefix + '.spliceevents.totals.tsv')
    totals = {row['event']: {s: int(row[s]) for s in samples} for row in total_rows if row['total'] == PSI_TOTALS[psi]}
    event_rows = []
    for row in rows:
        if row['alternative'] == 'other':  # the reads of an event's other ends or junctions
            continue
        if row['event'] not in totals:
            raise FlairInputDataError(f"event {row['event']} of {input_prefix}.spliceevents.counts.tsv is not in "
                                      f'{input_prefix}.spliceevents.totals.tsv')
        row['counts'] = {s: int(row[s]) for s in samples}
        row['totals'] = totals[row['event']]
        event_rows.append(row)
    return event_rows, samples


def _junctions(field):
    return frozenset(field.split(',')) if field else frozenset()


def _row_junctions(row):
    "an event row's junctions: (included, excluded, outer)"
    return (_junctions(row['junctions_included']), _junctions(row['junctions_excluded']),
            _junctions(row['outer_junctions']))


def group_by_shared_junctions(rows):
    """event rows grouped by events sharing junctions, linked through each other, in
    table order: an outlier is only explained by an event of the same junctions.  An
    event of no junctions, a start or end, is alone"""
    parent = {}

    def find(x):
        while parent.setdefault(x, x) != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for row in rows:
        event = ('event', row['event'])
        for junc in set().union(*_row_junctions(row)):
            parent[find(('junc', junc))] = find(event)
    groups = defaultdict(list)
    for row in rows:
        groups[find(('event', row['event']))].append(row)
    return list(groups.values())


def row_outliers(row, samples, event_support):
    """the outliers of an event row: per sample with one, its outlier line"""
    import scipy.stats as sps  # 0.7s to load
    psis = {s: row['counts'][s] / row['totals'][s] for s in samples if row['totals'][s] >= event_support}
    if len(psis) < MIN_SAMPLES:
        return {}
    vals = list(psis.values())
    med, dev = median(vals), sps.iqr(vals) / 2
    outliers = {}
    for s, psi in psis.items():
        delta = abs(psi - med)
        if delta >= MIN_DELTA_PSI and ((dev == 0 and psi != med) or delta > dev * MIN_DEVIATIONS):
            outliers[s] = [row['event'], row['alternative'], row['event_type'], row['gene'], row['gene_name'],
                           row['additional_genes'], row['additional_gene_names'], s, round(med, 6), round(dev, 6), round(psi, 6),
                           len(vals), f"{row['counts'][s]};{row['totals'][s]}", round(delta, 6),
                           round(delta / dev, 6) if dev != 0 else NO_DEVIATION]
    return outliers


def _by_deviation(lines):
    return sorted(lines, reverse=True, key=lambda x: x[::-1])


def _better_than(this, other):
    "is this outlier line's deviation enough better than the other's"
    return other[-1] / this[-1] < BETTER_DEVIATION_FRACTION


def _write_lines(fh, lines):
    for line in lines:
        fh.write('\t'.join(str(x) for x in line) + '\n')


# the events whose outliers are all kept, and those whose outliers a simpler event
# may explain, compared after them.  Alternative splice sites overlapping first and
# last exons are already left out when the events are found
SIMPLE_EVENTS = ('es', 'afe', 'ale', 'ir', 'alt5ss', 'alt3ss')
COMPLEX_EVENTS = ('ces', 'alt5', 'alt3')
# NOTE: left out of the filtered outliers on purpose, for now: start and end
# events, from read ends, give far too many false positive outliers.  Move them
# to SIMPLE_EVENTS when they are better (see the module docstring)
LEFT_OUT_EVENTS = ('tss', 'tts')


def _explaining(eventtype, juncs, es_juncs, site_juncs):
    """A complex exon skip's or alternative splice site's own junctions, and the
    outliers of the events that may explain it: for a complex exon skip, the exon
    skips holding its junctions, its exons' for its exclusion row, and for an
    alternative splice site, the events in site_juncs holding both its junctions"""
    inc, exc, outer = juncs
    if eventtype == 'ces':
        own = (inc or exc) | outer
        return own, [outliers for skip_juncs, outliers in es_juncs.items() if own <= skip_juncs]
    own = inc | exc
    return own, [outliers for key, by_juncs in site_juncs.items() if own <= key for outliers in by_juncs.values()]


def filter_outliers(etype_to_sig):
    """The outliers, more than FILTERED_MIN_DEVIATIONS out, of events sharing
    junctions, given by event type, then junctions (included, excluded, outer), then
    sample.  An outlier of a complex exon skip or alternative splice site is kept
    only when no simpler event of its junctions explains it as well: an exon skip,
    or a complex event not explained itself.  Starts and ends (LEFT_OUT_EVENTS)
    are left out, and an event type of none of the lists is an error"""
    unknown = set(etype_to_sig) - set(SIMPLE_EVENTS + COMPLEX_EVENTS + LEFT_OUT_EVENTS)
    if unknown:
        raise ValueError(f'event types {sorted(unknown)} not in SIMPLE_EVENTS, COMPLEX_EVENTS or LEFT_OUT_EVENTS')
    # the exon skips by their junctions, and the events by the junctions an
    # alternative splice site they explain splices: an exon skip's flanking junction
    # and its skipping one
    es_juncs, site_juncs, kept = {}, {}, []
    for eventtype in SIMPLE_EVENTS:
        for juncs, outliers in etype_to_sig.get(eventtype, {}).items():
            if eventtype == 'es':
                inc, exc, outer = juncs
                es_juncs[inc | exc | outer] = outliers
                for flanking, skipping in product(inc, outer):
                    site_juncs.setdefault(frozenset((flanking, skipping)), {})[juncs] = outliers
            kept.extend(outliers.values())
    for eventtype in COMPLEX_EVENTS:
        for juncs, outliers in etype_to_sig.get(eventtype, {}).items():
            own, explaining = _explaining(eventtype, juncs, es_juncs, site_juncs)
            if explaining:
                kept.extend(line for s, line in outliers.items()
                            if any(s not in other or _better_than(line, other[s]) for other in explaining))
            else:
                site_juncs.setdefault(own, {})[juncs] = outliers
                kept.extend(outliers.values())
    return [line for line in _by_deviation(kept) if line[-1] > FILTERED_MIN_DEVIATIONS]


def spliceevents_outliers(input_prefix, output, psi, event_support):
    """Write every outlier of the events of input_prefix, and those filtered by
    filter_outliers, by group of events sharing junctions"""
    event_rows, samples = read_events(input_prefix, psi)
    with open(output + '.spliceevents.outliers.tsv', 'w') as out, \
            open(output + '.spliceevents.outliers.filtered.tsv', 'w') as out_filtered:
        for fh in (out, out_filtered):
            fh.write('\t'.join(OUTLIER_COLUMNS) + '\n')
        for rows in group_by_shared_junctions(event_rows):
            lines, etype_to_sig = [], {}
            for row in rows:
                for s, line in row_outliers(row, samples, event_support).items():
                    lines.append(line)
                    etype_to_sig.setdefault(row['event_type'], {}).setdefault(_row_junctions(row), {})[s] = line
            _write_lines(out, _by_deviation(lines))
            _write_lines(out_filtered, filter_outliers(etype_to_sig))


def main():
    args = build_parser().parse_args()
    with cli.ErrorHandler():
        spliceevents_outliers(args.input, args.output, args.psi, args.event_support)


if __name__ == '__main__':
    main()
