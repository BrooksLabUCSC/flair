#!/usr/bin/env python3
"""
Outlier samples of flair spliceevents' events: samples whose PSI of an event row is
far from the median of all samples', from OUTPUT.spliceevents.counts.tsv and
OUTPUT.spliceevents.totals.tsv.
"""
import argparse
import csv
from collections import defaultdict
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

OUTLIER_COLUMNS = ['event', 'alternative', 'event_type', 'gene', 'sample', 'medianPSI', 'dev(IQR/2)', 'sample_val',
                   'tot_not_NA_samples', 'event_reads;total_locus_reads', 'delta_PSI_to_med', 'dev_from_med']
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
                             'deviations from the median and not explained by a simpler event')
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
    """the event rows of the counts table, by gene in table order, each with its
    per-sample counts and the per-sample totals of its PSI, and the samples"""
    header, rows = read_tsv(input_prefix + '.spliceevents.counts.tsv')
    samples = header[header.index('exons') + 1:]
    _, total_rows = read_tsv(input_prefix + '.spliceevents.totals.tsv')
    totals = {row['event']: {s: int(row[s]) for s in samples} for row in total_rows if row['total'] == PSI_TOTALS[psi]}
    gene_to_rows = defaultdict(list)
    for row in rows:
        if row['alternative'] == 'other':  # the reads of an event's other ends or junctions
            continue
        if row['event'] not in totals:
            raise FlairInputDataError(f"event {row['event']} of {input_prefix}.spliceevents.counts.tsv is not in "
                                      f'{input_prefix}.spliceevents.totals.tsv')
        row['counts'] = {s: int(row[s]) for s in samples}
        row['totals'] = totals[row['event']]
        gene_to_rows[row['gene']].append(row)
    return gene_to_rows, samples


def _junctions(field):
    return frozenset(field.split(',')) if field else frozenset()


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
            outliers[s] = [row['event'], row['alternative'], row['event_type'], row['gene'], s, round(med, 6), round(dev, 6), round(psi, 6),
                           len(vals), f"{row['counts'][s]};{row['totals'][s]}", round(delta, 6),
                           round(delta / dev, 6) if dev != 0 else NO_DEVIATION]
    return outliers


def _by_deviation(lines):
    return sorted(lines, reverse=True, key=lambda x: x[::-1])


def _better_than(this, other):
    "is this outlier line's deviation enough better than the other's"
    return other[-1] / this[-1] < BETTER_DEVIATION_FRACTION


def filter_outliers(etype_to_sig):  # noqa: C901 - FIXME: reduce complexity
    """A gene's outliers, by event type and then junctions, (included, excluded,
    outer), and sample, kept when a simpler event of the same junctions doesn't
    explain them as well: a complex exon skip by an exon skip, an alternative splice
    site by an exon skip or alternative first or last exon"""
    seen_junctions, seen_es_juncs, kept = {}, {}, []
    # alternative splice site events overlapping alternative first and last exons
    # are already left out when the events are found
    for et in ['es', 'afe', 'ale', 'ir', 'alt5ss', 'alt3ss']:
        for juncs, sample_to_vals in etype_to_sig.get(et, {}).items():
            if et != 'ir':
                inc_juncs, exc_juncs, outer_juncs = juncs
                for j1 in inc_juncs:
                    for j2 in outer_juncs:
                        seen_junctions.setdefault(frozenset((j1, j2)), {})[juncs] = sample_to_vals
                if et == 'es':
                    seen_es_juncs[(inc_juncs | exc_juncs) | outer_juncs] = sample_to_vals
            kept.extend(sample_to_vals.values())

    for et in ['ces', 'alt5', 'alt3']:
        for juncs, sample_to_vals in etype_to_sig.get(et, {}).items():
            inc_juncs, exc_juncs, outer_juncs = juncs
            good_s = set()
            has_overlap = False
            if et == 'ces':  # against the exon skips whose junctions hold its own
                # the exclusion row, by its exons' junctions, else the exon included
                all_juncs = (exc_juncs if len(inc_juncs) == 0 else inc_juncs) | outer_juncs
                for es_juncs, es_vals in seen_es_juncs.items():
                    if all_juncs <= es_juncs:
                        has_overlap = True
                        good_s.update(s for s in sample_to_vals
                                      if s not in es_vals or _better_than(sample_to_vals[s], es_vals[s]))
            else:
                all_juncs = inc_juncs | exc_juncs
                for jkey, others in seen_junctions.items():
                    if all_juncs <= jkey:
                        has_overlap = True
                        for other_vals in others.values():
                            good_s.update(s for s in sample_to_vals
                                          if s not in other_vals or _better_than(sample_to_vals[s], other_vals[s]))
            if not has_overlap:
                good_s = sample_to_vals.keys()
                seen_junctions.setdefault(all_juncs, {})[juncs] = sample_to_vals
            kept.extend(sample_to_vals[s] for s in good_s)
    return [line for line in _by_deviation(kept) if line[-1] > FILTERED_MIN_DEVIATIONS]


def spliceevents_outliers(input_prefix, output, psi, event_support):
    gene_to_rows, samples = read_events(input_prefix, psi)
    with open(output + '.spliceevents.outliers.tsv', 'w') as out, \
            open(output + '.spliceevents.outliers.filtered.tsv', 'w') as out_filtered:
        for fh in (out, out_filtered):
            fh.write('\t'.join(OUTLIER_COLUMNS) + '\n')
        # each gene's events alone, as an outlier is only explained by an event of its gene
        for rows in gene_to_rows.values():
            sig_events, etype_to_sig = [], {}
            for row in rows:
                outliers = row_outliers(row, samples, event_support)
                juncs = (_junctions(row['junctions_included']), _junctions(row['junctions_excluded']),
                         _junctions(row['outer_junctions']))
                for s, line in outliers.items():
                    sig_events.append(line)
                    etype_to_sig.setdefault(row['event_type'], {}).setdefault(juncs, {})[s] = line
            for line in _by_deviation(sig_events):
                out.write('\t'.join(str(x) for x in line) + '\n')
            for line in filter_outliers(etype_to_sig):
                out_filtered.write('\t'.join(str(x) for x in line) + '\n')


def main():
    args = build_parser().parse_args()
    with cli.ErrorHandler():
        spliceevents_outliers(args.input, args.output, args.psi, args.event_support)


if __name__ == '__main__':
    main()
