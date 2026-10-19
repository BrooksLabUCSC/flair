"""spliceevents_outliers: outlier samples of flair spliceevents' events"""
import pytest
from flair import FlairInputDataError
from flair.spliceevents_outliers import (row_outliers, filter_outliers, read_events, spliceevents_outliers,
                                         group_by_shared_junctions, NO_DEVIATION, OUTLIER_COLUMNS)

SAMPLES = ['s1', 's2', 's3', 's4', 's5', 's6']


def _row(counts, totals, event='es|chr1:101-200(+)', alternative='inc', etype='es', gene='G1'):
    return {'event': event, 'alternative': alternative, 'event_type': etype, 'gene': gene, 'gene_name': 'NAME1',
            'additional_genes': '', 'additional_gene_names': '',
            'counts': dict(zip(SAMPLES, counts)), 'totals': dict(zip(SAMPLES, totals))}


def test_a_sample_far_from_the_median_is_an_outlier():
    # PSIs 0.2, 0.2, 0.22, 0.18, 0.2 and 0.6
    outliers = row_outliers(_row([20, 20, 22, 18, 20, 60], [100] * 6), SAMPLES, 20)
    assert list(outliers) == ['s6']
    line = outliers['s6']
    assert dict(zip(OUTLIER_COLUMNS, line))['event_reads;total_locus_reads'] == '60;100'
    columns = dict(zip(OUTLIER_COLUMNS, line))
    assert (columns['medianPSI'], columns['sample_val'], columns['delta_PSI_to_med']) == (0.2, 0.6, 0.4)


def test_a_small_difference_is_not_an_outlier_however_many_deviations():
    # 0.09 from the median, under 0.1, though the others agree exactly
    assert row_outliers(_row([20, 20, 20, 20, 20, 29], [100] * 6), SAMPLES, 20) == {}


def test_samples_agreeing_exactly_give_no_deviation():
    outliers = row_outliers(_row([20, 20, 20, 20, 20, 50], [100] * 6), SAMPLES, 20)
    assert outliers['s6'][-1] == NO_DEVIATION


def test_a_sample_with_exactly_event_support_reads_is_used():
    # s6 has 20 reads, the --event_support: its PSI counts, as in the PSI tables
    outliers = row_outliers(_row([20, 20, 22, 18, 20, 12], [100] * 5 + [20]), SAMPLES, 20)
    assert list(outliers) == ['s6']


def test_fewer_than_five_samples_with_a_psi_give_no_outliers():
    # two samples under --event_support leave four
    assert row_outliers(_row([20, 20, 22, 18, 5, 12], [100, 100, 100, 100, 10, 19]), SAMPLES, 20) == {}


def _line(event, sample, dev_from_med):
    return [event, 'inc', 'x', 'G1', 'NAME1', '', '', sample, 0.2, 0.01, 0.6, 6, '60;100', 0.4, dev_from_med]


ES = (frozenset({'chr1:100-200', 'chr1:300-400'}), frozenset(), frozenset({'chr1:100-400'}))
# a complex exon skip of the same exon and outer junction
CES = (frozenset({'chr1:100-200', 'chr1:300-400'}), frozenset(), frozenset({'chr1:100-400'}))


def test_complex_exon_skip_explained_by_an_exon_skip_is_filtered():
    es, ces = _line('es', 's6', 20), _line('ces', 's6', 18)
    kept = filter_outliers({'es': {ES: {'s6': es}}, 'ces': {CES: {'s6': ces}}})
    assert kept == [es]


def test_complex_exon_skip_much_better_than_its_exon_skip_is_kept():
    # the exon skip's deviation is under 0.8 of the complex one's
    es, ces = _line('es', 's6', 10), _line('ces', 's6', 20)
    kept = filter_outliers({'es': {ES: {'s6': es}}, 'ces': {CES: {'s6': ces}}})
    assert kept == [ces, es]


# an alt3 of the exon skip's first junction, its 5' site spliced to the exon or past it:
# the exon skip's flanking and skipping junctions
ALT3 = (frozenset({'chr1:100-200'}), frozenset({'chr1:100-400'}), frozenset())


def test_alternative_site_explained_by_an_exon_skip_is_filtered():
    es, alt3 = _line('es', 's6', 20), _line('alt3', 's6', 18)
    kept = filter_outliers({'es': {ES: {'s6': es}}, 'alt3': {ALT3: {'s6': alt3}}})
    assert kept == [es]


def test_alternative_site_much_better_than_its_exon_skip_is_kept():
    es, alt3 = _line('es', 's6', 10), _line('alt3', 's6', 20)
    kept = filter_outliers({'es': {ES: {'s6': es}}, 'alt3': {ALT3: {'s6': alt3}}})
    assert kept == [alt3, es]


def test_alternative_site_explained_by_one_of_its_events_is_filtered_whatever_the_others():
    # a complex exon skip of two exons, between the same outer sites, also holds the
    # alt3's junctions but has no outlier of s6: the exon skip still explains it
    ces = (frozenset({'chr1:100-200', 'chr1:250-300', 'chr1:350-400'}), frozenset(), frozenset({'chr1:100-400'}))
    es, alt3 = _line('es', 's6', 20), _line('alt3', 's6', 18)
    ces_line = _line('ces', 's5', 20)
    kept = filter_outliers({'es': {ES: {'s6': es}}, 'ces': {ces: {'s5': ces_line}}, 'alt3': {ALT3: {'s6': alt3}}})
    assert kept == [es, ces_line]


def test_alternative_site_sharing_no_event_junctions_is_kept():
    # an alt5 of other junctions than the exon skip's
    alt5 = (frozenset({'chr1:500-600'}), frozenset({'chr1:550-600'}), frozenset())
    es, alt = _line('es', 's6', 20), _line('alt5', 's6', 18)
    kept = filter_outliers({'es': {ES: {'s6': es}}, 'alt5': {alt5: {'s6': alt}}})
    assert kept == [es, alt]


def test_starts_and_ends_are_left_out_of_the_filtered_outliers():
    # on purpose, for their many false positives: see LEFT_OUT_EVENTS
    tss = (frozenset({'101'}), frozenset({'201'}), frozenset())
    es = _line('es', 's6', 20)
    assert filter_outliers({'es': {ES: {'s6': es}}, 'tss': {tss: {'s6': _line('tss', 's6', 30)}}}) == [es]


def test_an_unknown_event_type_is_an_error():
    with pytest.raises(ValueError):
        filter_outliers({'new': {ES: {'s6': _line('new', 's6', 20)}}})


def test_complex_exon_skip_is_compared_with_the_exon_skips_inclusion_row():
    # by PSItot, an exon skip's rows needn't mirror each other: the exclusion row
    # has an outlier of s6, the inclusion row, the complex skip's match, has none
    es_exc_juncs = (frozenset(), ES[0], ES[2])
    es_exc, ces = _line('es', 's6', 20), _line('ces', 's6', 18)
    kept = filter_outliers({'es': {ES: {}, es_exc_juncs: {'s6': es_exc}}, 'ces': {CES: {'s6': ces}}})
    assert kept == [es_exc, ces]


def test_filtered_outliers_are_more_than_six_deviations_out():
    kept = filter_outliers({'ir': {(frozenset(), frozenset({'chr1:100-200'}), frozenset()):
                                   {'s5': _line('ir', 's5', 6), 's6': _line('ir', 's6', 6.5)}}})
    assert [line[OUTLIER_COLUMNS.index('sample')] for line in kept] == ['s6']


COUNTS_HEADER = ['event', 'alternative', 'event_type', 'gene', 'gene_name', 'additional_genes', 'additional_gene_names',
                 'junctions_included', 'junctions_excluded',
                 'outer_junctions', 'exons', 'novel_junctions'] + SAMPLES
TOTALS_HEADER = ['event', 'event_type', 'gene', 'gene_name', 'additional_genes', 'additional_gene_names', 'total'] + SAMPLES


def _write_tables(tmp_path, counts_rows, totals_rows):
    prefix = str(tmp_path / 'events')
    with open(prefix + '.spliceevents.counts.tsv', 'w') as fh:
        for row in [COUNTS_HEADER] + counts_rows:
            fh.write('\t'.join(str(x) for x in row) + '\n')
    with open(prefix + '.spliceevents.totals.tsv', 'w') as fh:
        for row in [TOTALS_HEADER] + totals_rows:
            fh.write('\t'.join(str(x) for x in row) + '\n')
    return prefix


IR = 'ir|chr1:101-200(+)'
IR_ROWS = [[IR, 'retained', 'ir', 'G1', 'NAME1', '', '', '', 'chr1:101-200', '', '', '', 20, 20, 22, 18, 20, 60],
           [IR, 'spliced', 'ir', 'G1', 'NAME1', '', '', 'chr1:101-200', '', '', '', '', 80, 80, 78, 82, 80, 40]]


def test_other_rows_are_left_out(tmp_path):
    prefix = _write_tables(tmp_path, IR_ROWS + [[IR, 'other', 'ir', 'G1', 'NAME1', '', '', '', '', '', '', '', 1, 1, 1, 1, 1, 1]],
                           [[IR, 'ir', 'G1', 'NAME1', '', '', 'event_reads'] + [100] * 6, [IR, 'ir', 'G1', 'NAME1', '', '', 'locus_reads'] + [200] * 6])
    rows, samples = read_events(prefix, 'junc')
    assert samples == SAMPLES
    assert [row['alternative'] for row in rows] == ['retained', 'spliced']


def test_an_event_without_totals_is_an_error(tmp_path):
    prefix = _write_tables(tmp_path, IR_ROWS, [])
    with pytest.raises(FlairInputDataError):
        read_events(prefix, 'junc')


def test_psi_tot_divides_by_the_locus_reads(tmp_path):
    # of the event's 100 reads s6 has 60 retained, against the others' 20: both rows
    # outliers by PSIjunc.  Of the locus's reads, 400 for s6 and 200 for the others,
    # its retained reads are 0.15 against 0.1, under 0.1 from the median, while its
    # spliced reads are 0.1 against 0.4: only the spliced row by PSItot
    prefix = _write_tables(tmp_path, IR_ROWS, [[IR, 'ir', 'G1', 'NAME1', '', '', 'event_reads'] + [100] * 6,
                                               [IR, 'ir', 'G1', 'NAME1', '', '', 'locus_reads'] + [200] * 5 + [400]])
    for psi, outliers in (('junc', {('retained', 's6'), ('spliced', 's6')}),
                          ('tot', {('spliced', 's6')})):
        out = str(tmp_path / psi)
        spliceevents_outliers(prefix, out, psi, 20)
        lines = open(out + '.spliceevents.outliers.tsv').read().splitlines()
        assert lines[0].split('\t') == OUTLIER_COLUMNS
        sample = OUTLIER_COLUMNS.index('sample')
        assert {(line.split('\t')[1], line.split('\t')[sample]) for line in lines[1:]} == outliers


def _event_row(event, included='', excluded='', outer=''):
    return {'event': event, 'junctions_included': included, 'junctions_excluded': excluded, 'outer_junctions': outer}


def test_events_sharing_junctions_are_compared_whatever_their_genes():
    es_inc = _event_row('es|chr1:201-300(+)', 'chr1:101-200,chr1:301-400', '', 'chr1:101-400')
    es_exc = _event_row('es|chr1:201-300(+)', '', 'chr1:101-200,chr1:301-400', 'chr1:101-400')
    # an alt3 sharing the exon skip's first intron, and through it an alt5 sharing none of the skip's
    alt3 = _event_row('alt3|chr1:101-250(+)', 'chr1:101-200', 'chr1:101-250')
    alt5 = _event_row('alt5|chr1:90-250(+)', 'chr1:90-250', 'chr1:101-250')
    apart = _event_row('ir|chr1:1001-1100(+)', 'chr1:1001-1100')
    tss1, tss2 = _event_row('tss|chr1:50-90(+)'), _event_row('tss|chr1:5000-5100(+)')
    groups = group_by_shared_junctions([es_inc, alt3, es_exc, apart, alt5, tss1, tss2])
    assert [[row['event'] for row in group] for group in groups] == [
        ['es|chr1:201-300(+)', 'alt3|chr1:101-250(+)', 'es|chr1:201-300(+)', 'alt5|chr1:90-250(+)'],
        ['ir|chr1:1001-1100(+)'], ['tss|chr1:50-90(+)'], ['tss|chr1:5000-5100(+)']]
