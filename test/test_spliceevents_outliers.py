"""spliceevents_outliers: outlier samples of flair spliceevents' events"""
import pytest
from flair import FlairInputDataError
from flair.spliceevents_outliers import (row_outliers, filter_outliers, read_events, spliceevents_outliers,
                                         NO_DEVIATION, OUTLIER_COLUMNS)

SAMPLES = ['s1', 's2', 's3', 's4', 's5', 's6']


def _row(counts, totals, eventname='inc_es-of-chr1:100-200', eventtype='es', gene='G1'):
    return {'eventname': eventname, 'eventtype': eventtype, 'gene': gene,
            'counts': dict(zip(SAMPLES, counts)), 'totals': dict(zip(SAMPLES, totals))}


def test_a_sample_far_from_the_median_is_an_outlier():
    # PSIs 0.2, 0.2, 0.22, 0.18, 0.2 and 0.6
    outliers = row_outliers(_row([20, 20, 22, 18, 20, 60], [100] * 6), SAMPLES, 20)
    assert list(outliers) == ['s6']
    line = outliers['s6']
    assert dict(zip(OUTLIER_COLUMNS, line))['event_reads;total_locus_reads'] == '60;100'
    assert line[4] == 0.2 and line[6] == 0.6 and line[9] == 0.4


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


def _line(eventname, sample, dev_from_med):
    return [eventname, 'x', 'G1', sample, 0.2, 0.01, 0.6, 6, '60;100', 0.4, dev_from_med]


ES = (frozenset({'chr1:100-200', 'chr1:300-400'}), frozenset(), frozenset({'chr1:100-400'}))
# a complex exon skip of the same exon and outer junction
CES = (frozenset({'chr1:100-200', 'chr1:300-400'}), frozenset(), frozenset({'chr1:100-400'}))


def test_complex_exon_skip_explained_by_an_exon_skip_is_filtered():
    es, ces = _line('inc_es', 's6', 20), _line('inc_ces', 's6', 18)
    kept = filter_outliers({'es': {ES: {'s6': es}}, 'ces': {CES: {'s6': ces}}})
    assert kept == [es]


def test_complex_exon_skip_much_better_than_its_exon_skip_is_kept():
    # the exon skip's deviation is under 0.8 of the complex one's
    es, ces = _line('inc_es', 's6', 10), _line('inc_ces', 's6', 20)
    kept = filter_outliers({'es': {ES: {'s6': es}}, 'ces': {CES: {'s6': ces}}})
    assert kept == [ces, es]


def test_filtered_outliers_are_more_than_six_deviations_out():
    kept = filter_outliers({'ir': {(frozenset(), frozenset({'chr1:100-200'}), frozenset()):
                                   {'s5': _line('retained_ir', 's5', 6), 's6': _line('retained_ir', 's6', 6.5)}}})
    assert [line[3] for line in kept] == ['s6']


COUNTS_HEADER = ['eventname', 'event', 'eventtype', 'gene', 'junctions_included', 'junctions_excluded',
                 'outer_junctions', 'exons'] + SAMPLES
TOTALS_HEADER = ['event', 'eventtype', 'gene', 'total'] + SAMPLES


def _write_tables(tmp_path, counts_rows, totals_rows):
    prefix = str(tmp_path / 'events')
    with open(prefix + '.spliceevents.counts.tsv', 'w') as fh:
        for row in [COUNTS_HEADER] + counts_rows:
            fh.write('\t'.join(str(x) for x in row) + '\n')
    with open(prefix + '.spliceevents.totals.tsv', 'w') as fh:
        for row in [TOTALS_HEADER] + totals_rows:
            fh.write('\t'.join(str(x) for x in row) + '\n')
    return prefix


IR = 'ir-of-chr1:100-200(+)-G1'
IR_ROWS = [[f'retained_{IR}', IR, 'ir', 'G1', '', 'chr1:100-200', '', '', 20, 20, 22, 18, 20, 60],
           [f'spliced_{IR}', IR, 'ir', 'G1', 'chr1:100-200', '', '', '', 80, 80, 78, 82, 80, 40]]


def test_other_rows_are_left_out(tmp_path):
    prefix = _write_tables(tmp_path, IR_ROWS + [[f'other_{IR}', IR, 'ir', 'G1', '', '', '', '', 1, 1, 1, 1, 1, 1]],
                           [[IR, 'ir', 'G1', 'event_reads'] + [100] * 6, [IR, 'ir', 'G1', 'locus_reads'] + [200] * 6])
    gene_to_rows, samples = read_events(prefix, 'junc')
    assert samples == SAMPLES
    assert [row['eventname'] for row in gene_to_rows['G1']] == [f'retained_{IR}', f'spliced_{IR}']


def test_an_event_without_totals_is_an_error(tmp_path):
    prefix = _write_tables(tmp_path, IR_ROWS, [])
    with pytest.raises(FlairInputDataError):
        read_events(prefix, 'junc')


def test_psi_tot_divides_by_the_locus_reads(tmp_path):
    # of the event's 100 reads s6 has 60 retained, against the others' 20: both rows
    # outliers by PSIjunc.  Of the locus's reads, 400 for s6 and 200 for the others,
    # its retained reads are 0.15 against 0.1, under 0.1 from the median, while its
    # spliced reads are 0.1 against 0.4: only the spliced row by PSItot
    prefix = _write_tables(tmp_path, IR_ROWS, [[IR, 'ir', 'G1', 'event_reads'] + [100] * 6,
                                               [IR, 'ir', 'G1', 'locus_reads'] + [200] * 5 + [400]])
    for psi, outliers in (('junc', {(f'retained_{IR}', 's6'), (f'spliced_{IR}', 's6')}),
                          ('tot', {(f'spliced_{IR}', 's6')})):
        out = str(tmp_path / psi)
        spliceevents_outliers(prefix, out, psi, 20)
        lines = open(out + '.spliceevents.outliers.tsv').read().splitlines()
        assert lines[0].split('\t') == OUTLIER_COLUMNS
        assert {(line.split('\t')[0], line.split('\t')[3]) for line in lines[1:]} == outliers
