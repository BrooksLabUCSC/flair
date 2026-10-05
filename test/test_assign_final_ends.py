from flair.isoform_data import Isoform, Junc
from flair.flair_transcriptome import assign_final_ends

JUNCS = (Junc(1200, 1300), Junc(1400, 1500))


def _line(name, start, end):
    "count_sam_transcripts' ends line for a read aligned from start to end, in genomic coordinates"
    return f'{name}\ttx\t0\t{1200 - start}\tNone\t1\t{end - 1500}\tNone\n'


def _run(tmp_path, read_ends, max_ends, extra_lines=(), frac_support=0.05):
    # a first-pass isoform with padded, normalized ends, and the reads assigned to it
    iso = Isoform('chr1', '+', JUNCS, 800, 2200)
    iso.name = 'tx'
    ends = tmp_path / 'ends.tsv'
    ends.write_text(''.join(_line(f'r{i}', s, e) for i, (s, e) in enumerate(read_ends)) + ''.join(extra_lines))
    read_map = tmp_path / 'map.txt'
    final = assign_final_ends({'tx': iso}, str(ends), str(read_map), max_ends=max_ends, end_window=100,
                              sjc_support=2, se_support=3, frac_support=frac_support)
    assigned = {}
    for line in ends.read_text().splitlines():
        assigned.setdefault(line.split('\t')[1], set()).add(line.split('\t')[0])
    return final, assigned, read_map.read_text()


def test_one_end_variant_takes_the_densest_read_ends(tmp_path):
    final, assigned, read_map = _run(tmp_path, [(1000, 2000), (1002, 2003), (1004, 1990), (900, 1700)], max_ends=1)
    (iso,) = final.values()
    assert iso.unpadded_ends == (1000, 2003)
    assert assigned == {iso.name: {'r0', 'r1', 'r2', 'r3'}}
    assert read_map.split('\t')[0] == iso.name


def test_read_ends_come_from_the_realignment(tmp_path):
    # a read running into the padding gets its end there
    final, _, _ = _run(tmp_path, [(850, 2150)], max_ends=1)
    assert next(iter(final.values())).unpadded_ends == (850, 2150)


def test_more_max_ends_splits_the_reads_between_end_groups(tmp_path):
    read_ends = [(1000, 2000)] * 3 + [(1100, 1800)] * 3 + [(1105, 1790)]
    final, assigned, _ = _run(tmp_path, read_ends, max_ends=2)
    by_ends = {iso.unpadded_ends: assigned[name] for name, iso in final.items()}
    assert by_ends == {(1000, 2000): {'r0', 'r1', 'r2'}, (1100, 1800): {'r3', 'r4', 'r5', 'r6'}}


def test_a_line_without_distances_still_counts(tmp_path):
    final, assigned, _ = _run(tmp_path, [(1000, 2000), (1002, 2003)], max_ends=1,
                              extra_lines=['odd\ttx\tNone\tNone\t10\tNone\tNone\t10\n'])
    (name,) = final
    assert assigned[name] == {'r0', 'r1', 'odd'}


def test_a_supported_chain_whose_end_variants_all_fail_reports_its_best(tmp_path):
    # 3 and 4 reads of 7 in the gene each fail frac_support 0.6; the chain's 7 pass
    read_ends = [(1000, 2000)] * 3 + [(1100, 1800)] * 3 + [(1105, 1790)]
    final, assigned, _ = _run(tmp_path, read_ends, max_ends=2, frac_support=0.6)
    reported = {(iso.start, iso.end): (iso.report_unsupported, assigned[name]) for name, iso in final.items()}
    # the reads stay with their variants, and only the best supported is reported
    assert reported == {(1000, 2000): (False, {'r0', 'r1', 'r2'}),
                        (1100, 1800): (True, {'r3', 'r4', 'r5', 'r6'})}


def test_tied_end_variants_report_the_longer(tmp_path):
    read_ends = [(1000, 2000)] * 3 + [(1100, 1800)] * 3
    final, _, _ = _run(tmp_path, read_ends, max_ends=2, frac_support=0.6)
    assert [(iso.start, iso.end) for iso in final.values() if iso.report_unsupported] == [(1000, 2000)]


def test_end_variants_passing_support_are_kept(tmp_path):
    read_ends = [(1000, 2000)] * 3 + [(1100, 1800)] * 3 + [(1105, 1790)]
    final, _, _ = _run(tmp_path, read_ends, max_ends=2, frac_support=0.4)
    assert len(final) == 2 and not any(iso.report_unsupported for iso in final.values())
