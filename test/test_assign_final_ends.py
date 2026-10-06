from types import SimpleNamespace
from flair.isoform_data import Isoform, Junc
from flair.flair_transcriptome import assign_final_ends, assign_single_exon_isoforms

JUNCS = (Junc(1200, 1300), Junc(1400, 1500))


def _line(name, start, end):
    "count_sam_transcripts' ends line for a read aligned from start to end, in genomic coordinates"
    return f'{name}\ttx\t0\t{1200 - start}\tNone\t1\t{end - 1500}\tNone\n'


def _run(tmp_path, read_ends, max_ends, extra_lines=(), frac_support=0.05, normalize_ends=False):
    # a first-pass isoform with padded, normalized ends, and the reads assigned to it
    iso = Isoform('chr1', '+', JUNCS, 800, 2200)
    iso.name = 'tx'
    # the ends it shares with its gene's isoforms with similar terminal splice sites
    iso.unpadded_ends = (900, 2100)
    ends = tmp_path / 'ends.tsv'
    ends.write_text(''.join(_line(f'r{i}', s, e) for i, (s, e) in enumerate(read_ends)) + ''.join(extra_lines))
    read_map = tmp_path / 'map.txt'
    final = assign_final_ends({'tx': iso}, str(ends), str(read_map), normalize_ends=normalize_ends,
                              max_ends=max_ends, end_window=100, sjc_support=2, frac_support=frac_support)
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


def test_a_supported_chain_whose_end_variants_all_fail_is_one_isoform(tmp_path):
    # 3 and 4 reads of 7 in the gene each fail frac_support 0.6; the chain's 7 pass,
    # so it is one isoform at the densest ends of all its reads, as with max_ends 1
    read_ends = [(1000, 2000)] * 3 + [(1100, 1800)] * 3 + [(1105, 1790)]
    final, assigned, _ = _run(tmp_path, read_ends, max_ends=2, frac_support=0.6)
    ((name, iso),) = final.items()
    assert (iso.start, iso.end) == (1100, 1800)
    assert assigned == {name: {f'r{i}' for i in range(7)}}


def test_end_variants_passing_support_are_kept(tmp_path):
    read_ends = [(1000, 2000)] * 3 + [(1100, 1800)] * 3 + [(1105, 1790)]
    final, assigned, _ = _run(tmp_path, read_ends, max_ends=2, frac_support=0.4)
    assert sorted(len(assigned[name]) for name in final) == [3, 4]


def test_normalize_ends_keeps_the_shared_ends_and_all_the_reads(tmp_path):
    read_ends = [(1000, 2000)] * 3 + [(1100, 1800)] * 3
    final, assigned, _ = _run(tmp_path, read_ends, max_ends=2, normalize_ends=True)
    assert list(final) == ['tx'] and final['tx'].unpadded_ends == (900, 2100)
    assert assigned == {'tx': {f'r{i}' for i in range(6)}}


def test_single_exon_isoforms_are_left_for_assign_single_exon_isoforms(tmp_path):
    iso = Isoform('chr1', '-', (), 1000, 2000, reads=[SimpleNamespace(name='r', start=1000, end=2000)])
    iso.name = 'se'
    ends = tmp_path / 'ends.tsv'
    ends.write_text('')
    final = assign_final_ends({'se': iso}, str(ends), None, normalize_ends=False, max_ends=1, end_window=100,
                              sjc_support=2, frac_support=0.05)
    assert final == {'se': iso} and ends.read_text() == ''


def _se_isoform(name, start, end, read_ends, cluster=None, gene='g'):
    "a single-exon first-pass isoform, with the reads it was built from"
    reads = [SimpleNamespace(name=f'{name}{i}', start=s, end=e) for i, (s, e) in enumerate(read_ends)]
    iso = Isoform('chr1', '-', (), start, end, reads=reads, gene_id=gene)
    iso.name = name
    iso.end_variant_cluster = cluster
    return iso


def _assign_se(isoforms, frac_support=0.05, spliced_reads=0):
    iso_to_counts, gene_to_tot = {}, {'g': [spliced_reads, spliced_reads, spliced_reads]}
    final = assign_single_exon_isoforms(isoforms, iso_to_counts, gene_to_tot, se_support=3, frac_support=frac_support)
    return final, iso_to_counts, gene_to_tot


def test_single_exon_isoform_keeps_its_first_pass_ends_and_reads():
    iso = _se_isoform('se', 1010, 1990, [(1010, 1990)] * 3 + [(1100, 1700)])
    final, counts, gene_to_tot = _assign_se([iso], spliced_reads=10)
    assert final == {'se': [iso]} and counts == {'se': [4, 4]}
    # its reads are full-length reads of its gene, with the spliced ones
    assert gene_to_tot == {'g': [10, 14, 14]}


def test_single_exon_isoform_keeps_only_reads_covering_more_than_half_of_it():
    iso = _se_isoform('se', 1000, 2000, [(1000, 2000)] * 3 + [(1000, 1400), (1600, 3000)])
    final, counts, _ = _assign_se([iso])
    assert [read.name for read in final['se'][0].reads] == ['se0', 'se1', 'se2'] and counts['se'] == [3, 3]


def test_single_exon_cluster_whose_end_variants_all_fail_is_one_isoform():
    # 3 and 4 reads of 7 fail frac_support 0.6, the cluster's 7 pass
    a = _se_isoform('a', 1000, 2000, [(1000, 2000)] * 3, cluster='c')
    b = _se_isoform('b', 1500, 3000, [(1500, 3000)] * 4, cluster='c')
    final, counts, _ = _assign_se([a, b], frac_support=0.6)
    ((merged,),) = final.values()
    # the longest read's ends; a's reads cover only a third of it
    assert list(final) == ['a'] and (merged.start, merged.end) == (1500, 3000)
    assert counts == {merged.name: [4, 4]}


def test_single_exon_end_variants_passing_support_are_kept():
    a = _se_isoform('a', 1000, 2000, [(1000, 2000)] * 3, cluster='c')
    b = _se_isoform('b', 1500, 3000, [(1500, 3000)] * 4, cluster='c')
    final, _, _ = _assign_se([a, b], frac_support=0.4)
    assert final == {'a': [a], 'b': [b]}


def test_single_exon_support_is_a_fraction_of_all_the_genes_spliced_reads():
    # 7 of the gene's 7 single-exon and 10 spliced reads: the variants' 3 and 4
    # fail frac_support 0.3, the cluster's 7 pass
    a = _se_isoform('a', 1000, 2000, [(1000, 2000)] * 3, cluster='c')
    b = _se_isoform('b', 1500, 3000, [(1500, 3000)] * 4, cluster='c')
    final, _, _ = _assign_se([a, b], frac_support=0.3, spliced_reads=10)
    assert len(final) == 1
