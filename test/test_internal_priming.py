from flair.isoform_data import check_intprim, internally_primed_end


def test_internal_priming_is_60_percent_a_in_the_20_bases_after_the_end():
    # SQANTI3's rule: at least 12 of the 20 bases A
    assert check_intprim('A' * 12 + 'C' * 8) == 12
    assert check_intprim('AC' * 10 + 'AAAA') == 0  # half A, then A past the 20 bases
    assert check_intprim('A' * 11 + 'G' * 9) == 0


def test_only_the_20_bases_after_the_end_count():
    assert check_intprim('C' * 20 + 'A' * 30) == 0
    assert check_intprim('A' * 15) == 0


def test_a_tail_clipped_past_a_rich_genome_is_not_internal_priming():
    after, aligned = 'A' * 14 + 'C' * 6, 'ACGT' * 5
    assert internally_primed_end(after, aligned, tail=0)
    assert not internally_primed_end(after, aligned, tail=20)


def test_a_read_ending_in_a_genomic_a_run_is_primed_whatever_is_clipped():
    # aligned through a genomic A run, A clipped after it: the rest of the primer
    after, aligned = 'GAAAGAGAGCATTACCCCAG', 'A' * 20
    assert internally_primed_end(after, aligned, tail=16)
    assert not internally_primed_end(after, 'A' * 15 + 'CGTCG', tail=0)
