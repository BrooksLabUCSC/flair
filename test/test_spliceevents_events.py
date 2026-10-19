"""
flair spliceevents' events, from each sample's corrected reads: which events are
found, and their rows' counts, PSIs and totals.

Reads are given as (junctions, start, end, number of reads), on one gene G1 on chr1,
with 0-based starts, and the tables' rows are keyed by name, event|alternative, 1-based.
"""
import csv
from flair import SeqRange
from flair.spliceevents_events import is_alternative, extract_terminal_info, group_ends, Thresholds
from flair.spliceevents_genes import process_region_events, read_event_genes, AnnotatedSplicing
from flair.spliceevents_reads import ReadEnds, ReadChain
from flair.spliceevents_tables import EVENT_ROW_COLUMNS

SAMPLES = ['s1', 's2', 's3']
CHROM = 'chr1'


def _sample_reads(chains, gene='G1', strand='+'):
    "one sample's ReadChains, all of gene"
    by_juncs = {}
    for juncs, start, end, count in chains:
        by_juncs.setdefault(tuple(juncs), []).extend([ReadEnds(start, end)] * count)
    return [ReadChain((gene,), juncs, strand, ends) for juncs, ends in by_juncs.items()]


def _events(tmp_path, *samples_chains, strand='+', feature_support=2, event_frac_of_tot=0.1,
            junc_frac_of_event=0.05, event_support=20, annot_afe_ss=None, annot_ale_ss=None, annotated=None):
    """process_region_events on the reads of each sample: the rows of each table, as
    {row name: {column: value}}, keyed by table (counts, PSIjunc, PSItot), and the
    totals, as {(event, total): {sample: count}}"""
    samples = SAMPLES[:len(samples_chains)]
    prefix = str(tmp_path / 'region')
    sample_chains = [_sample_reads(chains, strand=strand) for chains in samples_chains]
    # each gene's junctions, those of its reads, as a novel gene's
    gene_juncs = {}
    for chains in sample_chains:
        for chain in chains:
            for gene in chain.genes:
                gene_juncs.setdefault(gene, set()).update(chain.juncs)
    process_region_events(prefix, CHROM, samples, sample_chains,
                          {}, Thresholds(feature_support, event_frac_of_tot, junc_frac_of_event, event_support),
                          annot_afe_ss or {}, annot_ale_ss or {}, [SeqRange(CHROM, 0, 10 ** 9)], gene_juncs, annotated)
    columns = EVENT_ROW_COLUMNS + samples
    tables = {}
    for table in ('counts', 'PSIjunc', 'PSItot'):
        with open(f'{prefix}.spliceevents.{table}.tsv') as fh:
            tables[table] = {f"{row['event']}|{row['alternative']}": row
                             for row in csv.DictReader(fh, fieldnames=columns, delimiter='\t')}
    totals = {}
    with open(f'{prefix}.spliceevents.totals.tsv') as fh:
        for line in fh:
            event, _, _, _, _, _, total, *counts = line.rstrip('\n').split('\t')
            totals[(event, total)] = dict(zip(samples, (int(c) for c in counts)))
    return tables, totals


def _of_type(tables, eventtype):
    return {name: row for name, row in tables['counts'].items() if row['event_type'] == eventtype}


# exons 100-200, 300-400 and 500-600 on +: the middle one included or skipped
INC = [(200, 300), (400, 500)]
EXC = [(200, 500)]
ES = 'es|chr1:301-400(+)'


def test_a_skipped_exon_is_an_es_event_with_its_inclusion_and_exclusion(tmp_path):
    tables, totals = _events(tmp_path, [(INC, 100, 600, 30), (EXC, 100, 600, 10)],
                             [(INC, 100, 600, 10), (EXC, 100, 600, 30)])
    es = _of_type(tables, 'es')
    assert sorted(es) == [f'{ES}|exc', f'{ES}|inc']
    inc, exc = tables['counts'][f'{ES}|inc'], tables['counts'][f'{ES}|exc']
    assert (inc['s1'], inc['s2'], exc['s1'], exc['s2']) == ('30', '10', '10', '30')
    # 1-based, the introns as STAR writes them
    assert inc['exons'] == 'chr1:301-400' and exc['outer_junctions'] == 'chr1:201-500'
    assert inc['junctions_included'] == 'chr1:201-300,chr1:401-500'
    assert (tables['PSIjunc'][f'{ES}|inc']['s1'], tables['PSIjunc'][f'{ES}|inc']['s2']) == ('0.75', '0.25')
    assert totals[(ES, 'event_reads')] == {'s1': 40, 's2': 40}


def test_psitot_divides_by_all_reads_covering_the_event(tmp_path):
    # 40 more reads in s1 across the exon's locus, spliced only downstream of it: reads
    # of the locus, not of the event
    tables, totals = _events(tmp_path, [(INC, 100, 600, 30), (EXC, 100, 600, 10), ([(620, 700)], 150, 800, 40)],
                             [(INC, 100, 600, 10), (EXC, 100, 600, 30)])
    assert totals[(ES, 'locus_reads')]['s1'] == 80
    assert tables['PSItot'][f'{ES}|exc']['s1'] == '0.125'
    assert tables['PSIjunc'][f'{ES}|exc']['s1'] == '0.25'


def test_a_sample_with_fewer_event_reads_than_event_support_has_no_psi(tmp_path):
    tables, _ = _events(tmp_path, [(INC, 100, 600, 30), (EXC, 100, 600, 10)],
                        [(INC, 100, 600, 10), (EXC, 100, 600, 30)],
                        [(INC, 100, 600, 5), (EXC, 100, 600, 14)])
    assert tables['counts'][f'{ES}|inc']['s3'] == '5'
    assert tables['PSIjunc'][f'{ES}|inc']['s3'] == ''


def test_an_event_with_a_psi_in_only_one_sample_is_reported(tmp_path):
    # s2 under event_support leaves s1's PSI, 0.75, alone: an alternative exon, as in one sample
    tables, _ = _events(tmp_path, [(INC, 100, 600, 30), (EXC, 100, 600, 10)],
                        [(INC, 100, 600, 5), (EXC, 100, 600, 14)])
    assert tables['PSIjunc'][f'{ES}|inc']['s1'] == '0.75'
    assert tables['PSIjunc'][f'{ES}|inc']['s2'] == ''


def test_an_event_with_a_psi_in_only_one_sample_is_not_reported_when_it_is_constitutive(tmp_path):
    tables, _ = _events(tmp_path, [(INC, 100, 600, 300), (EXC, 100, 600, 2)],
                        [(INC, 100, 600, 5), (EXC, 100, 600, 2)], junc_frac_of_event=0.05)
    assert _of_type(tables, 'es') == {}


def test_an_alternative_has_a_psi_away_from_0_and_1_in_some_sample():
    # alike in every sample, but alternative
    assert is_alternative([0.5, 0.51], 0.05)
    assert is_alternative([0.5], 0.05)
    assert is_alternative([0.99, 0.5], 0.05)
    # used alone in one sample and not at all in another
    assert is_alternative([0.0, 1.0], 0.05)
    assert is_alternative([0.97, 0.03], 0.05)
    # differing between samples, but never alternative
    assert not is_alternative([0.01, 0.04], 0.05)
    assert not is_alternative([0.99], 0.05)
    assert not is_alternative([0.97, 0.99], 0.05)
    assert not is_alternative([], 0.05)


def test_a_skip_in_fewer_reads_than_feature_support_is_not_an_event(tmp_path):
    tables, _ = _events(tmp_path, [(INC, 100, 600, 30), (EXC, 100, 600, 1)],
                        [(INC, 100, 600, 30)])
    assert _of_type(tables, 'es') == {}


def test_an_exon_skipped_alike_in_every_sample_is_an_event(tmp_path):
    # whether its PSI differs between samples is for what reads the tables
    tables, _ = _events(tmp_path, [(INC, 100, 600, 30), (EXC, 100, 600, 10)],
                        [(INC, 100, 600, 30), (EXC, 100, 600, 10)])
    assert sorted(_of_type(tables, 'es')) == [f'{ES}|exc', f'{ES}|inc']


def test_an_exon_included_in_one_sample_and_skipped_in_the_other_is_an_event(tmp_path):
    tables, _ = _events(tmp_path, [(INC, 100, 600, 30)], [(EXC, 100, 600, 30)])
    assert sorted(_of_type(tables, 'es')) == [f'{ES}|exc', f'{ES}|inc']


def test_an_exon_skipped_in_too_few_reads_in_every_sample_is_not_an_event(tmp_path):
    # skipped in 3% and 4% of the reads, under junc_frac_of_event, though they differ
    tables, _ = _events(tmp_path, [(INC, 100, 600, 300), (EXC, 100, 600, 9)],
                        [(INC, 100, 600, 300), (EXC, 100, 600, 12)])
    assert _of_type(tables, 'es') == {}


def test_one_sample_needs_its_psi_between_junc_frac_of_event_and_one_minus_it(tmp_path):
    assert _of_type(_events(tmp_path, [(INC, 100, 600, 30), (EXC, 100, 600, 10)])[0], 'es')
    assert not _of_type(_events(tmp_path, [(INC, 100, 600, 300), (EXC, 100, 600, 2)], junc_frac_of_event=0.05)[0], 'es')


# exons 100-200, 300-400, 500-600 and 700-800 on +: the two middle ones included or skipped together
INC2 = [(200, 300), (400, 500), (600, 700)]
EXC2 = [(200, 700)]
CES = 'ces|chr1:201-700(+)'


def _ces_counts(tables):
    rows = _of_type(tables, 'ces')
    assert all(row['event'] == CES for row in rows.values())
    return {row['alternative']: {s: row[s] for s in SAMPLES[:2]} for row in rows.values()}


def test_two_exons_skipped_together_are_a_ces_event(tmp_path):
    tables, _ = _events(tmp_path, [(INC2, 100, 800, 30), (EXC2, 100, 800, 10)],
                        [(INC2, 100, 800, 10), (EXC2, 100, 800, 30)])
    assert _ces_counts(tables) == {'exc': {'s1': '10', 's2': '30'},
                                   '301-400,501-600': {'s1': '30', 's2': '10'}}


def test_a_ces_events_counts_do_not_depend_on_the_order_of_its_reads(tmp_path):
    # the skipping reads first, as when the bam has them first
    tables, _ = _events(tmp_path, [(EXC2, 100, 800, 10), (INC2, 100, 800, 30)],
                        [(EXC2, 100, 800, 30), (INC2, 100, 800, 10)])
    assert _ces_counts(tables) == {'exc': {'s1': '10', 's2': '30'},
                                   '301-400,501-600': {'s1': '30', 's2': '10'}}


def test_a_ces_events_inclusion_counts_the_samples_without_skipping(tmp_path):
    # s1 has no skipping reads, s2's come after s1 is counted
    tables, _ = _events(tmp_path, [(INC2, 100, 800, 30)],
                        [(INC2, 100, 800, 10), (EXC2, 100, 800, 30)])
    assert _ces_counts(tables)['301-400,501-600'] == {'s1': '30', 's2': '10'}


# from a donor at 200, acceptors at 300 or 350, then on to an exon at 500-600
NEAR = [(200, 300), (400, 500)]
FAR = [(200, 350), (400, 500)]


def test_two_acceptors_of_one_donor_are_an_alt3_event_on_plus(tmp_path):
    tables, _ = _events(tmp_path, [(NEAR, 100, 600, 30), (FAR, 100, 600, 10)],
                        [(NEAR, 100, 600, 10), (FAR, 100, 600, 30)])
    alt3 = _of_type(tables, 'alt3')
    # from the shared donor to the furthest acceptor, rows by their introns
    assert sorted(alt3) == ['alt3|chr1:201-350(+)|201-300', 'alt3|chr1:201-350(+)|201-350']
    assert alt3['alt3|chr1:201-350(+)|201-350']['s2'] == '30'
    assert _of_type(tables, 'alt5') == {}


def test_the_same_junctions_on_minus_are_an_alt5_event(tmp_path):
    # on -, 200 is the acceptor and 300 and 350 the donors
    tables, _ = _events(tmp_path, [(NEAR, 100, 600, 30), (FAR, 100, 600, 10)],
                        [(NEAR, 100, 600, 10), (FAR, 100, 600, 30)], strand='-')
    assert sorted(_of_type(tables, 'alt5')) == ['alt5|chr1:201-350(-)|201-300', 'alt5|chr1:201-350(-)|201-350']
    assert _of_type(tables, 'alt3') == {}


FARTHER = [(200, 380), (400, 500)]


def test_an_alternative_alike_in_every_sample_is_a_row_when_a_large_part_of_the_event(tmp_path):
    # the acceptor at 300, half the reads in both samples, doesn't vary but is a row
    tables, _ = _events(tmp_path, [(NEAR, 100, 600, 50), (FAR, 100, 600, 30), (FARTHER, 100, 600, 20)],
                        [(NEAR, 100, 600, 50), (FAR, 100, 600, 10), (FARTHER, 100, 600, 40)])
    assert sorted(_of_type(tables, 'alt3')) == ['alt3|chr1:201-380(+)|201-300', 'alt3|chr1:201-380(+)|201-350',
                                                'alt3|chr1:201-380(+)|201-380']


def test_alternative_sites_of_too_few_reads_are_not_an_event(tmp_path):
    # the acceptors at 350 and 380, of one read each, vary, but fewer than feature_support
    tables, _ = _events(tmp_path, [(NEAR, 100, 600, 30), (FAR, 100, 600, 1)],
                        [(NEAR, 100, 600, 30), (FARTHER, 100, 600, 1)])
    assert _of_type(tables, 'alt3') == {}


def test_a_ces_inclusion_alike_in_every_sample_is_a_row_when_a_large_part_of_the_event(tmp_path):
    # the first exon included alone, 30% of the reads in both samples, doesn't vary
    first_only = [(200, 300), (400, 700)]
    tables, _ = _events(tmp_path, [(INC2, 100, 800, 40), (first_only, 100, 800, 30), (EXC2, 100, 800, 30)],
                        [(INC2, 100, 800, 20), (first_only, 100, 800, 30), (EXC2, 100, 800, 50)])
    assert _ces_counts(tables) == {'exc': {'s1': '30', 's2': '50'},
                                   '301-400,501-600': {'s1': '40', 's2': '20'},
                                   '301-400': {'s1': '30', 's2': '30'}}


def test_an_intron_inside_another_reads_exon_is_retained(tmp_path):
    # the retaining reads are spliced downstream, as single-exon reads are left out
    spliced, retained = [(200, 300), (400, 500)], [(400, 500)]
    tables, _ = _events(tmp_path, [(spliced, 100, 600, 30), (retained, 100, 600, 10)],
                        [(spliced, 100, 600, 10), (retained, 100, 600, 30)])
    ir = _of_type(tables, 'ir')
    assert sorted(ir) == ['ir|chr1:201-300(+)|retained', 'ir|chr1:201-300(+)|spliced']
    assert ir['ir|chr1:201-300(+)|retained']['s1'] == '10'


def test_an_ir_events_other_row_is_the_junctions_near_its_own(tmp_path):
    # the acceptor at 320 too: its reads don't retain 200-300 either
    spliced, near, retained = [(200, 300), (400, 500)], [(200, 320), (400, 500)], [(400, 500)]
    tables, _ = _events(tmp_path, [(spliced, 100, 600, 30), (near, 100, 600, 10), (retained, 100, 600, 10)],
                        [(spliced, 100, 600, 10), (near, 100, 600, 30), (retained, 100, 600, 10)])
    other = _of_type(tables, 'ir')['ir|chr1:201-300(+)|other']
    assert (other['s1'], other['s2']) == ('10', '30')


def test_a_shifted_copy_of_an_intron_of_few_reads_is_not_retained(tmp_path):
    # 206-306, 2 reads in each sample, is under junc_frac_of_event of the reads over it
    spliced, shifted, retained = [(200, 300), (400, 500)], [(206, 306), (400, 500)], [(400, 500)]
    tables, _ = _events(tmp_path, [(spliced, 100, 600, 30), (shifted, 100, 600, 2), (retained, 100, 600, 10)],
                        [(spliced, 100, 600, 10), (shifted, 100, 600, 2), (retained, 100, 600, 30)])
    ir = _of_type(tables, 'ir')
    assert not any(name.startswith('ir|chr1:207-306(+)') for name in ir)
    assert (ir['ir|chr1:201-300(+)|other']['s1'], ir['ir|chr1:201-300(+)|other']['s2']) == ('2', '2')


def test_a_junction_of_few_reads_nested_in_a_long_intron_is_not_retained(tmp_path):
    # 300-350, 2 reads in each sample, inside the intron 200-500 of 40: under
    # junc_frac_of_event of the reads over it, though its sites are far from the long one's
    long_intron, nested, retained = [(200, 500)], [(300, 350), (400, 500)], [(400, 500)]
    tables, _ = _events(tmp_path, [(long_intron, 100, 600, 40), (nested, 250, 600, 2), (retained, 100, 600, 10)],
                        [(long_intron, 100, 600, 40), (nested, 250, 600, 2), (retained, 100, 600, 30)])
    assert not any(name.startswith('ir|chr1:301-350(+)') for name in _of_type(tables, 'ir'))


def test_first_exons_spliced_to_one_acceptor_are_an_afe_event_not_alt5(tmp_path):
    # first exons ending at 200 and 320, 120 apart, both spliced to 400
    first, second = [(200, 400), (500, 600)], [(320, 400), (500, 600)]
    tables, _ = _events(tmp_path, [(first, 100, 700, 30), (second, 250, 700, 10)],
                        [(first, 100, 700, 10), (second, 250, 700, 30)])
    assert sorted(_of_type(tables, 'afe')) == ['afe|chr1:201-400(+)|201-400', 'afe|chr1:201-400(+)|321-400']
    assert _of_type(tables, 'alt5') == {}


def test_an_internal_exon_whose_reads_are_truncated_in_one_sample_is_not_a_first_exon(tmp_path):
    # exons 100-200, 300-400 and 500-600: s1's reads are whole, while s2's start in
    # the middle exon or stop in it, so in s2 alone the middle exon looks like a
    # first exon.  Pooled, most reads of its junction have an exon before it
    whole, truncated_start, truncated_end = [(200, 300), (400, 500)], [(400, 500)], [(200, 300)]
    tables, _ = _events(tmp_path, [(whole, 100, 600, 30)],
                        [(truncated_end, 100, 380, 20), (truncated_start, 320, 600, 20)])
    assert _of_type(tables, 'afe') == {}


def _afe_with_internal_reads(tmp_path, annot_afe_ss=None):
    # first exons ending at 200 and 320, both spliced to 400, and 4 reads of 44 that
    # have an exon before the one ending at 320, 9%
    first, second, internal = [(200, 400), (500, 600)], [(320, 400), (500, 600)], [(200, 250), (320, 400), (500, 600)]
    return _events(tmp_path, [(first, 100, 700, 30), (second, 250, 700, 10), (internal, 100, 700, 2)],
                   [(first, 100, 700, 10), (second, 250, 700, 30), (internal, 100, 700, 2)], annot_afe_ss=annot_afe_ss)[0]


def test_a_first_exon_also_internal_in_some_reads_is_a_first_exon_only_if_annotated(tmp_path):
    assert _of_type(_afe_with_internal_reads(tmp_path), 'afe') == {}
    afe = _of_type(_afe_with_internal_reads(tmp_path, annot_afe_ss={'+': [320]}), 'afe')
    # 200-250, in 2 of 44 reads in each sample, is too rare to be one of its introns
    assert sorted(row['alternative'] for row in afe.values()) == ['201-400', '321-400']


def test_an_afe_events_reads_are_those_of_its_first_exons(tmp_path):
    # s1's 20 reads starting inside the exon at 400-500, after a first exon in s2's
    # reads, are truncated, not of a first exon: not in the event's total
    first, second, truncated = [(200, 400), (500, 600)], [(320, 400), (500, 600)], [(500, 600)]
    _, totals = _events(tmp_path, [(first, 100, 700, 30), (second, 250, 700, 10), (truncated, 450, 700, 20)],
                        [(first, 100, 700, 10), (second, 250, 700, 30)])
    assert totals[('afe|chr1:201-400(+)', 'event_reads')] == {'s1': 40, 's2': 40}


def test_first_exons_closer_than_min_terminal_junction_separation_are_not_an_afe_event(tmp_path):
    first, second = [(200, 400), (500, 600)], [(260, 400), (500, 600)]
    tables, _ = _events(tmp_path, [(first, 100, 700, 30), (second, 220, 700, 10)],
                        [(first, 100, 700, 10), (second, 220, 700, 30)])
    assert _of_type(tables, 'afe') == {}


def test_an_event_of_few_of_its_locus_reads_is_not_written(tmp_path):
    # the skip's 40 reads per sample under event_frac_of_tot (0.1) of the 440 covering it
    tables, _ = _events(tmp_path, [(INC, 100, 600, 30), (EXC, 100, 600, 10), ([(620, 700)], 150, 800, 400)],
                        [(INC, 100, 600, 10), (EXC, 100, 600, 30), ([(620, 700)], 150, 800, 400)])
    assert _of_type(tables, 'es') == {}
    tables, _ = _events(tmp_path, [(INC, 100, 600, 30), (EXC, 100, 600, 10), ([(620, 700)], 150, 800, 400)],
                        [(INC, 100, 600, 10), (EXC, 100, 600, 30), ([(620, 700)], 150, 800, 400)], event_frac_of_tot=0.05)
    assert _of_type(tables, 'es')


def test_an_exon_included_in_few_of_its_locus_reads_is_written(tmp_path):
    # included in 6 of the event's 100 reads, but under junc_frac_of_event of the 130
    # covering it: an alternative exon, the event 77% of its locus
    tables, _ = _events(tmp_path, [(INC, 100, 600, 6), (EXC, 100, 600, 94), ([(620, 700)], 150, 800, 30)],
                        [(INC, 100, 600, 6), (EXC, 100, 600, 94), ([(620, 700)], 150, 800, 30)])
    assert sorted(_of_type(tables, 'es')) == [f'{ES}|exc', f'{ES}|inc']


def test_last_exons_spliced_from_one_donor_are_an_ale_event_not_alt3(tmp_path):
    # from the donor at 400, last exons starting at 500 and 620
    near, far = [(200, 300), (400, 500)], [(200, 300), (400, 620)]
    tables, _ = _events(tmp_path, [(near, 100, 600, 30), (far, 100, 700, 10)],
                        [(near, 100, 600, 10), (far, 100, 700, 30)])
    assert sorted(_of_type(tables, 'ale')) == ['ale|chr1:401-620(+)|401-500', 'ale|chr1:401-620(+)|401-620']
    assert _of_type(tables, 'alt3') == {}


def test_donors_sharing_alternative_acceptors_are_an_alt5ss_event(tmp_path):
    # donors 200 and 250, each spliced to acceptors 300 and 350, s1 mostly from 200
    def chains(counts):
        return [([(donor, acceptor), (400, 500)], 100, 600, n) for (donor, acceptor), n in counts.items()]
    tables, _ = _events(tmp_path, chains({(200, 300): 40, (200, 350): 20, (250, 300): 20, (250, 350): 10}),
                        chains({(200, 300): 10, (200, 350): 20, (250, 300): 20, (250, 350): 40}))
    alt5ss = _of_type(tables, 'alt5ss')
    # each donor's row by its introns to the shared acceptors
    event = 'alt5ss|chr1:201-350(+)'
    assert sorted(alt5ss) == [f'{event}|201-300,201-350', f'{event}|251-300,251-350']
    assert (alt5ss[f'{event}|201-300,201-350']['s1'], alt5ss[f'{event}|201-300,201-350']['s2']) == ('60', '30')
    assert sorted(_of_type(tables, 'alt3ss')) == ['alt3ss|chr1:201-350(+)|201-300,251-300',
                                                  'alt3ss|chr1:201-350(+)|201-350,251-350']


def test_a_splice_site_in_too_few_reads_in_every_sample_is_counted_as_other(tmp_path):
    # a third acceptor at 380, in 2 reads of 42 in both samples, under junc_frac_of_event
    third = [(200, 380), (400, 500)]
    tables, _ = _events(tmp_path, [(NEAR, 100, 600, 30), (FAR, 100, 600, 10), (third, 100, 600, 2)],
                        [(NEAR, 100, 600, 10), (FAR, 100, 600, 30), (third, 100, 600, 2)])
    alt3 = _of_type(tables, 'alt3')
    assert 'alt3|chr1:201-350(+)|201-380' not in alt3
    other = alt3['alt3|chr1:201-350(+)|other']
    assert (other['s1'], other['s2']) == ('2', '2')


def test_a_chains_first_and_last_exons_and_ends_follow_its_strand():
    two = ((200, 300), (400, 500))
    one = ((600, 700),)
    chains = [{'G1': {two: [ReadEnds(100, 550)] * 2, one: [ReadEnds(580, 800)]}}]
    # on +, the first exon ends at the first junction's start and reads start at their starts
    _, afe, starts = extract_terminal_info(['s1'], chains, 'G1', '+', 'afe')
    assert afe == {200: {(200, 300): {'s1': 2}}, 600: {(600, 700): {'s1': 1}}}
    assert sorted(starts) == [(100, 's1'), (100, 's1'), (580, 's1')]
    # on -, it starts at the last junction's end, also of a chain of one junction, and reads start at their ends
    alljuncs, afe, starts = extract_terminal_info(['s1'], chains, 'G1', '-', 'afe')
    assert afe == {500: {(400, 500): {'s1': 2}}, 700: {(600, 700): {'s1': 1}}}
    assert sorted(starts) == [(550, 's1'), (550, 's1'), (800, 's1')]
    assert alljuncs == {(200, 300): {'s1': 2}, (400, 500): {'s1': 2}, (600, 700): {'s1': 1}}
    _, ale, _ = extract_terminal_info(['s1'], chains, 'G1', '-', 'ale')
    assert ale == {200: {(200, 300): {'s1': 2}}, 600: {(600, 700): {'s1': 1}}}


def test_read_ends_within_the_window_of_the_next_are_counted_at_their_median():
    ends = [(100, 's1'), (150, 's2'), (190, 's1'), (400, 's2')]
    assert group_ends(ends, 100) == {150: {'s1': 2, 's2': 1}, 400: {'s2': 1}}
    assert group_ends([], 100) == {}


def _annotated(gene_introns, strand='+', num_transcripts=None):
    "AnnotatedSplicing of genes' introns, {gene: [(start, end)]}"
    annotated = AnnotatedSplicing()
    annotated.num_transcripts = num_transcripts or {gene: 1 for gene in gene_introns}
    for gene, introns in gene_introns.items():
        for junc in introns:
            annotated.intron_genes[strand].setdefault(junc, set()).add(gene)
            for site in junc:
                annotated.site_genes[strand].setdefault(site, set()).add(gene)
    return annotated


ES_READS = ([(INC, 100, 600, 30), (EXC, 100, 600, 10)], [(INC, 100, 600, 10), (EXC, 100, 600, 30)])


def test_a_splicing_event_is_of_the_gene_annotated_with_most_of_its_introns_and_sites(tmp_path):
    # the reads were given G1, but GX has both inclusion introns, GY only the first
    annotated = _annotated({'GX': INC, 'GY': [(200, 300)]})
    tables, totals = _events(tmp_path, *ES_READS, annotated=annotated)
    assert {row['gene'] for row in _of_type(tables, 'es').values()} == {'GX'}
    assert totals[(ES, 'event_reads')] == {'s1': 40, 's2': 40}


def test_an_events_gene_name_is_the_annotations_empty_if_it_names_none(tmp_path):
    annotated = _annotated({'GX': INC})
    annotated.names = {'GX': 'GENEX'}
    tables, totals = _events(tmp_path, *ES_READS, annotated=annotated)
    assert {row['gene_name'] for row in _of_type(tables, 'es').values()} == {'GENEX'}
    annotated.names = {}
    tables, _ = _events(tmp_path, *ES_READS, annotated=annotated)
    assert {row['gene_name'] for row in _of_type(tables, 'es').values()} == {''}


def test_the_other_genes_sharing_an_events_splice_sites_are_its_additional_genes(tmp_path):
    # GX has both inclusion introns, 6 with the skipping junction's sites; GY the first
    # inclusion intron, 3; GZ the skipping junction's acceptor, 1
    annotated = _annotated({'GX': INC, 'GY': [(200, 300)], 'GZ': [(450, 500)]})
    annotated.names = {'GX': 'GENEX', 'GZ': 'GENEZ'}
    tables, totals = _events(tmp_path, *ES_READS, annotated=annotated)
    for row in _of_type(tables, 'es').values():
        assert (row['gene'], row['gene_name']) == ('GX', 'GENEX')
        assert (row['additional_genes'], row['additional_gene_names']) == ('GY,GZ', ',GENEZ')


def test_genes_annotated_alike_with_an_event_are_all_its_genes(tmp_path):
    tables, _ = _events(tmp_path, *ES_READS, annotated=_annotated({'GX': INC, 'GY': INC}))
    assert {row['gene'] for row in _of_type(tables, 'es').values()} == {'GX,GY'}


def test_a_tie_goes_to_the_gene_of_more_annotated_transcripts(tmp_path):
    # as a gene and a readthrough of it, annotated with the same introns
    annotated = _annotated({'GX': INC, 'GX-GZ': INC}, num_transcripts={'GX': 5, 'GX-GZ': 1})
    tables, _ = _events(tmp_path, *ES_READS, annotated=annotated)
    assert {row['gene'] for row in _of_type(tables, 'es').values()} == {'GX'}


def test_an_annotated_intron_counts_2_and_each_splice_site_of_another_1(tmp_path):
    # the event's introns are the inclusion introns 200-300 and 400-500 and the skipping
    # junction 200-500.  GY has 200-300, 2, and of 200-500 the site 200, 1: 3.  GX has
    # sites 200 and 300 of 200-300, 400 of 400-500, and 200 of 200-500: 4
    annotated = _annotated({'GY': [(200, 300)], 'GX': [(150, 300), (400, 450), (180, 200)]})
    tables, _ = _events(tmp_path, *ES_READS, annotated=annotated)
    assert {row['gene'] for row in _of_type(tables, 'es').values()} == {'GX'}


def test_an_event_of_no_annotated_intron_or_site_is_of_its_reads_gene(tmp_path):
    tables, _ = _events(tmp_path, *ES_READS, annotated=_annotated({'GX': [(1000, 1100)]}))
    assert {row['gene'] for row in _of_type(tables, 'es').values()} == {'G1'}


def test_an_event_of_no_annotation_is_of_the_gene_of_most_of_its_reads():
    chain_genes = {('+', ((200, 300),)): ('A',), ('+', ((200, 300), (400, 500))): ('A', 'B'),
                   ('+', ((400, 500),)): ('B',), ('+', ((700, 800),)): ('C',)}
    cluster_chains = [{((200, 300),): [ReadEnds(100, 350)] * 2, ((200, 300), (400, 500)): [ReadEnds(100, 550)],
                       ((400, 500),): [ReadEnds(350, 550)] * 2},
                      {((700, 800),): [ReadEnds(650, 850)] * 10}]
    # of the intron 200-300, A has 3 reads and B 1, the readthrough's counted for both
    assert read_event_genes({(200, 300)}, cluster_chains, chain_genes, '+') == ['A']
    # of both introns, 3 each: a tie; C's reads carry neither
    assert read_event_genes({(200, 300), (400, 500)}, cluster_chains, chain_genes, '+') == ['A', 'B']


def test_a_rows_introns_not_in_the_annotation_are_its_novel_junctions(tmp_path):
    # the inclusion introns annotated, the skipping junction not
    tables, _ = _events(tmp_path, *ES_READS, annotated=_annotated({'GX': INC}))
    assert tables['counts'][f'{ES}|inc']['novel_junctions'] == ''
    assert tables['counts'][f'{ES}|exc']['novel_junctions'] == 'chr1:201-500'
    # an alternative acceptor's novel intron
    tables, _ = _events(tmp_path, [([(200, 300)], 100, 400, 30), ([(200, 350)], 100, 400, 10)],
                        [([(200, 300)], 100, 400, 10), ([(200, 350)], 100, 400, 30)],
                        annotated=_annotated({'GX': [(200, 300)]}))
    alt3 = _of_type(tables, 'alt3')
    assert alt3['alt3|chr1:201-350(+)|201-300']['novel_junctions'] == ''
    assert alt3['alt3|chr1:201-350(+)|201-350']['novel_junctions'] == 'chr1:201-350'


def test_an_events_gene_counts_the_introns_of_its_other_row(tmp_path):
    # alt3 from the donor at 200 to 300 and 350, varying, and to 380 alike in both
    # samples and under MIN_STEADY_FRACTION of the reads, in its other row.  Of 200-300
    # and 200-350, GX and GY score 3 each; with 200-380, which only GX has, GX 5 and GY 4
    s1 = [([(200, 300)], 100, 400, 30), ([(200, 350)], 100, 400, 10), ([(200, 380)], 100, 400, 3)]
    s2 = [([(200, 300)], 100, 400, 10), ([(200, 350)], 100, 400, 30), ([(200, 380)], 100, 400, 3)]
    annotated = _annotated({'GX': [(200, 300), (200, 380)], 'GY': [(200, 350)]})
    tables, _ = _events(tmp_path, s1, s2, junc_frac_of_event=0.2, annotated=annotated)
    alt3 = _of_type(tables, 'alt3')
    assert 'alt3|chr1:201-350(+)|other' in alt3
    assert {row['gene'] for row in alt3.values()} == {'GX'}


def test_an_events_introns_count_for_its_gene_by_their_reads(tmp_path):
    # alt3 to 300, which GX has, in 60 reads, and to 350, which GY has, in 20: each
    # gene has one intron and the other's donor, a tie but for the reads
    s1 = [([(200, 300)], 100, 400, 40), ([(200, 350)], 100, 400, 10)]
    s2 = [([(200, 300)], 100, 400, 20), ([(200, 350)], 100, 400, 10)]
    annotated = _annotated({'GX': [(200, 300)], 'GY': [(200, 350)]})
    tables, _ = _events(tmp_path, s1, s2, annotated=annotated)
    assert {row['gene'] for row in _of_type(tables, 'alt3').values()} == {'GX'}
