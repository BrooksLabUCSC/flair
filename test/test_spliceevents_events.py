"""
flair spliceevents' events, from each sample's corrected reads: which events are
found, and their rows' counts, PSIs and totals.

Reads are given as (junctions, start, end, number of reads), on one gene G1 on chr1,
with 0-based starts, and the tables' rows are keyed by name, event|alternative, 1-based.
"""
import csv
from flair import SeqRange
from flair.isoform_data import ReadRec
from flair.flair_spliceevents import process_gene_to_events, psi_varies, EVENT_ROW_COLUMNS

SAMPLES = ['s1', 's2', 's3']
CHROM = 'chr1'


def _sample_reads(chains, gene='G1', strand='+'):
    "one sample's {gene: {junctions: [ReadRec]}}"
    gene_to_juncs = {gene: {}}
    n = 0
    for juncs, start, end, count in chains:
        for _ in range(count):
            n += 1
            gene_to_juncs[gene].setdefault(tuple(juncs), []).append(ReadRec(None, strand, (), start, end, f'r{n}'))
    return gene_to_juncs


def _events(tmp_path, *samples_chains, strand='+', feature_support=2, event_frac_of_tot=0.1,
            junc_frac_of_event=0.02, event_support=20, annot_afe_ss=None, annot_ale_ss=None):
    """process_gene_to_events on the reads of each sample: the rows of each table, as
    {row name: {column: value}}, keyed by table (counts, PSIjunc, PSItot), and the
    totals, as {(event, total): {sample: count}}"""
    samples = SAMPLES[:len(samples_chains)]
    prefix = str(tmp_path / 'region')
    process_gene_to_events(prefix, CHROM, samples, [_sample_reads(chains, strand=strand) for chains in samples_chains],
                           {'G1': strand}, feature_support, False, event_frac_of_tot, junc_frac_of_event, event_support,
                           annot_afe_ss or {}, annot_ale_ss or {}, [SeqRange(CHROM, 0, 10 ** 9)], {})
    columns = EVENT_ROW_COLUMNS + samples
    tables = {}
    for table in ('counts', 'PSIjunc', 'PSItot'):
        with open(f'{prefix}.spliceevents.{table}.tsv') as fh:
            tables[table] = {f"{row['event']}|{row['alternative']}": row
                             for row in csv.DictReader(fh, fieldnames=columns, delimiter='\t')}
    totals = {}
    with open(f'{prefix}.spliceevents.totals.tsv') as fh:
        for line in fh:
            event, _, _, total, *counts = line.rstrip('\n').split('\t')
            totals[(event, total)] = dict(zip(samples, (int(c) for c in counts)))
    return tables, totals


def _of_type(tables, eventtype):
    return {name: row for name, row in tables['counts'].items() if row['event_type'] == eventtype}


# exons 100-200, 300-400 and 500-600 on +: the middle one included or skipped
INC = [(200, 300), (400, 500)]
EXC = [(200, 500)]
ES = 'es|G1|chr1:301-400(+)'


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


def test_psis_vary_by_junc_frac_of_event_or_one_psi_is_away_from_0_and_1():
    assert psi_varies([0.5, 0.53], 0.02)
    assert not psi_varies([0.5, 0.51], 0.02)
    assert psi_varies([0.5], 0.02)
    assert not psi_varies([0.99], 0.02)
    assert not psi_varies([], 0.02)


def test_a_skip_in_fewer_reads_than_feature_support_is_not_an_event(tmp_path):
    tables, _ = _events(tmp_path, [(INC, 100, 600, 30), (EXC, 100, 600, 1)],
                        [(INC, 100, 600, 30)])
    assert _of_type(tables, 'es') == {}


def test_an_exon_skipped_alike_in_every_sample_is_not_an_event(tmp_path):
    # PSIs differ by less than junc_frac_of_event
    tables, _ = _events(tmp_path, [(INC, 100, 600, 30), (EXC, 100, 600, 10)],
                        [(INC, 100, 600, 30), (EXC, 100, 600, 10)], junc_frac_of_event=0.05)
    assert _of_type(tables, 'es') == {}


def test_one_sample_needs_its_psi_between_junc_frac_of_event_and_one_minus_it(tmp_path):
    assert _of_type(_events(tmp_path, [(INC, 100, 600, 30), (EXC, 100, 600, 10)])[0], 'es')
    assert not _of_type(_events(tmp_path, [(INC, 100, 600, 300), (EXC, 100, 600, 2)], junc_frac_of_event=0.05)[0], 'es')


# exons 100-200, 300-400, 500-600 and 700-800 on +: the two middle ones included or skipped together
INC2 = [(200, 300), (400, 500), (600, 700)]
EXC2 = [(200, 700)]
CES = 'ces|G1|chr1:201-700(+)'


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
    assert sorted(alt3) == ['alt3|G1|chr1:201-350(+)|201-300', 'alt3|G1|chr1:201-350(+)|201-350']
    assert alt3['alt3|G1|chr1:201-350(+)|201-350']['s2'] == '30'
    assert _of_type(tables, 'alt5') == {}


def test_the_same_junctions_on_minus_are_an_alt5_event(tmp_path):
    # on -, 200 is the acceptor and 300 and 350 the donors
    tables, _ = _events(tmp_path, [(NEAR, 100, 600, 30), (FAR, 100, 600, 10)],
                        [(NEAR, 100, 600, 10), (FAR, 100, 600, 30)], strand='-')
    assert sorted(_of_type(tables, 'alt5')) == ['alt5|G1|chr1:201-350(-)|201-300', 'alt5|G1|chr1:201-350(-)|201-350']
    assert _of_type(tables, 'alt3') == {}


def test_an_intron_inside_another_reads_exon_is_retained(tmp_path):
    # the retaining reads are spliced downstream, as single-exon reads are left out
    spliced, retained = [(200, 300), (400, 500)], [(400, 500)]
    tables, _ = _events(tmp_path, [(spliced, 100, 600, 30), (retained, 100, 600, 10)],
                        [(spliced, 100, 600, 10), (retained, 100, 600, 30)])
    ir = _of_type(tables, 'ir')
    assert sorted(ir) == ['ir|G1|chr1:201-300(+)|retained', 'ir|G1|chr1:201-300(+)|spliced']
    assert ir['ir|G1|chr1:201-300(+)|retained']['s1'] == '10'


def test_first_exons_spliced_to_one_acceptor_are_an_afe_event_not_alt5(tmp_path):
    # first exons ending at 200 and 320, 120 apart, both spliced to 400
    first, second = [(200, 400), (500, 600)], [(320, 400), (500, 600)]
    tables, _ = _events(tmp_path, [(first, 100, 700, 30), (second, 250, 700, 10)],
                        [(first, 100, 700, 10), (second, 250, 700, 30)])
    assert sorted(_of_type(tables, 'afe')) == ['afe|G1|chr1:201-400(+)|201-400', 'afe|G1|chr1:201-400(+)|321-400']
    assert _of_type(tables, 'alt5') == {}


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


def test_last_exons_spliced_from_one_donor_are_an_ale_event_not_alt3(tmp_path):
    # from the donor at 400, last exons starting at 500 and 620
    near, far = [(200, 300), (400, 500)], [(200, 300), (400, 620)]
    tables, _ = _events(tmp_path, [(near, 100, 600, 30), (far, 100, 700, 10)],
                        [(near, 100, 600, 10), (far, 100, 700, 30)])
    assert sorted(_of_type(tables, 'ale')) == ['ale|G1|chr1:401-620(+)|401-500', 'ale|G1|chr1:401-620(+)|401-620']
    assert _of_type(tables, 'alt3') == {}


def test_donors_sharing_alternative_acceptors_are_an_alt5ss_event(tmp_path):
    # donors 200 and 250, each spliced to acceptors 300 and 350, s1 mostly from 200
    def chains(counts):
        return [([(donor, acceptor), (400, 500)], 100, 600, n) for (donor, acceptor), n in counts.items()]
    tables, _ = _events(tmp_path, chains({(200, 300): 40, (200, 350): 20, (250, 300): 20, (250, 350): 10}),
                        chains({(200, 300): 10, (200, 350): 20, (250, 300): 20, (250, 350): 40}))
    alt5ss = _of_type(tables, 'alt5ss')
    # each donor's row by its introns to the shared acceptors
    event = 'alt5ss|G1|chr1:201-350(+)'
    assert sorted(alt5ss) == [f'{event}|201-300,201-350', f'{event}|251-300,251-350']
    assert (alt5ss[f'{event}|201-300,201-350']['s1'], alt5ss[f'{event}|201-300,201-350']['s2']) == ('60', '30')
    assert sorted(_of_type(tables, 'alt3ss')) == ['alt3ss|G1|chr1:201-350(+)|201-300,251-300',
                                                  'alt3ss|G1|chr1:201-350(+)|201-350,251-350']


def test_a_splice_site_used_alike_in_every_sample_is_counted_as_other(tmp_path):
    # a third acceptor at 380, in 4 reads of 44 in both samples
    third = [(200, 380), (400, 500)]
    tables, _ = _events(tmp_path, [(NEAR, 100, 600, 30), (FAR, 100, 600, 10), (third, 100, 600, 4)],
                        [(NEAR, 100, 600, 10), (FAR, 100, 600, 30), (third, 100, 600, 4)])
    alt3 = _of_type(tables, 'alt3')
    assert 'alt3|G1|chr1:201-350(+)|201-380' not in alt3
    other = alt3['alt3|G1|chr1:201-350(+)|other']
    assert (other['s1'], other['s2']) == ('4', '4')
