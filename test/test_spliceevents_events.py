"""
flair spliceevents' events, from each sample's corrected reads: which events are
found, and their rows' counts, PSIs and totals.

Reads are given as (junctions, start, end, number of reads), on one gene G1 on chr1.
"""
import csv
from flair import SeqRange
from flair.isoform_data import ReadRec
from flair.flair_spliceevents import process_gene_to_events, psi_varies

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
    columns = ['eventname', 'event', 'eventtype', 'gene', 'junctions_included', 'junctions_excluded',
               'outer_junctions', 'exons'] + samples
    tables = {}
    for table in ('counts', 'PSIjunc', 'PSItot'):
        with open(f'{prefix}.spliceevents.{table}.tsv') as fh:
            tables[table] = {row['eventname']: row for row in csv.DictReader(fh, fieldnames=columns, delimiter='\t')}
    totals = {}
    with open(f'{prefix}.spliceevents.totals.tsv') as fh:
        for line in fh:
            event, _, _, total, *counts = line.rstrip('\n').split('\t')
            totals[(event, total)] = dict(zip(samples, (int(c) for c in counts)))
    return tables, totals


def _of_type(tables, eventtype):
    return {name: row for name, row in tables['counts'].items() if row['eventtype'] == eventtype}


# exons 100-200, 300-400 and 500-600 on +: the middle one included or skipped
INC = [(200, 300), (400, 500)]
EXC = [(200, 500)]
ES = 'es-of-chr1:300-400'


def test_a_skipped_exon_is_an_es_event_with_its_inclusion_and_exclusion(tmp_path):
    tables, totals = _events(tmp_path, [(INC, 100, 600, 30), (EXC, 100, 600, 10)],
                             [(INC, 100, 600, 10), (EXC, 100, 600, 30)])
    es = _of_type(tables, 'es')
    assert sorted(es) == [f'exc_{ES}', f'inc_{ES}']
    inc, exc = tables['counts'][f'inc_{ES}'], tables['counts'][f'exc_{ES}']
    assert (inc['s1'], inc['s2'], exc['s1'], exc['s2']) == ('30', '10', '10', '30')
    assert inc['exons'] == 'chr1:300-400' and exc['outer_junctions'] == 'chr1:200-500'
    assert (tables['PSIjunc'][f'inc_{ES}']['s1'], tables['PSIjunc'][f'inc_{ES}']['s2']) == ('0.75', '0.25')
    assert totals[(ES, 'event_reads')] == {'s1': 40, 's2': 40}


def test_psitot_divides_by_all_reads_covering_the_event(tmp_path):
    # 40 more reads in s1 across the exon's locus, spliced only downstream of it: reads
    # of the locus, not of the event
    tables, totals = _events(tmp_path, [(INC, 100, 600, 30), (EXC, 100, 600, 10), ([(620, 700)], 150, 800, 40)],
                             [(INC, 100, 600, 10), (EXC, 100, 600, 30)])
    assert totals[(ES, 'locus_reads')]['s1'] == 80
    assert tables['PSItot'][f'exc_{ES}']['s1'] == '0.125'
    assert tables['PSIjunc'][f'exc_{ES}']['s1'] == '0.25'


def test_a_sample_with_fewer_event_reads_than_event_support_has_no_psi(tmp_path):
    tables, _ = _events(tmp_path, [(INC, 100, 600, 30), (EXC, 100, 600, 10)],
                        [(INC, 100, 600, 10), (EXC, 100, 600, 30)],
                        [(INC, 100, 600, 5), (EXC, 100, 600, 14)])
    assert tables['counts'][f'inc_{ES}']['s3'] == '5'
    assert tables['PSIjunc'][f'inc_{ES}']['s3'] == ''


def test_an_event_with_a_psi_in_only_one_sample_is_reported(tmp_path):
    # s2 under event_support leaves s1's PSI, 0.75, alone: an alternative exon, as in one sample
    tables, _ = _events(tmp_path, [(INC, 100, 600, 30), (EXC, 100, 600, 10)],
                        [(INC, 100, 600, 5), (EXC, 100, 600, 14)])
    assert tables['PSIjunc'][f'inc_{ES}']['s1'] == '0.75'
    assert tables['PSIjunc'][f'inc_{ES}']['s2'] == ''


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
CES = 'ces-relTo-chr1:200-700(+)-G1'


def _ces_counts(tables):
    rows = _of_type(tables, 'ces')
    return {name.split('_ces')[0]: {s: row[s] for s in SAMPLES[:2]} for name, row in rows.items()}


def test_two_exons_skipped_together_are_a_ces_event(tmp_path):
    tables, _ = _events(tmp_path, [(INC2, 100, 800, 30), (EXC2, 100, 800, 10)],
                        [(INC2, 100, 800, 10), (EXC2, 100, 800, 30)])
    assert _ces_counts(tables) == {'exc': {'s1': '10', 's2': '30'},
                                   'inc-of-chr1:300-400,chr1:500-600': {'s1': '30', 's2': '10'}}


def test_a_ces_events_counts_do_not_depend_on_the_order_of_its_reads(tmp_path):
    # the skipping reads first, as when the bam has them first
    tables, _ = _events(tmp_path, [(EXC2, 100, 800, 10), (INC2, 100, 800, 30)],
                        [(EXC2, 100, 800, 30), (INC2, 100, 800, 10)])
    assert _ces_counts(tables) == {'exc': {'s1': '10', 's2': '30'},
                                   'inc-of-chr1:300-400,chr1:500-600': {'s1': '30', 's2': '10'}}


def test_a_ces_events_inclusion_counts_the_samples_without_skipping(tmp_path):
    # s1 has no skipping reads, s2's come after s1 is counted
    tables, _ = _events(tmp_path, [(INC2, 100, 800, 30)],
                        [(INC2, 100, 800, 10), (EXC2, 100, 800, 30)])
    assert _ces_counts(tables)['inc-of-chr1:300-400,chr1:500-600'] == {'s1': '30', 's2': '10'}


# from a donor at 200, acceptors at 300 or 350, then on to an exon at 500-600
NEAR = [(200, 300), (400, 500)]
FAR = [(200, 350), (400, 500)]


def test_two_acceptors_of_one_donor_are_an_alt3_event_on_plus(tmp_path):
    tables, _ = _events(tmp_path, [(NEAR, 100, 600, 30), (FAR, 100, 600, 10)],
                        [(NEAR, 100, 600, 10), (FAR, 100, 600, 30)])
    alt3 = _of_type(tables, 'alt3')
    assert sorted(alt3) == ['chr1:200-300_alt3-relTo-chr1:200(+)-G1', 'chr1:200-350_alt3-relTo-chr1:200(+)-G1']
    assert alt3['chr1:200-350_alt3-relTo-chr1:200(+)-G1']['s2'] == '30'
    assert _of_type(tables, 'alt5') == {}


def test_the_same_junctions_on_minus_are_an_alt5_event(tmp_path):
    # on -, 200 is the acceptor and 300 and 350 the donors
    tables, _ = _events(tmp_path, [(NEAR, 100, 600, 30), (FAR, 100, 600, 10)],
                        [(NEAR, 100, 600, 10), (FAR, 100, 600, 30)], strand='-')
    assert sorted(_of_type(tables, 'alt5')) == ['chr1:200-300_alt5-relTo-chr1:200(-)-G1',
                                                'chr1:200-350_alt5-relTo-chr1:200(-)-G1']
    assert _of_type(tables, 'alt3') == {}


def test_an_intron_inside_another_reads_exon_is_retained(tmp_path):
    # the retaining reads are spliced downstream, as single-exon reads are left out
    spliced, retained = [(200, 300), (400, 500)], [(400, 500)]
    tables, _ = _events(tmp_path, [(spliced, 100, 600, 30), (retained, 100, 600, 10)],
                        [(spliced, 100, 600, 10), (retained, 100, 600, 30)])
    ir = _of_type(tables, 'ir')
    assert sorted(ir) == ['retained_ir-of-chr1:200-300(+)-G1', 'spliced_ir-of-chr1:200-300(+)-G1']
    assert ir['retained_ir-of-chr1:200-300(+)-G1']['s1'] == '10'


def test_first_exons_spliced_to_one_acceptor_are_an_afe_event_not_alt5(tmp_path):
    # first exons ending at 200 and 320, 120 apart, both spliced to 400
    first, second = [(200, 400), (500, 600)], [(320, 400), (500, 600)]
    tables, _ = _events(tmp_path, [(first, 100, 700, 30), (second, 250, 700, 10)],
                        [(first, 100, 700, 10), (second, 250, 700, 30)])
    assert sorted(_of_type(tables, 'AFE')) == ['chr1:200_AFE-(+)-G1', 'chr1:320_AFE-(+)-G1']
    assert _of_type(tables, 'alt5') == {}


def test_first_exons_closer_than_min_terminal_junction_separation_are_not_an_afe_event(tmp_path):
    first, second = [(200, 400), (500, 600)], [(260, 400), (500, 600)]
    tables, _ = _events(tmp_path, [(first, 100, 700, 30), (second, 220, 700, 10)],
                        [(first, 100, 700, 10), (second, 220, 700, 30)])
    assert _of_type(tables, 'AFE') == {}


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
    assert sorted(_of_type(tables, 'ALE')) == ['chr1:500_ALE-(+)-G1', 'chr1:620_ALE-(+)-G1']
    assert _of_type(tables, 'alt3') == {}


def test_donors_sharing_alternative_acceptors_are_an_alt5ss_event(tmp_path):
    # donors 200 and 250, each spliced to acceptors 300 and 350, s1 mostly from 200
    def chains(counts):
        return [([(donor, acceptor), (400, 500)], 100, 600, n) for (donor, acceptor), n in counts.items()]
    tables, _ = _events(tmp_path, chains({(200, 300): 40, (200, 350): 20, (250, 300): 20, (250, 350): 10}),
                        chains({(200, 300): 10, (200, 350): 20, (250, 300): 20, (250, 350): 40}))
    alt5ss = _of_type(tables, 'alt5ss')
    event = 'alt5ss-relTo-chr1:300,chr1:350(+)-G1'
    assert sorted(alt5ss) == [f'chr1:200_{event}', f'chr1:250_{event}']
    assert (alt5ss[f'chr1:200_{event}']['s1'], alt5ss[f'chr1:200_{event}']['s2']) == ('60', '30')
    assert sorted(_of_type(tables, 'alt3ss')) == ['chr1:300_alt3ss-relTo-chr1:200,chr1:250(+)-G1',
                                                  'chr1:350_alt3ss-relTo-chr1:200,chr1:250(+)-G1']


def test_a_splice_site_used_alike_in_every_sample_is_counted_as_other(tmp_path):
    # a third acceptor at 380, in 4 reads of 44 in both samples
    third = [(200, 380), (400, 500)]
    tables, _ = _events(tmp_path, [(NEAR, 100, 600, 30), (FAR, 100, 600, 10), (third, 100, 600, 4)],
                        [(NEAR, 100, 600, 10), (FAR, 100, 600, 30), (third, 100, 600, 4)])
    alt3 = _of_type(tables, 'alt3')
    assert 'chr1:200-380_alt3-relTo-chr1:200(+)-G1' not in alt3
    other = alt3['other_alt3-relTo-chr1:200(+)-G1']
    assert (other['s1'], other['s2']) == ('4', '4')
