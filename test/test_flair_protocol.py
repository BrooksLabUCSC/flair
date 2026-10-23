"""
misc/flair_protocol.py, the protocol from one manifest: its checks and the steps it
runs.  A prototype, so these run only with FLAIR_PROTOCOL_TESTS=1 (make
test-protocol), not in the default suite.
"""
import importlib.util
import io
import os
import pytest

pytestmark = pytest.mark.skipif(os.environ.get('FLAIR_PROTOCOL_TESTS') != '1',
                                reason='the protocol script is a prototype: set FLAIR_PROTOCOL_TESTS=1')

_spec = importlib.util.spec_from_file_location(
    'flair_protocol', os.path.join(os.path.dirname(__file__), '..', 'misc', 'flair_protocol.py'))
protocol = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(protocol)


def _manifest(tmp_path, columns, rows):
    "a manifest of rows, each sample's bam (indexed) and vcf made empty files"
    lines = ['\t'.join(columns)]
    for row in rows:
        values = dict(zip(columns, row))
        for col in ('bam', 'vcf'):
            if col in values:
                path = tmp_path / values[col]
                path.touch()
                if col == 'bam':
                    (tmp_path / (values[col] + '.bai')).touch()
                values[col] = str(path)
        lines.append('\t'.join(values[c] for c in columns))
    path = tmp_path / 'manifest.tsv'
    path.write_text('\n'.join(lines) + '\n')
    return str(path)


def _args(tmp_path, manifest, *options):
    for name in ('genome.fa', 'genes.gtf'):
        (tmp_path / name).touch()
    return protocol.parse_args(['-m', manifest, '-g', str(tmp_path / 'genome.fa'), '-f', str(tmp_path / 'genes.gtf'),
                                '-o', str(tmp_path / 'out'), *options])


def _plan(tmp_path, manifest, *options):
    "the steps of a run, as {name: step}, in order"
    args = _args(tmp_path, manifest, *options)
    samples = protocol.read_manifest(args.manifest)
    patients = protocol.check(args, samples, which=lambda tool: tool)
    return args, {step.name: step for step in protocol.build_steps(args, samples, patients)}


CONDITIONS = [(f'{c}{n}', f'{c}{n}.bam', c) for c in ('ctl', 'trt') for n in (1, 2, 3)]
PATIENTS = [('P1N', 'P1N.bam', 'P1', 'normal', 'P1N.vcf'), ('P1T', 'P1T.bam', 'P1', 'tumor', 'P1T.vcf'),
            ('P2T', 'P2T.bam', 'P2', 'tumor', 'P2T.vcf')]
PATIENT_COLUMNS = ('sample', 'bam', 'patient', 'tissue', 'vcf')


def test_the_base_run_is_transcriptome_combine_quantify_with_junctions_from_the_reads(tmp_path):
    manifest = _manifest(tmp_path, ('sample', 'bam'), [('s1', 's1.bam'), ('s2', 's2.bam')])
    args, steps = _plan(tmp_path, manifest)
    assert list(steps) == ['junctions.s1', 'junctions.s2', 'junctions.pool', 'transcriptome.s1', 'transcriptome.s2',
                           'combine.manifest', 'combine', 'quantify.manifest', 'quantify']
    cmd = steps['transcriptome.s1'].cmd
    assert cmd[cmd.index('--sample_name') + 1] == 's1'
    assert cmd[cmd.index('--junction_bed') + 1] == os.path.join(args.output, 'junctions', 'all.bed')
    assert '--normalize_ends' in steps['quantify'].cmd


def test_a_given_junction_file_is_used_and_intron_prospector_is_not_run(tmp_path):
    manifest = _manifest(tmp_path, ('sample', 'bam'), [('s1', 's1.bam')])
    (tmp_path / 'SJ.out.tab').touch()
    _, steps = _plan(tmp_path, manifest, '--junction_star', str(tmp_path / 'SJ.out.tab'))
    assert not any(name.startswith('junctions') for name in steps)
    assert '--junction_star' in steps['transcriptome.s1'].cmd


def test_fusion_isoforms_are_combined_with_the_rest(tmp_path):
    manifest = _manifest(tmp_path, ('sample', 'bam'), [('s1', 's1.bam')])
    args, steps = _plan(tmp_path, manifest, '--fusion', '--junction_bed', str(tmp_path / 'genes.gtf'))
    assert 'fusion.s1' in steps
    steps['combine.manifest'].write()
    beds = open(os.path.join(args.output, 'combine', 'manifest.txt')).read().split()
    assert [os.path.basename(b) for b in beds] == ['s1.isoforms.bed', 's1.fusions.isoforms.bed']


def test_conditions_run_diffexp_on_isoform_and_event_counts(tmp_path):
    manifest = _manifest(tmp_path, ('sample', 'bam', 'condition'), CONDITIONS)
    args, steps = _plan(tmp_path, manifest, '--spliceevents', '--stat_test', 'conditions')
    assert [n for n in steps if n.startswith('diffexp')] == ['diffexp.isoforms', 'diffexp.spliceevents']
    cmd = steps['diffexp.isoforms'].cmd
    assert cmd[cmd.index('--condition_a') + 1:cmd.index('--condition_a') + 4:2] == ['ctl', 'trt']
    # diffexp reads the event counts' conditions from a sample info file beside them
    steps['spliceevents.sample_info'].write()
    info = open(os.path.join(args.output, 'spliceevents', 'flair.spliceevents.sample_info.tsv')).read().splitlines()
    assert info[:2] == ['sample_id\tcondition\tbatch', 'ctl1\tctl\tb1']


def test_too_few_replicates_for_diffexp_is_an_error(tmp_path):
    manifest = _manifest(tmp_path, ('sample', 'bam', 'condition'), CONDITIONS[:5])
    with pytest.raises(protocol.ProtocolError, match='at least 3'):
        _plan(tmp_path, manifest, '--stat_test', 'conditions')


def test_patients_phase_pairs_and_tumor_only_samples_and_test_only_pairs(tmp_path):
    manifest = _manifest(tmp_path, PATIENT_COLUMNS, PATIENTS)
    args, steps = _plan(tmp_path, manifest, '--spliceevents', '--stat_test', 'patients', '--haplotype')
    # the pair: the normal's variants phased with both samples' reads, both haplotagged
    phase = steps['phase.P1.phased.vcf.gz'].cmd
    assert phase[-3:] == [os.path.join(str(tmp_path), 'P1N.vcf'), os.path.join(str(tmp_path), 'P1N.bam'),
                          os.path.join(str(tmp_path), 'P1T.bam')]
    assert {'haplotag.P1N', 'haplotag.P1T', 'haplotag.P2T'} <= set(steps)
    pair = steps['isoalleles.P1T'].cmd
    assert pair[pair.index('--normal_name') + 1] == 'P1N' and '--tumor_vcf' in pair
    # the tumor-only patient: its own variants, one sample, no somatic variants
    single = steps['isoalleles.P2T'].cmd
    assert '--normal_bam' not in single and '--tumor_vcf' not in single
    assert '--generate_map' in steps['quantify'].cmd
    # only the pair is tested, on every counts table
    usage = [n for n in steps if n.startswith('usage')]
    assert usage == ['usage.P1T.isoforms', 'usage.P1T.spliceevents', 'usage.P1T.isoalleles']
    assert steps['usage.P1T.isoforms'].cmd[2:4] == ['P1T', 'P1N']
    report = protocol.summary(args, protocol.patients_of(protocol.read_manifest(args.manifest)))
    assert 'tumor-only patients, reported but not tested: P2' in report


@pytest.mark.parametrize('columns, rows, options, message', [
    (('sample', 'bam'), [('s1', 's1.bam'), ('s1', 's2.bam')], [], 'unique'),
    (('sample', 'bam'), [('s1', 's1.bam')], ['--stat_test', 'patients'], 'patient column'),
    (('sample', 'bam', 'condition'), CONDITIONS, ['--stat_test', 'conditions', '--haplotype'], 'vcf column'),
    (PATIENT_COLUMNS, [('a', 'a.bam', 'P1', 'normal', 'a.vcf'), ('b', 'b.bam', 'P1', 'normal', 'b.vcf')],
     ['--stat_test', 'patients'], 'two normals'),
    (PATIENT_COLUMNS, PATIENTS, ['--haplotype'], '--stat_test patients'),
])
def test_manifest_and_option_errors(tmp_path, columns, rows, options, message):
    manifest = _manifest(tmp_path, columns, rows)
    with pytest.raises(protocol.ProtocolError, match=message):
        _plan(tmp_path, manifest, *options)


def test_a_rerun_skips_the_steps_already_done(tmp_path):
    manifest = _manifest(tmp_path, ('sample', 'bam'), [('s1', 's1.bam')])
    args = _args(tmp_path, manifest)
    out = tmp_path / 'out' / 'x.txt'
    steps = [protocol.Step('write', None, [str(out)], protocol._write_lines(str(out), ['a']))]
    log = io.StringIO()
    protocol.run(args, steps, out=log)
    protocol.run(args, steps, out=log)
    assert log.getvalue().splitlines() == ['[run] write', '[done] write']
