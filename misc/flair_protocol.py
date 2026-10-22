#!/usr/bin/env python3
"""Run the FLAIR protocol (docs/source/flair-protocol.rst) from one manifest.

PROTOTYPE: not yet documented or run by the test suite.

The base run is, per sample, flair transcriptome, then flair combine of every
sample's isoforms and flair quantify of every sample's reads on them.  Options add:

  --fusion        flair fusion per sample, its fusion isoforms combined and
                  quantified with the rest
  --spliceevents  flair spliceevents over every sample's reads
  --stat_test conditions
                  flair diffexp between the two conditions of the manifest, on the
                  isoform counts and, with --spliceevents, the event counts
  --stat_test patients
                  diff_iso_usage of each patient's tumors against its normal, on
                  the isoform, event and isoalleles counts.  A patient of tumors
                  only is run through every step and reported, but not tested
  --haplotype     (with --stat_test patients) whatshap phasing of the normal's
                  variants, with its and its tumors' reads, haplotagging, and
                  flair isoalleles of each tumor and its normal; a tumor-only
                  patient's tumors are phased with their own variants and run as
                  one sample

The manifest is tab-separated, with a header naming its columns:

  sample     a unique name, given as each step's sample name
  bam        the reads aligned to the genome, sorted and indexed
  condition  for --stat_test conditions
  patient    for --stat_test patients
  tissue     for --stat_test patients: tumor or normal
  vcf        for --haplotype: the sample's variant calls

Each step writes a marker in OUTDIR/.done when it finishes, and a rerun skips the
steps that have one, so a run that stopped can be continued; --force runs every
step again.  --dry_run prints the commands without running them.  Each step's
output is in OUTDIR/logs/<step>.log.

Examples:

  Isoforms of every sample, combined and quantified, with junctions from the reads
  (intron-prospector), from a manifest of sample and bam:

    sample   bam
    s1       s1.genomealigned.bam
    s2       s2.genomealigned.bam

    flair_protocol.py -m manifest.tsv -g genome.fa -f annotation.gtf -o run

  Two conditions of 3 or more replicates each, with fusions and splice events,
  tested with flair diffexp (a condition column), with short-read junctions:

    sample   bam                    condition
    ctl1     ctl1.genomealigned.bam ctl
    ...
    trt3     trt3.genomealigned.bam trt

    flair_protocol.py -m manifest.tsv -g genome.fa -f annotation.gtf -o run \\
        --junction_star SJ.out.tab --fusion --spliceevents \\
        --stat_test conditions --condition_a ctl --condition_b trt

  Tumors against their patient's normal (patient and tissue columns), with the
  reads phased and haplotagged and isoform usage by haplotype and somatic variant
  (a vcf column).  P2 has a tumor only, so it is run and reported, not tested:

    sample  bam       patient tissue vcf
    P1N     P1N.bam   P1      normal P1N.vcf.gz
    P1T     P1T.bam   P1      tumor  P1T.vcf.gz
    P2T     P2T.bam   P2      tumor  P2T.vcf.gz

    flair_protocol.py -m manifest.tsv -g genome.fa -f annotation.gtf -o run \\
        --stat_test patients --haplotype --spliceevents -t 8

  The commands a run would execute, without running them:

    flair_protocol.py ... --dry_run
"""
import argparse
import collections
import csv
import gzip
import os
import shlex
import shutil
import subprocess
import sys

MANIFEST_COLUMNS = ('sample', 'bam', 'condition', 'patient', 'tissue', 'vcf')
TISSUES = ('tumor', 'normal')
# diffexp's needs: replicates of each condition, a batch name distinct from conditions
MIN_REPLICATES = 3
BATCH = 'b1'
NO_CONDITION = 'none'
PREFIX = 'flair'


class ProtocolError(Exception):
    pass


Step = collections.namedtuple('Step', ('name', 'cmd', 'outputs', 'write'))
Step.__new__.__defaults__ = ((), None)
Step.__doc__ = """One step: cmd, a command to run, or write, a function writing its
outputs (a manifest, a pooled file), and the outputs it must leave"""


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-m', '--manifest', required=True, help='tab-separated manifest, with a header (above)')
    parser.add_argument('-g', '--genome', required=True, help='FASTA of the genome, indexed')
    parser.add_argument('-f', '--gtf', required=True, help='annotation GTF')
    parser.add_argument('-o', '--output', required=True, help='output directory')
    junctions = parser.add_mutually_exclusive_group()
    junctions.add_argument('--junction_bed',
                           help='splice junctions as a BED of supporting reads (default: intron-prospector of '
                                'every BAM, pooled)')
    junctions.add_argument('--junction_star', help='splice junctions as a STAR SJ.out.tab')
    parser.add_argument('-t', '--threads', type=int, default=4, help='threads per step (default: %(default)s)')
    parser.add_argument('--fusion', action='store_true', help='add flair fusion after transcriptome')
    parser.add_argument('--spliceevents', action='store_true', help='add flair spliceevents')
    parser.add_argument('--stat_test', choices=('conditions', 'patients'),
                        help='test between conditions (flair diffexp) or each patient\'s tumors against its normal '
                             '(diff_iso_usage)')
    parser.add_argument('--condition_a', help='with --stat_test conditions, the reference condition')
    parser.add_argument('--condition_b', help='with --stat_test conditions, the condition compared with it')
    parser.add_argument('--haplotype', action='store_true',
                        help='phase and haplotag the reads with whatshap and run flair isoalleles (needs '
                             '--stat_test patients and a vcf for every sample)')
    parser.add_argument('--trust_ends', action='store_true',
                        help='the reads\' ends are reliable: flair quantify --trust_ends rather than --normalize_ends')
    parser.add_argument('--dry_run', action='store_true', help='print the commands without running them')
    parser.add_argument('--force', action='store_true', help='run every step, even those already done')
    parser.add_argument('--flair', default='flair', help='the flair command (default: %(default)s)')
    return parser.parse_args(argv)


def read_manifest(path):
    "the manifest's rows, as dicts of its columns, without blank and '#' lines"
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rt') as fh:
        lines = [line for line in fh if line.strip() and not line.startswith('#')]
    if not lines:
        raise ProtocolError(f'{path} is empty')
    reader = csv.DictReader(lines, delimiter='\t')
    unknown = set(reader.fieldnames) - set(MANIFEST_COLUMNS)
    if unknown:
        raise ProtocolError(f'{path}: unknown columns {", ".join(sorted(unknown))}; the columns are '
                            f'{", ".join(MANIFEST_COLUMNS)}, tab-separated')
    rows = []
    for n, row in enumerate(reader, 2):
        if None in row or any(v is None for v in row.values()):
            raise ProtocolError(f'{path} line {n}: expected {len(reader.fieldnames)} tab-separated fields')
        rows.append({k: v.strip() for k, v in row.items()})
    return rows


def _needed_columns(args):
    needed = ['sample', 'bam']
    if args.stat_test == 'conditions':
        needed.append('condition')
    if args.stat_test == 'patients':
        needed += ['patient', 'tissue']
    if args.haplotype:
        needed.append('vcf')
    return needed


def patients_of(samples):
    """patient -> (normal sample or None, [tumor samples]), in manifest order"""
    patients = collections.OrderedDict()
    for s in samples:
        normal, tumors = patients.setdefault(s['patient'], [None, []])
        if s['tissue'] == 'normal':
            if normal is not None:
                raise ProtocolError(f'patient {s["patient"]} has two normals, {normal["sample"]} and {s["sample"]}')
            patients[s['patient']][0] = s
        else:
            tumors.append(s)
    return {p: (normal, tumors) for p, (normal, tumors) in patients.items()}


def conditions_to_compare(args, samples):
    "the two conditions diffexp compares, condition_a the reference, as diffexp chooses them"
    counts = collections.Counter(s['condition'] for s in samples)
    if args.condition_a or args.condition_b:
        if not (args.condition_a and args.condition_b):
            raise ProtocolError('give both --condition_a and --condition_b')
        pair = (args.condition_a, args.condition_b)
    elif len(counts) == 2:
        pair = tuple(sorted(counts))
    else:
        raise ProtocolError(f'the manifest has {len(counts)} conditions ({", ".join(sorted(counts))}); flair diffexp '
                            'compares two: name them with --condition_a and --condition_b')
    for condition in pair:
        if counts[condition] < MIN_REPLICATES:
            raise ProtocolError(f'condition {condition} has {counts[condition]} samples; flair diffexp needs at least '
                                f'{MIN_REPLICATES} of each condition')
        if condition.isdigit() or condition == BATCH:
            raise ProtocolError(f'condition {condition}: flair diffexp needs condition names that are not numbers '
                                f'and not {BATCH}, the batch name')
    if len(counts) != 2:
        raise ProtocolError('flair diffexp needs a counts matrix of exactly the two conditions it compares: '
                            f'the manifest has {", ".join(sorted(counts))}')
    return pair


def _check_files(args, samples):
    "the inputs exist, and every BAM is indexed"
    for s in samples:
        files = [s['bam']] + ([s['vcf']] if args.haplotype else [])
        for path in files:
            if not os.path.exists(path):
                raise ProtocolError(f'sample {s["sample"]}: {path} does not exist')
        if not any(os.path.exists(s['bam'] + ext) for ext in ('.bai', '.csi')) and \
                not os.path.exists(s['bam'][:-len('.bam')] + '.bai'):
            raise ProtocolError(f'sample {s["sample"]}: {s["bam"]} is not indexed (samtools index)')
    for path in (args.genome, args.gtf, args.junction_bed, args.junction_star):
        if path and not os.path.exists(path):
            raise ProtocolError(f'{path} does not exist')


def _check_tools(args, which):
    "the programs the run's steps call are in PATH, but for a dry run"
    tools = [args.flair]
    if not (args.junction_bed or args.junction_star):
        tools.append('intron-prospector')
    if args.stat_test == 'patients':
        tools.append('diff_iso_usage')
    if args.haplotype:
        tools += ['whatshap', 'samtools']
    missing = [t for t in tools if not which(t)]
    if missing and not args.dry_run:
        raise ProtocolError(f'not in PATH: {", ".join(missing)}')


def check(args, samples, which=shutil.which):
    """Every check of the options and manifest, before any step runs.  Returns the
    patients, of --stat_test patients"""
    needed = _needed_columns(args)
    for column in needed:
        missing = [s.get('sample', '?') for s in samples if not s.get(column)]
        if missing:
            raise ProtocolError(f'the manifest needs a {column} column for this run, filled for every sample: '
                                f'missing for {", ".join(missing)}')
    names = [s['sample'] for s in samples]
    dups = sorted(n for n, c in collections.Counter(names).items() if c > 1)
    if dups:
        raise ProtocolError(f'sample names must be unique: {", ".join(dups)}')
    if args.haplotype and args.stat_test != 'patients':
        raise ProtocolError('--haplotype needs --stat_test patients: haplotypes are only compared within a person')
    _check_files(args, samples)
    patients = None
    if args.stat_test == 'patients':
        bad = sorted({s['tissue'] for s in samples} - set(TISSUES))
        if bad:
            raise ProtocolError(f'tissue must be tumor or normal, not {", ".join(bad)}')
        patients = patients_of(samples)
        for patient, (normal, tumors) in patients.items():
            if not tumors:
                raise ProtocolError(f'patient {patient} has a normal and no tumor')
    if args.stat_test == 'conditions':
        args.condition_a, args.condition_b = conditions_to_compare(args, samples)
    _check_tools(args, which)
    return patients


def _out(args, *parts):
    return os.path.join(args.output, *parts)


def _junction_options(args):
    "the junction options of transcriptome and spliceevents"
    if args.junction_star:
        return ['--junction_star', args.junction_star]
    bed = args.junction_bed or _out(args, 'junctions', 'all.bed')
    return ['--junction_bed', bed, '--junction_support', '2']


def _write_lines(path, lines):
    def write():
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, 'w') as fh:
            fh.writelines(line + '\n' for line in lines)
    return write


def _pool(path, parts):
    def write():
        with open(path, 'w') as out:
            for part in parts:
                with open(part) as fh:
                    shutil.copyfileobj(fh, out)
    return write


def _indexed_vcf(src, dst):
    "a bgzipped, indexed copy of a VCF, dst ending .vcf.gz"
    def write():
        import pysam
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        plain = dst[:-len('.gz')]
        opener = gzip.open if src.endswith('.gz') else open
        with opener(src, 'rb') as fh, open(plain, 'wb') as out:
            shutil.copyfileobj(fh, out)
        pysam.tabix_index(plain, preset='vcf', force=True)
    return write


def _index_vcf(path):
    def write():
        import pysam
        pysam.tabix_index(path, preset='vcf', force=True)
    return write


def junction_steps(args, samples):
    if args.junction_bed or args.junction_star:
        return []
    steps = [Step(f'junctions.{s["sample"]}',
                  ['intron-prospector', f'--genome-fasta={args.genome}',
                   f'--intron-bed6={_out(args, "junctions", s["sample"] + ".bed")}', '-C', '0.0', s['bam']],
                  [_out(args, 'junctions', s['sample'] + '.bed')]) for s in samples]
    parts = [_out(args, 'junctions', s['sample'] + '.bed') for s in samples]
    steps.append(Step('junctions.pool', None, [_out(args, 'junctions', 'all.bed')],
                      _pool(_out(args, 'junctions', 'all.bed'), parts)))
    return steps


def transcriptome_steps(args, samples):
    steps = []
    for s in samples:
        name = s['sample']
        steps.append(Step(f'transcriptome.{name}',
                          [args.flair, 'transcriptome', '-b', s['bam'], '-g', args.genome, '-f', args.gtf,
                           *_junction_options(args), '--sample_name', name,
                           '-o', _out(args, 'transcriptome', name), '-t', str(args.threads)],
                          [_out(args, 'transcriptome', name + '.isoforms.bed')]))
        if args.fusion:
            steps.append(Step(f'fusion.{name}',
                              [args.flair, 'fusion', '-g', args.genome, '-f', args.gtf, '-b', s['bam'],
                               '--sample_name', name, '-o', _out(args, 'fusion', name), '-t', str(args.threads)],
                              [_out(args, 'fusion', name + '.fusions.isoforms.bed')]))
    return steps


def combined_bed(args):
    return _out(args, 'combine', PREFIX + '.combined.isoforms.bed')


def combine_quantify_steps(args, samples):
    beds = [_out(args, 'transcriptome', s['sample'] + '.isoforms.bed') for s in samples]
    if args.fusion:
        beds += [_out(args, 'fusion', s['sample'] + '.fusions.isoforms.bed') for s in samples]
    combine_manifest = _out(args, 'combine', 'manifest.txt')
    quantify_manifest = _out(args, 'quantify', 'manifest.tsv')
    quantify_rows = ['\t'.join([s['sample'], s.get('condition') or NO_CONDITION, BATCH, s['bam']]) for s in samples]
    quantify = [args.flair, 'quantify', '--manifest', quantify_manifest, '-g', args.genome,
                '--isoform_bed', combined_bed(args), '-o', _out(args, 'quantify', PREFIX), '-t', str(args.threads),
                '--trust_ends' if args.trust_ends else '--normalize_ends']
    if args.haplotype:
        quantify.append('--generate_map')
    return [Step('combine.manifest', None, [combine_manifest], _write_lines(combine_manifest, beds)),
            Step('combine', [args.flair, 'combine', '--manifest', combine_manifest, '-g', args.genome, '-f', args.gtf,
                             '-o', _out(args, 'combine', PREFIX)], [combined_bed(args)]),
            Step('quantify.manifest', None, [quantify_manifest], _write_lines(quantify_manifest, quantify_rows)),
            Step('quantify', quantify, [_out(args, 'quantify', PREFIX + '.counts.tsv')])]


def spliceevents_counts(args):
    return _out(args, 'spliceevents', PREFIX + '.spliceevents.counts.tsv')


def spliceevents_steps(args, samples):
    if not args.spliceevents:
        return []
    manifest = _out(args, 'spliceevents', 'manifest.tsv')
    steps = [Step('spliceevents.manifest', None, [manifest],
                  _write_lines(manifest, ['\t'.join([s['sample'], s['bam']]) for s in samples])),
             Step('spliceevents', [args.flair, 'spliceevents', '-m', manifest, '-g', args.genome, '-f', args.gtf,
                                   *_junction_options(args), '-o', _out(args, 'spliceevents', PREFIX),
                                   '-t', str(args.threads)],
                  [spliceevents_counts(args)])]
    if args.stat_test == 'conditions':
        # diffexp reads the conditions from the sample info file beside the counts
        info = _out(args, 'spliceevents', PREFIX + '.spliceevents.sample_info.tsv')
        steps.append(Step('spliceevents.sample_info', None, [info],
                          _write_lines(info, ['sample_id\tcondition\tbatch'] +
                                       ['\t'.join([s['sample'], s['condition'], BATCH]) for s in samples])))
    return steps


def read_map(args, sample):
    return _out(args, 'quantify', f'{PREFIX}.{sample}.read.map.tsv')


def _phase_and_tag(args, phased, vcf, phasing_bams, tagged):
    """steps phasing vcf with the reads of phasing_bams into phased, and haplotagging
    each of tagged, samples, with it"""
    steps = [Step(f'phase.{os.path.basename(phased)}',
                  ['whatshap', 'phase', '--ignore-read-groups', '--reference', args.genome, '-o', phased, vcf,
                   *phasing_bams], [phased]),
             Step(f'phase.index.{os.path.basename(phased)}', None, [phased + '.tbi'], _index_vcf(phased))]
    for s in tagged:
        bam = _out(args, 'phasing', s['sample'] + '.haplotagged.bam')
        steps += [Step(f'haplotag.{s["sample"]}',
                       ['whatshap', 'haplotag', '--ignore-read-groups', '--reference', args.genome, '-o', bam,
                        phased, s['bam']], [bam]),
                  Step(f'haplotag.index.{s["sample"]}', ['samtools', 'index', bam], [bam + '.bai'])]
    return steps


def isoalleles_prefix(args, sample):
    return _out(args, 'isoalleles', sample)


def haplotype_steps(args, patients):
    if not args.haplotype:
        return []
    steps = []
    tagged = lambda s: _out(args, 'phasing', s['sample'] + '.haplotagged.bam')  # noqa: E731
    common = lambda: ['--isoform_bed', combined_bed(args), '-g', args.genome]  # noqa: E731
    for patient, (normal, tumors) in patients.items():
        if normal is not None:
            phased = _out(args, 'phasing', f'{patient}.phased.vcf.gz')
            steps += _phase_and_tag(args, phased, normal['vcf'], [normal['bam']] + [t['bam'] for t in tumors],
                                    [normal] + tumors)
            for t in tumors:
                tumor_vcf = _out(args, 'phasing', f'{t["sample"]}.calls.vcf.gz')
                steps += [Step(f'tumor_vcf.{t["sample"]}', None, [tumor_vcf + '.tbi'],
                               _indexed_vcf(t['vcf'], tumor_vcf)),
                          Step(f'isoalleles.{t["sample"]}',
                               [args.flair, 'isoalleles', *common(), '--bam', tagged(t),
                                '--read_map', read_map(args, t['sample']), '--sample_name', t['sample'],
                                '--normal_bam', tagged(normal), '--normal_read_map', read_map(args, normal['sample']),
                                '--normal_name', normal['sample'], '--phased_vcf', phased, '--tumor_vcf', tumor_vcf,
                                '-o', isoalleles_prefix(args, t['sample']), '-t', str(args.threads)],
                               [isoalleles_prefix(args, t['sample']) + '.isoalleles.counts.tsv'])]
        else:  # tumor-only: each tumor phased with its own variants, one sample
            for t in tumors:
                phased = _out(args, 'phasing', f'{t["sample"]}.phased.vcf.gz')
                steps += _phase_and_tag(args, phased, t['vcf'], [t['bam']], [t])
                steps.append(Step(f'isoalleles.{t["sample"]}',
                                  [args.flair, 'isoalleles', *common(), '--bam', tagged(t),
                                   '--read_map', read_map(args, t['sample']), '--phased_vcf', phased,
                                   '--sample_name', t['sample'], '-o', isoalleles_prefix(args, t['sample']),
                                   '-t', str(args.threads)],
                                  [isoalleles_prefix(args, t['sample']) + '.isoalleles.counts.tsv']))
    return steps


def _diffexp_steps(args):
    "flair diffexp of the isoform counts, and of the event counts with --spliceevents"
    tables = [('isoforms', _out(args, 'quantify', PREFIX + '.counts.tsv'))]
    if args.spliceevents:
        tables.append(('spliceevents', spliceevents_counts(args)))
    steps = []
    for label, counts in tables:
        out_dir = _out(args, 'diffexp', label)
        steps.append(Step(f'diffexp.{label}',
                          [args.flair, 'diffexp', '--counts_matrix', counts, '--condition_a', args.condition_a,
                           '--condition_b', args.condition_b, '-o', out_dir, '-t', str(args.threads),
                           '--overwrite_output'], [out_dir]))
    return steps


def _usage_steps(args, patients):
    """diff_iso_usage of each tumor against its patient's normal, on every counts
    table of the run; a tumor-only patient is reported, not tested"""
    steps = []
    for patient, (normal, tumors) in patients.items():
        if normal is None:
            continue
        for t in tumors:
            tables = [('isoforms', _out(args, 'quantify', PREFIX + '.counts.tsv'))]
            if args.spliceevents:
                tables.append(('spliceevents', spliceevents_counts(args)))
            if args.haplotype:
                tables.append(('isoalleles', isoalleles_prefix(args, t['sample']) + '.isoalleles.counts.tsv'))
            for label, counts in tables:
                usage = _out(args, 'stats', f'{patient}.{t["sample"]}.{label}.usage.tsv')
                steps.append(Step(f'usage.{t["sample"]}.{label}',
                                  ['diff_iso_usage', counts, t['sample'], normal['sample'], usage], [usage]))
    return steps


def stat_steps(args, patients):
    if args.stat_test == 'conditions':
        return _diffexp_steps(args)
    if args.stat_test == 'patients':
        return _usage_steps(args, patients)
    return []


def build_steps(args, samples, patients=None):
    "every step of the run, in order"
    return (junction_steps(args, samples) + transcriptome_steps(args, samples) +
            combine_quantify_steps(args, samples) + spliceevents_steps(args, samples) +
            haplotype_steps(args, patients) + stat_steps(args, patients))


def run(args, steps, out=sys.stdout):
    "run the steps not yet done, in order"
    done_dir, log_dir = _out(args, '.done'), _out(args, 'logs')
    if not args.dry_run:
        os.makedirs(done_dir, exist_ok=True)
        os.makedirs(log_dir, exist_ok=True)
    for step in steps:
        marker = os.path.join(done_dir, step.name)
        if os.path.exists(marker) and not args.force:
            print(f'[done] {step.name}', file=out)
            continue
        if args.dry_run:
            print(shlex.join(step.cmd) if step.cmd else f'# write {", ".join(step.outputs)}', file=out)
            continue
        print(f'[run] {step.name}', file=out, flush=True)
        for path in step.outputs:
            os.makedirs(os.path.dirname(path) or '.', exist_ok=True)
        if step.cmd:
            log = os.path.join(log_dir, step.name + '.log')
            with open(log, 'w') as fh:
                result = subprocess.run(step.cmd, stdout=fh, stderr=subprocess.STDOUT)
            if result.returncode != 0:
                raise ProtocolError(f'step {step.name} failed (exit {result.returncode}); see {log}')
        else:
            step.write()
        missing = [p for p in step.outputs if not os.path.exists(p)]
        if missing:
            raise ProtocolError(f'step {step.name} finished without writing {", ".join(missing)}')
        open(marker, 'w').close()


def summary(args, patients):
    lines = [f'combined isoforms: {combined_bed(args)}',
             f'isoform counts: {_out(args, "quantify", PREFIX + ".counts.tsv")}']
    if args.spliceevents:
        lines.append(f'event counts: {spliceevents_counts(args)}')
    if args.stat_test == 'conditions':
        lines.append(f'diffexp: {_out(args, "diffexp")}')
    if args.stat_test == 'patients':
        lines.append(f'tumor/normal usage tests: {_out(args, "stats")}')
        untested = [p for p, (normal, _) in patients.items() if normal is None]
        if untested:
            lines.append(f'tumor-only patients, reported but not tested: {", ".join(untested)}')
    if args.haplotype:
        lines.append(f'isoalleles: {_out(args, "isoalleles")}')
    return lines


def main(argv=None):
    args = parse_args(argv)
    try:
        samples = read_manifest(args.manifest)
        patients = check(args, samples)
        run(args, build_steps(args, samples, patients))
    except ProtocolError as ex:
        print(f'flair_protocol: {ex}', file=sys.stderr)
        return 1
    if not args.dry_run:
        print('\n'.join(summary(args, patients)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
