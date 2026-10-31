"""
flair isoalleles' core: variants from a phased VCF and a tumor's calls, haplotagged
reads grouped by haplotype and somatic set, and each group's sequence and protein.

The fixture is one isoform on chr1, + strand, of exons [100, 200), [400, 500) and
[700, 900), whose transcript is a 50 base UTR, ATG, 80 GCT codons and TAA.  Codon 20
(A21) has its middle base at chr1:412 (1-based), a germline het C>A (p.A21D) phased
on H2, and codon 30 (A31) at chr1:442, a somatic C>T (p.A31V).
"""
import gzip
import subprocess
import sys
from types import SimpleNamespace
import pysam
import pytest
from flair import FlairInputDataError
from flair.flair_bed import FlairBed
from flair.flair_isoalleles import check_options, isoalleles
from flair.isoallele_genes import GUARANTEED, TERMINAL, GeneOutput, process_gene, site_class
from flair.isoallele_sequences import build_transcript, call_cds, frameshift_residue, protein_change
from flair.isoallele_variants import (HOMOZYGOUS, PHASED_HET, UNPHASED_HET, ReadAlleles, Variant, germline_variants,
                                      parse_variant)
from flair.read_map_index import ReadMapIndex

EXONS = [(100, 200), (400, 500), (700, 900)]
TRANSCRIPT = 'C' * 50 + 'ATG' + 'GCT' * 80 + 'TAA' + 'C' * 104
GERMLINE = Variant('chr1', 412, 'C', 'A')
SOMATIC = Variant('chr1', 442, 'C', 'T')
CIGAR = '100M200N100M200N200M'


def _genome_seq():
    seq = ['G'] * 1000
    for (start, end), offset in zip(EXONS, (0, 100, 200)):
        seq[start:end] = TRANSCRIPT[offset:offset + end - start]
    return ''.join(seq)


@pytest.fixture
def genome(tmp_path):
    path = str(tmp_path / 'genome.fa')
    with open(path, 'w') as fh:
        fh.write('>chr1\n' + _genome_seq() + '\n')
    pysam.faidx(path)
    with pysam.FastaFile(path) as fasta:
        yield fasta


def _isoform():
    blocks = ','.join(str(e - s) for s, e in EXONS), ','.join(str(s - 100) for s, _ in EXONS)
    row = ['chr1', '100', '900', 'iso1', '0', '+', '150', '796', '0', '3', *blocks,
           'G1', '', 'G1', '', '', 'PRO', 'basic', '', '', '', '', '', 'FLP00000001']
    return FlairBed.parse(row)


def _transcript_index(gpos0):
    "a 0-based genome position's index in the transcript"
    offset = 0
    for s, e in EXONS:
        if s <= gpos0 < e:
            return offset + gpos0 - s
        offset += e - s
    raise ValueError(gpos0)


def _read(header, name, alleles=(), hp=None, delete=None):
    """a full-length read of the isoform carrying the alleles given, with HP tag hp,
    and with the base at genome position delete (0-based) deleted"""
    seq = list(TRANSCRIPT)
    for key in alleles:
        seq[_transcript_index(key.pos - 1)] = key.alt
    read = pysam.AlignedSegment(header)
    read.query_name, read.reference_id, read.reference_start, read.mapping_quality = name, 0, 100, 60
    if delete is None:
        read.cigarstring = CIGAR
    else:
        i = delete - 400
        read.cigarstring = f'100M200N{i}M1D{99 - i}M200N200M'
        del seq[_transcript_index(delete)]
    read.query_sequence = ''.join(seq)
    if hp is not None:
        read.set_tag('HP', hp)
        read.set_tag('PS', 100)
    return read


def _bam(tmp_path, name, reads_spec):
    "a sorted, indexed bam of reads_spec, (count, alleles, hp, delete), and its read map"
    header = pysam.AlignmentHeader.from_dict({'SQ': [{'SN': 'chr1', 'LN': 1000}]})
    unsorted, path = str(tmp_path / f'{name}.unsorted.bam'), str(tmp_path / f'{name}.bam')
    names = []
    with pysam.AlignmentFile(unsorted, 'wb', header=header) as bam:
        for n, (count, alleles, hp, delete) in enumerate(reads_spec):
            for k in range(count):
                names.append(f'{name}_{n}_{k}')
                bam.write(_read(header, names[-1], alleles, hp, delete))
    pysam.sort('-o', path, unsorted)
    pysam.index(path)
    read_map = str(tmp_path / f'{name}.read.map.tsv')
    with open(read_map, 'w') as fh:
        fh.write('name\treads\niso1\t' + ','.join(names) + '\n')
    return path, read_map


def _vcf(tmp_path, name, records):
    "a bgzipped, indexed single-sample VCF of records (pos, ref, alt, GT, PS or None)"
    path = str(tmp_path / f'{name}.vcf')
    with open(path, 'w') as fh:
        fh.write('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n'
                 '##FORMAT=<ID=GT,Number=1,Type=String,Description="genotype">\n'
                 '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="phase set">\n'
                 '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n')
        for pos, ref, alt, gt, ps in records:
            fmt, val = ('GT:PS', f'{gt}:{ps}') if ps is not None else ('GT', gt)
            fh.write(f'chr1\t{pos}\t.\t{ref}\t{alt}\t50\tPASS\t.\t{fmt}\t{val}\n')
    return pysam.tabix_index(path, preset='vcf', force=True)


def _ctx(tmp_path, genome, samples, phased_vcf, *, normal=None, somatic_mode=None, tumor_vcf=None,
         min_normal_depth=10, keep_labels=False, germline_vcf=None):
    "a worker context of samples, {name: (bam, read map)}, as flair_isoalleles makes one"
    order = [normal] + [s for s in samples if s != normal] if normal else list(samples)
    return SimpleNamespace(
        genome=genome, samples=order, output_samples=list(samples), normal=normal,
        bams={s: pysam.AlignmentFile(bam) for s, (bam, _) in samples.items()},
        map_readers={s: ReadMapIndex(rmap, str(tmp_path)).open_reader() for s, (_, rmap) in samples.items()},
        phased_vcf=pysam.VariantFile(phased_vcf), phased_vcf_path=phased_vcf,
        germline_vcf=pysam.VariantFile(germline_vcf) if germline_vcf else None, germline_vcf_path=germline_vcf,
        somatic_mode=somatic_mode, tumor_vcf=pysam.VariantFile(tumor_vcf) if tumor_vcf else None,
        tumor_vcf_path=tumor_vcf, somatic_vcf=None, somatic_vcf_path=None, variant_keys=[],
        min_normal_depth=min_normal_depth, carrier_min=0.9, min_tagged=10, keep_labels=keep_labels)


def _run(ctx):
    out = GeneOutput()
    process_gene(ctx, 'G1', [_isoform()], out)
    return out


def _sequences(out):
    "group -> (productivity, protein change, somatic protein change)"
    return {row[2]: (row[5], row[8], row[9]) for row in out.sequences}


# sequences and proteins

def test_a_missense_variant_changes_the_protein(genome):
    iso = _isoform()
    ref = call_cds(iso, *build_transcript(genome, iso, [])[:3], False)
    alt = call_cds(iso, *build_transcript(genome, iso, [GERMLINE])[:3], False)
    assert ref[0] == alt[0] == 'PRO' and ref[2] == alt[2] == (150, 796)
    assert protein_change(ref[1], alt[1]) == 'p.A21D'


def test_a_variant_whose_ref_disagrees_with_the_genome_is_not_applied(genome):
    iso = _isoform()
    assert build_transcript(genome, iso, [Variant('chr1', 412, 'G', 'A')])[3] == [Variant('chr1', 412, 'G', 'A')]


def test_protein_changes_are_named_as_hgvs_does():
    assert protein_change('MAAK', 'MAK') == 'p.A3del'  # shifted to the C terminus
    assert protein_change('MAKL', 'MAKKL') == 'p.K3dup'
    assert protein_change('MAKL', 'MA') == 'p.K3*'
    assert protein_change('MAKL', 'MAKLQ') == 'p.*5ext'
    assert protein_change('MAKL', 'MAQR', 2) == 'p.K3fs'
    # a missense before the frameshift is named on its own
    assert protein_change('MAKLW', 'MRKQS', 3) == 'p.A2R;p.L4fs'


def test_the_frameshift_is_at_the_residue_of_the_first_indel_shifting_the_frame(genome):
    iso = _isoform()
    # one base deleted in codon 30 (the residue of index 30), then two inserted after it
    deletion, insertion = Variant('chr1', 441, 'GC', 'G'), Variant('chr1', 470, 'C', 'CAA')
    assert frameshift_residue(iso, [deletion]) == 30
    assert frameshift_residue(iso, [deletion, insertion]) == 30
    assert frameshift_residue(iso, [Variant('chr1', 441, 'GCTG', 'G')]) is None


# variants and reads

def test_germline_variants_by_genotype(tmp_path):
    path = _vcf(tmp_path, 'g', [(300, 'A', 'G', '1/1', None), (412, 'C', 'A', '0|1', 100), (600, 'T', 'C', '0/1', None)])
    with pysam.VariantFile(path) as vcf:
        kinds = [(g.key.pos, g.kind, g.haplotype) for g in germline_variants(vcf.fetch('chr1', 0, 1000))]
    assert kinds == [(300, HOMOZYGOUS, None), (412, PHASED_HET, 'H2'), (600, UNPHASED_HET, None)]


def test_germline_calls_phased_by_a_vcf_of_only_their_hets(tmp_path):
    calls = _vcf(tmp_path, 'calls', [(300, 'A', 'G', '1/1', None), (412, 'C', 'A', '0/1', None),
                                     (600, 'T', 'C', '0/1', None)])
    phased = _vcf(tmp_path, 'phased', [(412, 'C', 'A', '1|0', 100)])
    with pysam.VariantFile(calls) as calls_vcf, pysam.VariantFile(phased) as phased_vcf:
        kinds = [(g.key.pos, g.kind, g.haplotype, g.phase_set)
                 for g in germline_variants(calls_vcf.fetch('chr1', 0, 1000), list(phased_vcf.fetch('chr1', 0, 1000)))]
    assert kinds == [(300, HOMOZYGOUS, None, None), (412, PHASED_HET, 'H1', 100), (600, UNPHASED_HET, None, None)]


def test_a_variant_is_given_with_or_without_its_alleles():
    assert parse_variant('chr7:140753336:A>T') == parse_variant('chr7:140753336A>T') == ('chr7', 140753336, 'A', 'T')
    assert parse_variant('chr7:140753336') == ('chr7', 140753336, None, None)
    with pytest.raises(FlairInputDataError):
        parse_variant('chr7')


def test_a_reads_allele_at_a_site():
    header = pysam.AlignmentHeader.from_dict({'SQ': [{'SN': 'chr1', 'LN': 1000}]})
    carrier, deleted = ReadAlleles(_read(header, 'a', [SOMATIC])), ReadAlleles(_read(header, 'b', delete=441))
    assert (carrier.allele(SOMATIC), carrier.allele(GERMLINE)) == ('alt', 'ref')
    assert deleted.allele(SOMATIC) is None


def test_site_classes():
    iso = _isoform()
    assert site_class(iso, SOMATIC) == GUARANTEED  # an internal exon
    assert site_class(iso, Variant('chr1', 197, 'C', 'T')) == GUARANTEED  # next to the first splice site
    assert site_class(iso, Variant('chr1', 120, 'C', 'T')) == TERMINAL  # the 5' UTR
    assert site_class(iso, Variant('chr1', 300, 'G', 'T')) is None  # an intron


# one sample: haplotypes

def test_one_samples_reads_by_haplotype_with_each_haplotypes_protein(tmp_path, genome):
    s1 = _bam(tmp_path, 's1', [(12, [], 1, None), (12, [GERMLINE], 2, None), (3, [], None, None)])
    phased = _vcf(tmp_path, 'phased', [(412, 'C', 'A', '0|1', 100)])
    out = _run(_ctx(tmp_path, genome, {'s1': s1}, phased))
    assert out.counts == [['iso1', 'iso1|H1', [12]], ['iso1', 'iso1|H2', [12]], ['iso1', 'iso1|untagged', [3]]]
    assert _sequences(out) == {'H1': ('PRO', 'none', 'none'), 'H2': ('PRO', 'p.A21D', 'none')}
    # the unchanged protein keeps the isoform's id
    assert {row[2]: row[7] for row in out.sequences}['H1'] == 'FLP00000001'


def test_an_unphased_het_is_placed_by_its_reads(tmp_path, genome):
    s1 = _bam(tmp_path, 's1', [(12, [], 1, None), (12, [GERMLINE], 2, None)])
    unphased = _vcf(tmp_path, 'unphased', [(412, 'C', 'A', '0/1', None)])
    out = _run(_ctx(tmp_path, genome, {'s1': s1}, unphased))
    assert _sequences(out)['H2'][1] == 'p.A21D'
    assert ['G1', 'iso1', 'chr1:412:C>A', 'germline unphased het, placed by reads', 'H2', ''] in out.variants


# a pair: somatic sets

def _pair(tmp_path, genome, **kwargs):
    tumor = _bam(tmp_path, 'tumor', [(10, [], 1, None), (10, [GERMLINE, SOMATIC], 2, None),
                                     (3, [GERMLINE], 2, None), (2, [GERMLINE], 2, 441)])
    normal = _bam(tmp_path, 'normal', [(12, [], 1, None), (12, [GERMLINE], 2, None)])
    phased = _vcf(tmp_path, 'phased', [(412, 'C', 'A', '0|1', 100)])
    tumor_vcf = _vcf(tmp_path, 'tumor', [(412, 'C', 'A', '0/1', None), (442, 'C', 'T', '0/1', None)])
    return _run(_ctx(tmp_path, genome, {'tumor': tumor, 'normal': normal}, phased, normal='normal',
                     somatic_mode='tumor_vcf', tumor_vcf=tumor_vcf, **kwargs))


def test_a_tumors_somatic_variant_makes_a_group_of_its_reads(tmp_path, genome):
    out = _pair(tmp_path, genome)
    assert out.somatic_sets == [['G1', 'S1', 'chr1:442:C>T', 10]]
    assert out.counts == [['iso1', 'iso1|H1', [10, 12]], ['iso1', 'iso1|H2', [3, 12]],
                          ['iso1', 'iso1|H2-S1', [10, 0]], ['iso1', 'iso1|H2-S?', [2, 0]]]
    assert _sequences(out) == {'H1': ('PRO', 'none', 'none'), 'H2': ('PRO', 'p.A21D', 'none'),
                               'H2-S1': ('PRO', 'p.A21D;p.A31V', 'p.A31V')}
    assert ['G1', 'chr1', 'tumor', 0] not in out.genes
    assert [row[:2] for row in out.genes] == [['G1', 'tumor'], ['G1', 'normal']]


def test_a_tumor_call_the_normal_makes_but_did_not_phase_is_germline(tmp_path, genome):
    # the tumor's call at 442 is also the normal's, a het it didn't phase: germline, not somatic
    calls = _vcf(tmp_path, 'calls', [(412, 'C', 'A', '0/1', None), (442, 'C', 'T', '0/1', None)])
    out = _pair(tmp_path, genome, germline_vcf=calls)
    assert out.somatic_sets == [] and out.excluded == []


def test_a_somatic_callers_calls_are_somatic_whatever_the_normals_reads(tmp_path, genome):
    tumor = _bam(tmp_path, 'tumor', [(10, [], 1, None), (10, [GERMLINE, SOMATIC], 2, None)])
    normal = _bam(tmp_path, 'normal', [(12, [], 1, None), (12, [GERMLINE], 2, None)])
    phased = _vcf(tmp_path, 'phased', [(412, 'C', 'A', '0|1', 100)])
    somatic = _vcf(tmp_path, 'somatic', [(442, 'C', 'T', '0/1', None)])
    ctx = _ctx(tmp_path, genome, {'tumor': tumor, 'normal': normal}, phased, normal='normal',
               somatic_mode='somatic_vcf', min_normal_depth=30)
    ctx.somatic_vcf, ctx.somatic_vcf_path = pysam.VariantFile(somatic), somatic
    out = _run(ctx)
    # not excluded for its 24 normal reads, but warned of
    assert out.somatic_sets == [['G1', 'S1', 'chr1:442:C>T', 10]] and out.excluded == []
    assert any(row[2] == 'LOW NORMAL COVERAGE AT SOMATIC SITE' for row in out.warnings)


def test_a_tumor_only_call_of_too_few_normal_reads_is_excluded(tmp_path, genome):
    out = _pair(tmp_path, genome, min_normal_depth=30)
    assert out.somatic_sets == [] and [row[:3] for row in out.excluded] == [['G1', 'chr1:442:C>T', 24]]
    assert {row[1] for row in out.counts} == {'iso1|H1', 'iso1|H2'}


def test_reads_are_labelled_by_isoform_and_group(tmp_path, genome):
    out = _pair(tmp_path, genome, keep_labels=True)
    labels = dict(out.labels['tumor'])
    assert labels['tumor_1_0'] == 'iso1|H2-S1' and labels['tumor_3_0'] == 'iso1|H2-S?'


# the whole command

def _bed_file(tmp_path):
    path = str(tmp_path / 'isoforms.bed')
    with open(path, 'w') as fh:
        fh.write('\t'.join(map(str, _isoform().toRow())) + '\n')
    return path


def _table(path):
    with open(path) as fh:
        return [line.rstrip('\n').split('\t') for line in fh]


@pytest.mark.parametrize('threads', [1, 2])
def test_isoalleles_of_a_pair(tmp_path, genome, threads):
    tumor = _bam(tmp_path, 'tumor', [(10, [], 1, None), (10, [GERMLINE, SOMATIC], 2, None)])
    normal = _bam(tmp_path, 'normal', [(12, [], 1, None), (12, [GERMLINE], 2, None)])
    phased = _vcf(tmp_path, 'phased', [(412, 'C', 'A', '0|1', 100)])
    tumor_vcf = _vcf(tmp_path, 'tumor', [(412, 'C', 'A', '0/1', None), (442, 'C', 'T', '0/1', None)])
    out = str(tmp_path / 'out')
    isoalleles(isoform_bed=_bed_file(tmp_path), genome=genome.filename.decode(), bam=tumor[0], read_map=tumor[1],
               phased_vcf=phased, normal_bam=normal[0], normal_read_map=normal[1], tumor_vcf=tumor_vcf,
               output=out, threads=threads, generate_map=True, label_bams=True, temp_dir=str(tmp_path))
    assert _table(out + '.isoalleles.counts.tsv') == [['gene_id', 'isoform_id', 'tumor', 'normal'],
                                                      ['iso1', 'iso1|H1', '10', '12'],
                                                      ['iso1', 'iso1|H2', '0', '12'],
                                                      ['iso1', 'iso1|H2-S1', '10', '0']]
    sequences = {row[2]: row[8] for row in _table(out + '.isoalleles.sequences.tsv')[1:]}
    assert sequences == {'H1': 'none', 'H2': 'p.A21D', 'H2-S1': 'p.A21D;p.A31V'}
    assert [row[3] for row in _table(out + '.isoalleles.bed')] == ['iso1|H1', 'iso1|H2', 'iso1|H2-S1']
    with pysam.AlignmentFile(out + '.isoalleles.tumor.bam') as bam:
        assert sorted({read.get_tag('IA') for read in bam}) == ['iso1|H1', 'iso1|H2-S1']
    assert _table(out + '.isoalleles.normal.read.map.tsv')[1][0] == 'iso1|H1'


def test_the_flair_command_runs_one_sample(tmp_path, genome):
    # through flair's own command line, with its global options, as a user runs it
    s1 = _bam(tmp_path, 's1', [(12, [], 1, None), (12, [GERMLINE], 2, None)])
    phased = _vcf(tmp_path, 'phased', [(412, 'C', 'A', '0|1', 100)])
    out = str(tmp_path / 'out')
    subprocess.run([sys.executable, '-m', 'flair.flair_cli', '--compress', 'isoalleles',
                    '--isoform_bed', _bed_file(tmp_path), '-g', genome.filename.decode(), '--bam', s1[0],
                    '--read_map', s1[1], '--phased_vcf', phased, '--sample_name', 's1', '-t', '1', '-o', out,
                    '--temp_dir', str(tmp_path)], check=True)
    with gzip.open(out + '.isoalleles.counts.tsv.gz', 'rt') as fh:
        assert fh.read().splitlines() == ['gene_id\tisoform_id\ts1', 'iso1\tiso1|H1\t12', 'iso1\tiso1|H2\t12']


def test_a_pair_needs_its_somatic_variants(tmp_path):
    with pytest.raises(FlairInputDataError, match='needs the somatic variants'):
        check_options(bam='t.bam', sample_name=None, normal_bam='n.bam', normal_read_map='n.map', normal_name=None,
                      somatic={'tumor_vcf': None, 'somatic_vcf': None, 'variant': None}, paths=[])
    with pytest.raises(FlairInputDataError, match='are of a tumor and normal pair'):
        check_options(bam='t.bam', sample_name=None, normal_bam=None, normal_read_map=None, normal_name=None,
                      somatic={'tumor_vcf': 't.vcf.gz', 'somatic_vcf': None, 'variant': None}, paths=[])
