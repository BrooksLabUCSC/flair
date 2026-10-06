import io
from types import SimpleNamespace
from flair.isoform_data import Isoform
from flair.flair_transcriptome import assign_single_exon_reads


def _isoform(name, start, end, strand='+'):
    iso = Isoform('chr1', strand, (), start, end)
    iso.name = name
    return iso


def _read(name, start, end, strand='+'):
    return SimpleNamespace(name=name, chrom='chr1', start=start, end=end, strand=strand)


def _assign(reads, isoforms, **opts):
    ends_fh = io.StringIO()
    assign_single_exon_reads(reads, {iso.name: iso for iso in isoforms}, ends_fh,
                             trust_ends=opts.get('trust_ends', False), trust_strand=opts.get('trust_strand', False))
    assigned = {line.split('\t')[0]: line.split('\t')[1] for line in ends_fh.getvalue().splitlines()}
    return assigned


def test_read_goes_to_the_overlapping_isoform_with_the_closest_ends():
    assigned = _assign([_read('r1', 1000, 2000), _read('r2', 1010, 2950)],
                       [_isoform('short', 1000, 2000), _isoform('long', 1000, 3000)])
    assert assigned == {'r1': 'short', 'r2': 'long'}


def test_read_must_overlap_half_of_the_isoform_and_half_of_itself():
    assigned = _assign([_read('in', 1100, 1900), _read('short', 1000, 1400), _read('mostly_outside', 1500, 4000)],
                       [_isoform('iso', 1000, 2000)])
    assert assigned == {'in': 'iso'}


def test_ends_are_written_from_the_transcript_5_prime_end():
    ends_fh = io.StringIO()
    assign_single_exon_reads([_read('r', 1100, 1950, '-')], {'iso': _isoform('iso', 1000, 2000, '-')}, ends_fh,
                             trust_ends=False, trust_strand=False)
    assert ends_fh.getvalue().split('\t')[4] == '50' and ends_fh.getvalue().rstrip().split('\t')[7] == '100'


def test_trust_strand_and_trust_ends():
    isoforms = [_isoform('iso', 1000, 2000, '+')]
    assert _assign([_read('r', 1000, 2000, '-')], isoforms, trust_strand=True) == {}
    assert _assign([_read('near', 1040, 1960), _read('far', 1100, 2000)], isoforms, trust_ends=True) == {'near': 'iso'}
