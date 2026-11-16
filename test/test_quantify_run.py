"""flair quantify run end to end on a small genome, with overlapping genes and reads
made from the genome: each read counted once, for the isoform it is full-length for."""
import random
import pysam
import pytest
from flair import FlairInputDataError
from flair.flair_bed import FlairBed
from flair.pycbio.hgdata.bed import BedBlock
from flair.flair_quantify import quantify

_RNG = random.Random(1)
GENOME = ''.join(_RNG.choice('ACGT') for _ in range(6000))
_HEADER = pysam.AlignmentHeader.from_dict({'SQ': [{'SN': 'chr1', 'LN': len(GENOME)}]})

# G1 and G2 overlap, their isoforms sharing the junction 1700-2000; G3 and G4 overlap,
# single-exon isoforms one inside the other
ISOFORMS = [('A1', 'G1', [(1000, 1200), (1500, 1700), (2000, 2200)]),
            ('B1', 'G2', [(1500, 1700), (2000, 2200), (2600, 2800)]),
            ('S1', 'G3', [(3000, 3500)]),
            ('S2', 'G4', [(3100, 3400)])]


def _cigar(exons):
    "the cigar of a read aligned over exons"
    cigar = []
    for i, (start, end) in enumerate(exons):
        if i > 0:
            cigar.append((pysam.CREF_SKIP, start - exons[i - 1][1]))
        cigar.append((pysam.CMATCH, end - start))
    return cigar


# name, exons covered, genome positions mutated
READS = [('a_full', [(1000, 1200), (1500, 1700), (2000, 2200)], ()),
         ('b_full', [(1500, 1700), (2000, 2200), (2600, 2800)], ()),
         # mismatches next to A1's first splice site: not cleanly aligned there, so
         # realigned, and counted, a splice site window allowing 2 mistakes
         ('a_realigned', [(1000, 1200), (1500, 1700), (2000, 2200)], (1196, 1197)),
         # its start 150 bases into A1's first exon: counted, but not with --trust_ends
         ('a_short_start', [(1150, 1200), (1500, 1700), (2000, 2200)], ()),
         # the shared junction only: full-length for neither A1 nor B1
         ('shared_only', [(1550, 1700), (2000, 2150)], ()),
         ('s1_read', [(3010, 3490)], ()),
         ('s2_read', [(3110, 3390)], ()),
         # less than half of either single-exon isoform
         ('too_short', [(3200, 3300)], ())]


def _inputs(tmp_path):
    genome = tmp_path / 'genome.fa'
    genome.write_text('>chr1\n' + GENOME + '\n')
    pysam.faidx(str(genome))
    bed = tmp_path / 'isoforms.bed'
    with open(bed, 'w') as fh:
        for name, gene_id, exons in ISOFORMS:
            FlairBed('chr1', exons[0][0], exons[-1][1], name, strand='+', blocks=[BedBlock(s, e) for s, e in exons],
                     gene_id=gene_id).write(fh)
    bam = str(tmp_path / 'reads.bam')
    with pysam.AlignmentFile(bam, 'wb', header=_HEADER) as out:
        for name, exons, mutated in sorted(READS, key=lambda read: read[1][0][0]):
            seg = pysam.AlignedSegment(_HEADER)
            seg.query_name, seg.reference_id, seg.reference_start = name, 0, exons[0][0]
            seg.cigartuples = _cigar(exons)
            seq = ''.join(GENOME[s:e] for s, e in exons)
            offsets = {}
            pos = 0
            for s, e in exons:
                for g in range(s, e):
                    offsets[g] = pos
                    pos += 1
            seq = list(seq)
            for g in mutated:
                seq[offsets[g]] = 'A' if seq[offsets[g]] != 'A' else 'C'
            seg.query_sequence = ''.join(seq)
            out.write(seg)
    pysam.index(bam)
    manifest = tmp_path / 'manifest.tsv'
    manifest.write_text(f's1\tA\tb0\t{bam}\n')
    return str(manifest), str(genome), str(bed)


def _run(tmp_path, threads, trust_ends=False):
    manifest, genome, bed = _inputs(tmp_path)
    output = str(tmp_path / 'out')
    quantify(manifest=manifest, genome=genome, isoform_bed=bed, output=output, threads=threads,
             trust_ends=trust_ends, generate_map=True, normalize_ends=False,
             temp_dir=str(tmp_path / 'tmp'))
    counts = [line.rstrip('\n').split('\t') for line in open(output + '.isoforms.counts.tsv')]
    read_map = {}
    lines = open(output + '.isoforms.s1.read.map.tsv').read().splitlines()
    assert lines[0] == 'name\treads'
    for line in lines[1:]:
        name, reads = line.rstrip('\n').split('\t')
        read_map.setdefault(name, []).extend(reads.split(','))
    return counts, {name: sorted(reads) for name, reads in read_map.items()}


@pytest.mark.parametrize('threads', [1, 2])
def test_each_read_counted_once_for_its_full_length_isoform(tmp_path, threads):
    counts, read_map = _run(tmp_path, threads)
    assert counts == [['gene_id', 'isoform_id', 's1'], ['G1', 'A1', '3'], ['G2', 'B1', '1'], ['G3', 'S1', '1'],
                      ['G4', 'S2', '1']]
    assert read_map == {'A1': ['a_full', 'a_realigned', 'a_short_start'], 'B1': ['b_full'],
                        'S1': ['s1_read'], 'S2': ['s2_read']}


def test_trust_ends_needs_ends_near_the_isoforms(tmp_path):
    _, read_map = _run(tmp_path, 1, trust_ends=True)
    assert read_map['A1'] == ['a_full', 'a_realigned']


def _subset_bam(src_bam, out_bam, read_names):
    "a BAM of just the named reads of another, standing in for a second sample"
    with pysam.AlignmentFile(src_bam, 'rb') as inp, \
         pysam.AlignmentFile(out_bam, 'wb', template=inp) as out:
        for seg in [seg for seg in inp if seg.query_name in read_names]:
            out.write(seg)
    pysam.index(out_bam)


def _read_map(output, sample):
    lines = open(f'{output}.isoforms.{sample}.read.map.tsv').read().splitlines()
    assert lines[0] == 'name\treads'
    return {name: sorted(reads.split(','))
            for name, reads in (line.split('\t') for line in lines[1:])}


def test_read_map_per_sample_matches_its_counts_column(tmp_path):
    """Each sample gets its own read map, holding that sample's reads and no other's,
    and the counts column is the number of reads the map names."""
    manifest, genome, bed = _inputs(tmp_path)
    s1_bam = str(tmp_path / 'reads.bam')
    s2_bam = str(tmp_path / 'reads.s2.bam')
    _subset_bam(s1_bam, s2_bam, {'a_full', 's1_read'})
    with open(manifest, 'w') as fh:
        fh.write(f's1\tA\tb0\t{s1_bam}\ns2\tB\tb0\t{s2_bam}\n')

    output = str(tmp_path / 'out')
    quantify(manifest=manifest, genome=genome, isoform_bed=bed, output=output, threads=1,
             trust_ends=False, generate_map=True, normalize_ends=False, temp_dir=str(tmp_path / 'tmp'))

    counts = [line.split('\t') for line in open(output + '.isoforms.counts.tsv').read().splitlines()]
    assert counts[0] == ['gene_id', 'isoform_id', 's1', 's2']
    maps = {'s1': _read_map(output, 's1'), 's2': _read_map(output, 's2')}
    assert maps['s1'] == {'A1': ['a_full', 'a_realigned', 'a_short_start'], 'B1': ['b_full'],
                          'S1': ['s1_read'], 'S2': ['s2_read']}
    assert maps['s2'] == {'A1': ['a_full'], 'S1': ['s1_read']}
    for _, isoform_id, s1_count, s2_count in counts[1:]:
        assert int(s1_count) == len(maps['s1'].get(isoform_id, []))
        assert int(s2_count) == len(maps['s2'].get(isoform_id, []))


def _bam_with_repeated_id(src_bam, out_bam, copied_read, repeated_id):
    "a BAM with one read's alignment also present under another read's id"
    unsorted = out_bam + '.unsorted.bam'
    with pysam.AlignmentFile(src_bam, 'rb') as inp, \
         pysam.AlignmentFile(unsorted, 'wb', template=inp) as out:
        for seg in inp:
            out.write(seg)
            if seg.query_name == copied_read:
                repeat = pysam.AlignedSegment.fromstring(seg.to_string(), inp.header)
                repeat.query_name = repeated_id
                out.write(repeat)
    pysam.sort('-o', out_bam, unsorted)
    pysam.index(out_bam)


def test_repeated_read_id_is_an_error(tmp_path):
    """Two reads of one id at unrelated loci, as a BAM merged from runs gives.  Both
    are assigned without realignment, so the sample's read map is what catches it."""
    manifest, genome, bed = _inputs(tmp_path)
    dup_bam = str(tmp_path / 'reads.dup.bam')
    _bam_with_repeated_id(str(tmp_path / 'reads.bam'), dup_bam, 's1_read', 'a_full')
    with open(manifest, 'w') as fh:
        fh.write(f's1\tA\tb0\t{dup_bam}\n')

    with pytest.raises(FlairInputDataError, match=r'1 read ids name more than one read .*a_full'):
        quantify(manifest=manifest, genome=genome, isoform_bed=bed, output=str(tmp_path / 'out'),
                 threads=1, trust_ends=False, generate_map=True, normalize_ends=False,
                 temp_dir=str(tmp_path / 'tmp'))
