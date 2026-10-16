"""
flair fusion passes transcriptome the synthetic-genome junctions with 2 or more
reads that are canonical, and for a breakpoint with no canonical junction across it,
a non-canonical one only if most reads crossing the breakpoint agree on it.
"""
import pysam
from flair.pycbio.hgdata.bed import Bed
from flair.flair_fusion import select_junctions, count_breakpoint_reads

BP = {'A--B': 1000, 'C--D': 5000}


def _junc(chrom, start, end, reads, motif, strand='+'):
    return Bed(chrom, start, end, name=f'sj1_{motif}', score=reads, strand=strand)


def _selected(junctions, crossing):
    return [(row[0], int(row[1]), int(row[2])) for row in select_junctions(junctions, BP, crossing)]


def test_canonical_breakpoint_junction_kept_without_share():
    juncs = [_junc('A--B', 900, 1200, 2, 'GT/AG')]
    assert _selected(juncs, {'A--B': 20}) == [('A--B', 900, 1200)]


def test_noncanonical_breakpoint_junction_needs_most_crossing_reads():
    "TAOK3--RPS2 in MCF7: 3 of 21 reads crossing the breakpoint on one gc/gg junction"
    scattered = [_junc('A--B', 900, 1200, 3, 'gc/gg', strand='.')]
    agreed = [_junc('A--B', 900, 1200, 27, 'tt/ga', strand='.')]
    assert _selected(scattered, {'A--B': 21}) == []
    assert _selected(agreed, {'A--B': 29}) == [('A--B', 900, 1200)]


def test_noncanonical_junction_not_used_when_a_canonical_one_crosses():
    juncs = [_junc('A--B', 900, 1200, 2, 'GT/AG'), _junc('A--B', 950, 1100, 5, 'gc/gg', strand='.')]
    assert _selected(juncs, {'A--B': 7}) == [('A--B', 900, 1200)]


def test_junctions_away_from_breakpoint():
    "canonical junctions anywhere are kept, non-canonical ones off the breakpoint are not"
    juncs = [_junc('C--D', 100, 400, 4, 'GT/AG'), _junc('C--D', 600, 800, 4, 'ga/ag', strand='.'),
             _junc('C--D', 4800, 5300, 1, 'GT/AG')]
    assert _selected(juncs, {'C--D': 4}) == [('C--D', 100, 400)]


def test_count_breakpoint_reads(tmp_path):
    bam = str(tmp_path / 'synth.bam')
    header = {'HD': {'VN': '1.6'}, 'SQ': [{'SN': 'A--B', 'LN': 5000}, {'SN': 'C--D', 'LN': 9000}]}
    with pysam.AlignmentFile(bam, 'wb', header=header) as out:
        for name, ref, start, cigar, flag in (('across', 0, 800, '100M300N100M', 0),
                                              ('before', 0, 500, '200M', 0),
                                              ('supplementary', 0, 800, '100M300N100M', 2048),
                                              ('other', 1, 4900, '200M', 16)):
            a = pysam.AlignedSegment(out.header)
            a.query_name, a.flag, a.reference_id, a.reference_start, a.cigarstring = name, flag, ref, start, cigar
            a.query_sequence = 'A' * a.query_alignment_length
            out.write(a)
    assert count_breakpoint_reads(bam, BP) == {'A--B': 1, 'C--D': 1}
