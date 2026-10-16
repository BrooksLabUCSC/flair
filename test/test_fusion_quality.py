"""
flair fusion --quality: genome-aligned chimeric reads with a piece below the minimum
MAPQ are left out of fusion calling.
"""
import pysam
from flair.identify_prelim_fusions import is_low_quality


def _write_bam(path, records):
    "records: (name, flag, reference position, MAPQ)"
    header = {'HD': {'VN': '1.6'}, 'SQ': [{'SN': 'chr12', 'LN': 200000}]}
    with pysam.AlignmentFile(path, 'wb', header=header) as out:
        for name, flag, start, mapq in records:
            a = pysam.AlignedSegment(out.header)
            a.query_name, a.flag, a.reference_id, a.reference_start = name, flag, 0, start
            a.cigarstring = '100M'
            a.query_sequence = 'A' * 100
            a.mapping_quality = mapq
            out.write(a)


def test_reads_with_a_low_quality_piece(tmp_path):
    "TAOK3--RPS2 in MCF7: a short partner piece at MAPQ 1 beside a MAPQ 60 one"
    bam = str(tmp_path / 'chim.bam')
    _write_bam(bam, [('low_partner', 0, 1000, 60), ('low_partner', 2048, 5000, 1),
                     ('both_good', 0, 1000, 60), ('both_good', 2048, 5000, 60),
                     ('at_threshold', 0, 1000, 10), ('at_threshold', 2048, 5000, 60),
                     ('low_secondary', 0, 1000, 60), ('low_secondary', 256, 5000, 0)])

    def low_reads(min_mapq):
        with pysam.AlignmentFile(bam) as aligns:
            return {a.query_name for a in aligns if is_low_quality(a, min_mapq)}
    assert low_reads(10) == {'low_partner'}
    assert low_reads(0) == set()
