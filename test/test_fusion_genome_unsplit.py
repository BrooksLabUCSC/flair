"""
flair fusion leaves out reads chimeric in the transcriptome alignment but aligned to
the genome in one piece, keeping those whose genome alignment leaves a partner-sized
stretch clipped.
"""
import re
import pysam
from flair.identify_prelim_fusions import genome_unsplit_reads, chimeric_read_names

MIN_CLIP = 40


def _write_bam(path, records):
    "records: (name, flag, cigar, SA tag or None)"
    header = {'HD': {'VN': '1.6'}, 'SQ': [{'SN': 'chr6', 'LN': 200000}]}
    with pysam.AlignmentFile(path, 'wb', header=header) as out:
        for name, flag, cigar, sa in records:
            a = pysam.AlignedSegment(out.header)
            a.query_name, a.flag = name, flag
            if not flag & 4:
                a.reference_id, a.reference_start, a.cigarstring = 0, 1000, cigar
                a.mapping_quality = 60
            # hard-clipped bases are not in the record
            a.query_sequence = 'A' * (sum(int(n) for n, op in re.findall(r'(\d+)([MIS=X])', cigar)) if cigar else 300)
            if sa:
                a.set_tag('SA', sa)
            out.write(a)


def test_reads_whole_in_the_genome(tmp_path):
    bam = str(tmp_path / 'genome.bam')
    _write_bam(bam, [
        ('whole', 0, '300M', None),
        ('spliced_whole', 0, '150M5000N150M', None),  # one spliced alignment, as two overlapping genes
        ('small_clips', 0, '20S270M10S', None),
        ('small_hard_clip', 16, '39H261M', None),
        ('clipped_partner', 0, '211S89M', None),  # MYO6--SENP6: the MYO6 exon 1 piece left clipped
        ('clipped_3prime', 16, '250M50S', None),
        ('genome_chimeric', 0, '200M100S', 'chr6,50000,+,200S100M,60,0;'),
        ('unmapped', 4, None, None),
        ('not_asked', 0, '300M', None),
    ])
    names = {'whole', 'spliced_whole', 'small_clips', 'small_hard_clip', 'clipped_partner',
             'clipped_3prime', 'genome_chimeric', 'unmapped'}
    assert genome_unsplit_reads(bam, names, MIN_CLIP) == {'whole', 'spliced_whole', 'small_clips', 'small_hard_clip'}


def test_chimeric_read_names(tmp_path):
    bam = str(tmp_path / 'transcriptome.bam')
    _write_bam(bam, [('chimeric', 0, '200M100S', 'chr6,50000,+,200S100M,60,0;'),
                     ('chimeric', 2048, '200H100M', 'chr6,1000,+,200M100S,60,0;'),
                     ('whole', 0, '300M', None)])
    assert chimeric_read_names(bam) == {'chimeric'}
