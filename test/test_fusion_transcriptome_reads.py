"""
flair fusion aligns to annotated transcripts only the reads that can show a fusion the
genome alignment misses: unmapped, chimeric, or with a partner-sized stretch clipped.
A read aligned to the genome in one piece is left out, as when chimeric in the
transcript alignment it is consecutive sequence put on two annotated genes.
"""
import re
import pysam
from flair.identify_prelim_fusions import scan_genome_alignment

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


RECORDS = [
    ('whole', 0, '300M', None),
    ('spliced_whole', 0, '150M5000N150M', None),  # one spliced alignment, as two overlapping genes
    ('small_clips', 0, '20S270M10S', None),
    ('clipped_partner', 0, '211S89M', None),  # MYO6--SENP6: the MYO6 exon 1 piece left clipped
    ('clipped_3prime', 16, '250M50S', None),
    ('genome_chimeric', 0, '200M100S', 'chr6,50000,+,200S100M,60,0;'),
    ('genome_chimeric', 2048, '200H100M', 'chr6,1000,+,200M100S,60,0;'),
    ('whole_secondary', 256, '100M', None),
    ('unmapped', 4, None, None),
]


def test_reads_to_align_to_transcripts(tmp_path):
    bam, fa, chim = str(tmp_path / 'genome.bam'), str(tmp_path / 'reads.fa'), str(tmp_path / 'chim.bam')
    _write_bam(bam, RECORDS)
    needed = scan_genome_alignment(bam, MIN_CLIP, chim_bam=chim, reads_fa=fa)
    assert needed == {'clipped_partner', 'clipped_3prime', 'genome_chimeric', 'unmapped'}
    # each read once, from its primary record, oriented as sequenced
    seqs = dict(zip(*[iter(line.strip() for line in open(fa))] * 2))
    assert sorted(seqs) == ['>clipped_3prime', '>clipped_partner', '>genome_chimeric', '>unmapped']
    assert seqs['>clipped_3prime'] == 'T' * 300 and seqs['>genome_chimeric'] == 'A' * 300
    # both alignments of the chimeric read, and nothing else
    with pysam.AlignmentFile(chim) as aligns:
        assert [(a.query_name, a.is_supplementary) for a in aligns] == [('genome_chimeric', False), ('genome_chimeric', True)]


def test_without_outputs(tmp_path):
    bam = str(tmp_path / 'genome.bam')
    _write_bam(bam, RECORDS)
    assert scan_genome_alignment(bam, MIN_CLIP) == {'clipped_partner', 'clipped_3prime', 'genome_chimeric', 'unmapped'}
