"""
tests of run_checks.check_annotation, the annotation checked where it is parsed
"""
import pytest
import pysam
from flair import FlairInputDataError, SeqRange
from flair.gtf_io import gtf_data_parser
from flair.run_checks import check_annotation

GTF_LINES = ('{chrom}\tHAVANA\ttranscript\t1000\t2000\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
             '{chrom}\tHAVANA\texon\t1000\t2000\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n')


def _genome(tmp_path, names):
    "a FASTA of named sequences, indexed as a run's genome is"
    fasta = tmp_path / 'genome.fa'
    fasta.write_text(''.join(f'>{name}\n' + ('ACGT' * 1000) + '\n' for name in names))
    pysam.faidx(str(fasta))
    return str(fasta)


def _annotation(tmp_path, chrom, name='annot.gtf'):
    gtf = tmp_path / name
    gtf.write_text(GTF_LINES.format(chrom=chrom))
    return gtf_data_parser(str(gtf))


def test_annotation_of_the_genome(tmp_path):
    check_annotation(_annotation(tmp_path, 'chr1'), _genome(tmp_path, ['chr1', 'chr2']))


def test_annotation_named_differently(tmp_path):
    with pytest.raises(FlairInputDataError, match='has none of the sequences of the genome'):
        check_annotation(_annotation(tmp_path, '1'), _genome(tmp_path, ['chr1']))


def test_annotation_without_transcripts(tmp_path):
    empty = tmp_path / 'empty.gtf'
    empty.write_text('')
    with pytest.raises(FlairInputDataError, match='has no transcripts'):
        check_annotation(gtf_data_parser(str(empty)), _genome(tmp_path, ['chr1']))


def test_annotation_of_no_transcript_in_the_regions(tmp_path):
    "a parse restricted to regions keeps only their transcripts, so having none is not an error"
    gtf = tmp_path / 'annot.gtf'
    gtf.write_text(GTF_LINES.format(chrom='chr1'))
    gtf_data = gtf_data_parser(str(gtf), regions=[SeqRange('chr1', 5000, 6000)])
    check_annotation(gtf_data, _genome(tmp_path, ['chr1']), region_filtered=True)
