"""flair variantquant: the SNVs of vcf files, and the reads showing them"""
import pysam
import pytest
from flair import FlairInputDataError
from flair.flair_variantquant import (NO_ALT, add_vcf_var, combine_vcf_files, parse_single_bam_read,
                                      read_vars_to_genome_pos_counts)


def _vcf(path, *records):
    "a vcf of records, (chrom, pos, id, ref, alt)"
    with open(path, 'w') as out:
        out.write('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        for record in records:
            out.write('\t'.join(record + ('.', 'PASS', '.')) + '\n')
    return str(path)


def test_each_alt_of_a_multi_allelic_site_is_an_snv(tmp_path):
    vcf = _vcf(tmp_path / 's1.vcf', ('chr1', '100', '.', 'A', 'G,T'))
    assert combine_vcf_files([vcf]) == {('chr1', '100', 'A'): {'G', 'T'}}


def test_indels_and_longer_substitutions_are_left_out(tmp_path, caplog):
    vcf = _vcf(tmp_path / 's1.vcf', ('chr1', '100', '.', 'A', 'G,AT'), ('chr1', '200', '.', 'AT', 'A'),
               ('chr1', '300', '.', 'AC', 'GT'), ('chr1', '400', '.', 'C', '*'))
    assert combine_vcf_files([vcf]) == {('chr1', '100', 'A'): {'G'}}
    assert '4 vcf alleles that are not SNVs' in caplog.text


def test_a_site_of_two_vcfs_has_the_alts_of_both_whatever_its_ids(tmp_path):
    vcf1 = _vcf(tmp_path / 's1.vcf', ('chr1', '100', 'rs1', 'A', 'G'))
    vcf2 = _vcf(tmp_path / 's2.vcf', ('chr1', '100', '.', 'A', 'T'))
    assert combine_vcf_files([vcf1, vcf2]) == {('chr1', '100', 'A'): {'G', 'T'}}


class _Temp:
    "the temp file of parse_single_bam_read's records, in memory"
    def __init__(self):
        self.lines = []

    def open(self, refname):
        return self

    def write(self, line):
        self.lines.append(line)


def _read(name, seq, md):
    "a read aligned at chr1:1-10 (1-based), its differences from the genome given by md"
    header = pysam.AlignmentHeader.from_dict({'SQ': [{'SN': 'chr1', 'LN': 1000}]})
    read = pysam.AlignedSegment(header)
    read.query_name, read.reference_id, read.reference_start = name, 0, 0
    read.cigarstring, read.query_sequence = '10M', seq
    read.set_tag('MD', md)
    return read


def test_a_read_showing_any_alt_of_a_multi_allelic_site_is_modified():
    # the genome is AAAAAAAAAA; the site is chr1:5 A>G,T
    sites = {5: ('A', ['G', 'T'], 'G1')}
    temp = _Temp()
    for name, seq, md in (('g', 'AAAAGAAAAA', '4A5'), ('t', 'AAAATAAAAA', '4A5'), ('c', 'AAAACAAAAA', '4A5'),
                          ('ref', 'AAAAAAAAAA', '10')):
        parse_single_bam_read(_read(name, seq, md), temp, sites, 0)
    status = {line.split('\t')[1].split('__')[1]: line.rstrip().split('|')[-1] for line in temp.lines}
    assert status == {'g': 'G', 't': 'T', 'c': NO_ALT, 'ref': NO_ALT}


def test_vcfs_disagreeing_on_a_sites_ref_are_an_error(tmp_path):
    vcf1 = _vcf(tmp_path / 's1.vcf', ('chr1', '100', '.', 'A', 'G'))
    vcf2 = _vcf(tmp_path / 's2.vcf', ('chr1', '100', '.', 'C', 'T'))
    with pytest.raises(FlairInputDataError, match='chr1:100'):
        combine_vcf_files([vcf1, vcf2])


def test_each_allele_of_a_site_is_counted_of_the_reads_covering_it(tmp_path):
    # chr1:5 A>G,T: of 4 reads, one G, two T, one A; G has no read in the second sample
    tempdir = tmp_path / 'tmp'
    tempdir.mkdir()
    with open(tempdir / 'chr1.txt', 'w') as out:
        for n, allele in enumerate(['G', 'T', 'T', NO_ALT]):
            out.write(f'chr1\t{n % 2}__r{n}\t5|G1|{allele}\n')
    vcfvars = add_vcf_var({}, 'chr1', 'A', ['G', 'T'], 5, 'G1')
    output = str(tmp_path / 'out')
    read_vars_to_genome_pos_counts(['chr1'], str(tempdir) + '/', output, [['s1'], ['s2']], 1, False, vcfvars)
    rows = [line.rstrip('\n').split('\t') for line in open(output + '.var.counts.tsv')]
    assert rows == [['varpos', 'ref', 'alt', 'gene', 'transcript', 's1', 's2'],
                    ['chr1:5', 'A', 'G', 'G1', '', '1;1', '2;0'],
                    ['chr1:5', 'A', 'T', 'G1', '', '1;1', '1;1']]
