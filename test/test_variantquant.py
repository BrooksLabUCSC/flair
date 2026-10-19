"""flair variantquant: the SNVs of vcf files, and the reads showing them"""
import gzip
from types import SimpleNamespace
import pysam
import pytest
from flair import FlairInputDataError
from flair.flair_variantquant import (NO_ALT, VAR_BIN_SIZE, add_vcf_var, combine_vcf_files, parse_all_bam_files,
                                      parse_single_bam_read, quantvarpos, quantvarpos_cmd, read_bases_at,
                                      read_vars_to_genome_pos_counts, samples_and_site_vcfs, site_positions,
                                      sites_in_span)


def _vcf(path, *records, filters=None, gzipped=False):
    "a vcf of records, (chrom, pos, id, ref, alt), their FILTERs filters, by default PASS"
    filters = filters or ['PASS'] * len(records)
    with (gzip.open if gzipped else open)(path, 'wt') as out:
        out.write('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        for record, vcf_filter in zip(records, filters):
            out.write('\t'.join(record + ('.', vcf_filter, '.')) + '\n')
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

    def write(self, line):
        self.lines.append(line)


def _read(name, seq, md=None, cigar='10M', start=0):
    "a read aligned at chr1:start + 1 (1-based), its differences from the genome given by md"
    header = pysam.AlignmentHeader.from_dict({'SQ': [{'SN': 'chr1', 'LN': 3 * VAR_BIN_SIZE}]})
    read = pysam.AlignedSegment(header)
    read.query_name, read.reference_id, read.reference_start = name, 0, start
    read.cigarstring, read.query_sequence = cigar, seq
    if md is not None:
        read.set_tag('MD', md)
    return read


def test_a_read_showing_any_alt_of_a_multi_allelic_site_is_modified():
    # the genome is AAAAAAAAAA; the site is chr1:5 A>G,T
    sites = add_vcf_var({}, 'chr1', 'A', ['G', 'T'], 5, 'G1')
    temp = _Temp()
    for name, seq, md in (('g', 'AAAAGAAAAA', '4A5'), ('t', 'AAAATAAAAA', '4A5'), ('c', 'AAAACAAAAA', '4A5'),
                          ('ref', 'AAAAAAAAAA', '10')):
        parse_single_bam_read(_read(name, seq, md), temp, sites, site_positions(sites), 0)
    status = {line.split('\t')[1].split('__')[1]: line.rstrip().split('|')[-1] for line in temp.lines}
    assert status == {'g': 'G', 't': 'T', 'c': NO_ALT, 'ref': NO_ALT}


def test_an_snv_and_an_indel_at_one_position_are_not_a_ref_disagreement(tmp_path):
    # a deletion of the A after chr1:100's T
    vcf1 = _vcf(tmp_path / 's1.vcf', ('chr1', '100', '.', 'T', 'C'))
    vcf2 = _vcf(tmp_path / 's2.vcf', ('chr1', '100', '.', 'TA', 'T'))
    assert combine_vcf_files([vcf1, vcf2]) == {('chr1', '100', 'T'): {'C'}}
    assert combine_vcf_files([vcf2, vcf1]) == {('chr1', '100', 'T'): {'C'}}


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
    read_vars_to_genome_pos_counts([str(tempdir / 'chr1.txt')], output, [['s1'], ['s2']], 1, False, vcfvars)
    rows = [line.rstrip('\n').split('\t') for line in open(output + '.var.counts.tsv')]
    assert rows == [['varpos', 'ref', 'alt', 'gene', 'transcript', 's1', 's2'],
                    ['chr1:5', 'A', 'G', 'G1', '', '1;1', '2;0'],
                    ['chr1:5', 'A', 'T', 'G1', '', '1;1', '1;1']]


def _bam(path, reads):
    "an unindexed bam of reads, (name, seq, md), as _read makes them"
    header = {'SQ': [{'SN': 'chr1', 'LN': 1000}]}
    with pysam.AlignmentFile(str(path), 'wb', header=header) as bam:
        for read in reads:
            bam.write(_read(*read))
    return str(path)


SITES_READS = [('g', 'AAAAGAAAAA', '4A5'), ('t', 'AAAATAAAAA', '4A5'), ('ref', 'AAAAAAAAAA', '10')]


def _pos_ref(path):
    "a --pos_ref file of the site chr1:5 A>G,T, as annotate_group_vcf_vars writes it"
    path.write_text('chr1\t0\t5\tA\tG,T\tG1\n')
    return str(path)


def _quantvarpos(tmp_path, **kwargs):
    quantvarpos(bam=_bam(tmp_path / 's1.bam', SITES_READS), pos_ref=_pos_ref(tmp_path / 'sites.tsv'),
                output=str(tmp_path / 'out'), min_coverage=1, output_all=False, **kwargs)
    return open(tmp_path / 'out.var.counts.tsv').read()


def test_a_rerun_to_the_same_output_counts_its_reads_once(tmp_path):
    # the first run's intermediate files are kept, where the second run puts its own
    temp_dir = str(tmp_path / 'tmp')
    first = _quantvarpos(tmp_path, keep_intermediate=True, temp_dir=temp_dir)
    assert first.splitlines()[1:] == ['chr1:5\tA\tG\tG1\t\t2;1', 'chr1:5\tA\tT\tG1\t\t2;1']
    assert _quantvarpos(tmp_path, keep_intermediate=True, temp_dir=temp_dir) == first


def test_samples_read_in_parallel_are_as_read_one_at_a_time(tmp_path):
    vcfvars = add_vcf_var({}, 'chr1', 'A', ['G', 'T'], 5, 'G1')
    samples = [[f's{n}', _bam(tmp_path / f's{n}.bam', SITES_READS[n:])] for n in range(3)]
    results = []
    for threads in (1, 2):
        tempdir = tmp_path / f'tmp{threads}'
        tempdir.mkdir()
        paths = parse_all_bam_files(samples, str(tempdir) + '/', vcfvars, threads)
        results.append([open(path).read() for path in paths])
    assert results[0] == results[1]
    assert [len(lines.splitlines()) for lines in results[0]] == [3, 2, 1]


def test_keeping_intermediate_files_needs_a_temp_dir():
    args = SimpleNamespace(keep_intermediate=True, temp_dir=None, sample_name=None, bam='s1.bam')
    with pytest.raises(FlairInputDataError, match='--temp_dir'):
        quantvarpos_cmd(args)


def test_a_reads_bases_at_sites_are_found_past_clips_insertions_deletions_and_introns():
    # 2 bases soft clipped, chr1:1-3 aligned, an insertion, 4-5, 6-7 deleted, 8-9, 10-19 an
    # intron, 20-21 aligned: bases of the query TT ACG G CA TG GA
    read = _read('r', 'TTACGGCATGGA', cigar='2S3M1I2M2D2M10N2M')
    sites = [1, 3, 4, 6, 9, 15, 21]
    assert list(read_bases_at(read, sites)) == [(1, 'A'), (3, 'G'), (4, 'C'), (9, 'G'), (21, 'A')]


def test_the_sites_of_a_read_crossing_a_bin_are_of_both_bins():
    vcfvars = {}
    for pos in (VAR_BIN_SIZE - 5, VAR_BIN_SIZE, VAR_BIN_SIZE + 1, VAR_BIN_SIZE + 50):
        add_vcf_var(vcfvars, 'chr1', 'A', ['G'], pos, '')
    positions = site_positions(vcfvars)
    # the read covers VAR_BIN_SIZE - 9 to VAR_BIN_SIZE + 10, 1-based
    b = VAR_BIN_SIZE
    assert sites_in_span(positions, 'chr1', b - 10, b + 10) == [b - 5, b, b + 1]
    # a read ending at a bin's first position, 1-based, has its site
    assert sites_in_span(positions, 'chr1', VAR_BIN_SIZE - 10, VAR_BIN_SIZE) == [VAR_BIN_SIZE - 5, VAR_BIN_SIZE]


def test_only_records_passing_their_filters_are_sites(tmp_path):
    vcf = _vcf(tmp_path / 's1.vcf', ('chr1', '100', '.', 'A', 'G'), ('chr1', '200', '.', 'A', 'G'),
               ('chr1', '300', '.', 'A', 'G'), filters=['PASS', 'RefCall', '.'])
    assert sorted(combine_vcf_files([vcf])) == [('chr1', '100', 'A'), ('chr1', '300', 'A')]


def test_a_vcf_of_its_first_five_columns_has_every_record_a_site(tmp_path):
    # as REDIportal's sites, with no FILTER
    vcf = tmp_path / 'sites.vcf'
    vcf.write_text('chr1\t100\t.\tA\tG\nchr1\t200\t.\tT\tC\n')
    assert sorted(combine_vcf_files([str(vcf)])) == [('chr1', '100', 'A'), ('chr1', '200', 'T')]


def test_a_gzipped_vcf_is_read(tmp_path):
    vcf = _vcf(tmp_path / 's1.vcf.gz', ('chr1', '100', '.', 'A', 'G'), gzipped=True)
    assert combine_vcf_files([vcf]) == {('chr1', '100', 'A'): {'G'}}


def _manifest(tmp_path, *lines):
    path = tmp_path / 'manifest.tsv'
    path.write_text(''.join(line + '\n' for line in lines))
    return str(path)


def _files(tmp_path, *names):
    for name in names:
        (tmp_path / name).write_text('')
    return [str(tmp_path / name) for name in names]


def test_a_manifest_of_bams_has_its_sites_from_a_vcf(tmp_path):
    bam1, bam2, vcf = _files(tmp_path, 's1.bam', 's2.bam', 'sites.vcf')
    manifest = _manifest(tmp_path, '# sample, bam', f's1\t{bam1}', '', f's2\t{bam2}')
    assert samples_and_site_vcfs(manifest, None, None, None, vcf) == ([['s1', bam1], ['s2', bam2]], [vcf])
    assert samples_and_site_vcfs(manifest, None, None, 'sites.tsv', None) == ([['s1', bam1], ['s2', bam2]], None)
    with pytest.raises(FlairInputDataError, match='give --vcf or --pos_ref'):
        samples_and_site_vcfs(manifest, None, None, None, None)


def test_a_manifest_of_each_samples_vcf_takes_no_other_sites(tmp_path):
    bam1, bam2, vcf1, vcf2 = _files(tmp_path, 's1.bam', 's2.bam', 's1.vcf', 's2.vcf')
    manifest = _manifest(tmp_path, f's1\t{bam1}\t{vcf1}', f's2\t{bam2}\t{vcf2}')
    assert samples_and_site_vcfs(manifest, None, None, None, None) == ([['s1', bam1], ['s2', bam2]], [vcf1, vcf2])
    with pytest.raises(FlairInputDataError, match='--vcf and --pos_ref can not be given'):
        samples_and_site_vcfs(manifest, None, None, None, vcf1)


def test_a_manifest_with_a_vcf_on_some_lines_is_an_error(tmp_path):
    bam1, bam2, vcf = _files(tmp_path, 's1.bam', 's2.bam', 'sites.vcf')
    # one vcf for all samples on the first line, as variantquant used to take
    for lines in ([f's1\t{bam1}\t{vcf}', f's2\t{bam2}'], [f's1\t{bam1}\t{vcf}', f's2\t{bam2}\t']):
        with pytest.raises(FlairInputDataError, match='give --vcf'):
            samples_and_site_vcfs(_manifest(tmp_path, *lines), None, None, None, None)


def test_a_manifest_naming_a_sample_twice_or_a_missing_file_is_an_error(tmp_path):
    bam1, = _files(tmp_path, 's1.bam')
    with pytest.raises(FlairInputDataError, match='more than once: s1'):
        samples_and_site_vcfs(_manifest(tmp_path, f's1\t{bam1}', f's1\t{bam1}'), None, None, None, 'sites.vcf')
    with pytest.raises(FlairInputDataError, match=r'manifest.tsv:2: bam file .*s2.bam does not exist'):
        samples_and_site_vcfs(_manifest(tmp_path, f's1\t{bam1}', f's2\t{tmp_path}/s2.bam'), None, None, None,
                              'sites.vcf')


def test_a_bams_sample_is_named_for_it_or_by_sample_name(tmp_path):
    bam, = _files(tmp_path, 'ULC0422N.bam')
    assert samples_and_site_vcfs(None, bam, None, None, 'sites.vcf') == ([['ULC0422N', bam]], ['sites.vcf'])
    assert samples_and_site_vcfs(None, bam, 'tumor', None, 'sites.vcf') == ([['tumor', bam]], ['sites.vcf'])
    with pytest.raises(FlairInputDataError, match='give --vcf or --pos_ref'):
        samples_and_site_vcfs(None, bam, None, None, None)


def test_a_vcfs_sites_without_an_isoform_bed_have_no_genes(tmp_path):
    vcf = _vcf(tmp_path / 'sites.vcf', ('chr1', '5', '.', 'A', 'G,T'))
    quantvarpos(bam=_bam(tmp_path / 's1.bam', SITES_READS), vcf=vcf, output=str(tmp_path / 'out'), min_coverage=1,
                output_all=False, keep_intermediate=False)
    assert open(tmp_path / 'out.var.counts.tsv').read().splitlines() == [
        'varpos\tref\talt\tgene\ttranscript\ts1', 'chr1:5\tA\tG\t\t\t2;1', 'chr1:5\tA\tT\t\t\t2;1']


def test_records_of_kept_flags_are_sites_too(tmp_path):
    vcf = _vcf(tmp_path / 's1.vcf', ('chr1', '100', '.', 'A', 'G'), ('chr1', '200', '.', 'A', 'G'),
               ('chr1', '300', '.', 'A', 'G'), ('chr1', '400', '.', 'A', 'G'),
               filters=['PASS', 'RefCall', 'NoCall', 'RefCall;LowQual'])

    def sites(keep_flags):
        return sorted(int(pos) for _, pos, _ in combine_vcf_files([vcf], keep_flags))

    assert sites({'RefCall'}) == [100, 200]
    # a record of more than one flag needs them all kept
    assert sites({'RefCall', 'LowQual'}) == [100, 200, 400]
    assert sites({'all'}) == [100, 200, 300, 400]
