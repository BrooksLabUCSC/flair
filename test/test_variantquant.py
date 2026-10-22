"""flair variantquant: the SNVs of vcf files, and the reads showing them"""
import gzip
from types import SimpleNamespace
import pysam
import pytest
from flair import FlairInputDataError
from flair.flair_variantquant import (NO_ALT, VAR_BIN_SIZE, add_vcf_var, check_sites_against_genome, exons_of_gtf,
                                      name_sites_by_genes, parse_all_bam_files, parse_single_bam_read,
                                      quantvarpos, quantvarpos_cmd, read_bases_at, read_vars_to_genome_pos_counts,
                                      read_vcf_sites, samples_and_site_vcfs, site_exons, site_positions,
                                      sites_in_span)


def _vcf(path, *records, filters=None, gzipped=False):
    "a vcf of records, (chrom, pos, id, ref, alt), their FILTERs filters, by default PASS"
    filters = filters or ['PASS'] * len(records)
    with (gzip.open if gzipped else open)(path, 'wt') as out:
        out.write('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        for record, vcf_filter in zip(records, filters):
            out.write('\t'.join(record + ('.', vcf_filter, '.')) + '\n')
    return str(path)


def combine_vcf_files(vcfs, keep_flags=frozenset()):
    "(chrom, pos, ref) -> alts of read_vcf_sites' sites"
    return {(chrom, str(pos), ref): set(alts) for (chrom, _), sites in read_vcf_sites(vcfs, keep_flags).items()
            for pos, (ref, alts, _) in sites.items()}


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
            out.write(f'chr1\t{n % 2}__r{n}\t.\t.\t.\t5|G1|{allele}\n')
    vcfvars = add_vcf_var({}, 'chr1', 'A', ['G', 'T'], 5, 'G1')
    output = str(tmp_path / 'out')
    read_vars_to_genome_pos_counts([str(tempdir / 'chr1.txt')], output, [['s1'], ['s2']], 1, False, vcfvars)
    rows = [line.rstrip('\n').split('\t') for line in open(output + '.var.counts.tsv')]
    assert rows == [['var_name', 'is_var', 'varpos', 'ref', 'alt', 'strand', 'phase_set', 'haplotype', 'gene', 's1', 's2'],
                    ['chr1:5:A>G', 'yes', 'chr1:5', 'A', 'G', '.', '.', '.', 'G1', '1', '0'],
                    ['chr1:5:A>G', 'no', 'chr1:5', 'A', 'G', '.', '.', '.', 'G1', '1', '2'],
                    ['chr1:5:A>T', 'yes', 'chr1:5', 'A', 'T', '.', '.', '.', 'G1', '1', '1'],
                    ['chr1:5:A>T', 'no', 'chr1:5', 'A', 'T', '.', '.', '.', 'G1', '1', '1']]


def _bam(path, reads):
    "an unindexed bam of reads, (name, seq, md), as _read makes them"
    header = {'SQ': [{'SN': 'chr1', 'LN': 1000}]}
    with pysam.AlignmentFile(str(path), 'wb', header=header) as bam:
        for read in reads:
            bam.write(_read(*read))
    return str(path)


SITES_READS = [('g', 'AAAAGAAAAA', '4A5'), ('t', 'AAAATAAAAA', '4A5'), ('ref', 'AAAAAAAAAA', '10')]


def _sites_vcf(path):
    "a vcf of the site chr1:5 A>G,T"
    return _vcf(path, ('chr1', '5', '.', 'A', 'G,T'))


def _quantvarpos(tmp_path, **kwargs):
    quantvarpos(bam=_bam(tmp_path / 's1.bam', SITES_READS), vcf=_sites_vcf(tmp_path / 'sites.vcf'),
                output=str(tmp_path / 'out'), min_coverage=1, output_all=False, **kwargs)
    return open(tmp_path / 'out.var.counts.tsv').read()


def test_a_rerun_to_the_same_output_counts_its_reads_once(tmp_path):
    # the first run's intermediate files are kept, where the second run puts its own
    temp_dir = str(tmp_path / 'tmp')
    first = _quantvarpos(tmp_path, keep_intermediate=True, temp_dir=temp_dir)
    assert first.splitlines()[1:] == ['chr1:5:A>G\tyes\tchr1:5\tA\tG\t.\t.\t.\t\t1',
                                      'chr1:5:A>G\tno\tchr1:5\tA\tG\t.\t.\t.\t\t2',
                                      'chr1:5:A>T\tyes\tchr1:5\tA\tT\t.\t.\t.\t\t1',
                                      'chr1:5:A>T\tno\tchr1:5\tA\tT\t.\t.\t.\t\t2']
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
    assert samples_and_site_vcfs(manifest, None, None, vcf) == ([['s1', bam1], ['s2', bam2]], [vcf])
    with pytest.raises(FlairInputDataError, match='give --vcf'):
        samples_and_site_vcfs(manifest, None, None, None)


def test_a_manifest_of_each_samples_vcf_takes_no_other_sites(tmp_path):
    bam1, bam2, vcf1, vcf2 = _files(tmp_path, 's1.bam', 's2.bam', 's1.vcf', 's2.vcf')
    manifest = _manifest(tmp_path, f's1\t{bam1}\t{vcf1}', f's2\t{bam2}\t{vcf2}')
    assert samples_and_site_vcfs(manifest, None, None, None) == ([['s1', bam1], ['s2', bam2]], [vcf1, vcf2])
    with pytest.raises(FlairInputDataError, match='--vcf can not be given'):
        samples_and_site_vcfs(manifest, None, None, vcf1)


def test_a_manifest_with_a_vcf_on_some_lines_is_an_error(tmp_path):
    bam1, bam2, vcf = _files(tmp_path, 's1.bam', 's2.bam', 'sites.vcf')
    # one vcf for all samples on the first line, as variantquant used to take
    for lines in ([f's1\t{bam1}\t{vcf}', f's2\t{bam2}'], [f's1\t{bam1}\t{vcf}', f's2\t{bam2}\t']):
        with pytest.raises(FlairInputDataError, match='give --vcf'):
            samples_and_site_vcfs(_manifest(tmp_path, *lines), None, None, None)


def test_a_manifest_naming_a_sample_twice_or_a_missing_file_is_an_error(tmp_path):
    bam1, = _files(tmp_path, 's1.bam')
    with pytest.raises(FlairInputDataError, match='more than once: s1'):
        samples_and_site_vcfs(_manifest(tmp_path, f's1\t{bam1}', f's1\t{bam1}'), None, None, 'sites.vcf')
    with pytest.raises(FlairInputDataError, match=r'manifest.tsv:2: bam file .*s2.bam does not exist'):
        samples_and_site_vcfs(_manifest(tmp_path, f's1\t{bam1}', f's2\t{tmp_path}/s2.bam'), None, None,
                              'sites.vcf')


def test_a_bams_sample_is_named_for_it_or_by_sample_name(tmp_path):
    bam, = _files(tmp_path, 'ULC0422N.bam')
    assert samples_and_site_vcfs(None, bam, None, 'sites.vcf') == ([['ULC0422N', bam]], ['sites.vcf'])
    assert samples_and_site_vcfs(None, bam, 'tumor', 'sites.vcf') == ([['tumor', bam]], ['sites.vcf'])
    with pytest.raises(FlairInputDataError, match='give --vcf'):
        samples_and_site_vcfs(None, bam, None, None)


def test_a_vcfs_sites_without_an_isoform_bed_have_no_genes(tmp_path):
    vcf = _vcf(tmp_path / 'sites.vcf', ('chr1', '5', '.', 'A', 'G,T'))
    quantvarpos(bam=_bam(tmp_path / 's1.bam', SITES_READS), vcf=vcf, output=str(tmp_path / 'out'), min_coverage=1,
                output_all=False, keep_intermediate=False)
    assert open(tmp_path / 'out.var.counts.tsv').read().splitlines() == ['var_name\tis_var\tvarpos\tref\talt\tstrand\tphase_set\thaplotype\tgene\ts1',
                                                                         'chr1:5:A>G\tyes\tchr1:5\tA\tG\t.\t.\t.\t\t1',
                                                                         'chr1:5:A>G\tno\tchr1:5\tA\tG\t.\t.\t.\t\t2',
                                                                         'chr1:5:A>T\tyes\tchr1:5\tA\tT\t.\t.\t.\t\t1',
                                                                         'chr1:5:A>T\tno\tchr1:5\tA\tT\t.\t.\t.\t\t2']


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


def _isoform_bed(path, *isoforms):
    "a flair isoform bed of isoforms, (name, gene, start, [(exon start, end), ...]) on chr1, + strand"
    with open(path, 'w') as out:
        for name, gene, start, exons in isoforms:
            sizes = ''.join(f'{end - s},' for s, end in exons)
            starts = ''.join(f'{s - start},' for s, _ in exons)
            fields = ['chr1', start, exons[-1][1], name, 20, '+', start, start, '0,0,0', len(exons), sizes, starts,
                      gene, '', '', 20, 1.0, 'PRO', 'basic', '', '', 's1,', '', '', '']
            out.write('\t'.join(str(f) for f in fields) + '\n')
    return str(path)


def _site_genes(exons, *positions):
    "the genes name_sites_by_genes gives sites at positions on chr1, by position"
    vcfvars = {}
    for pos in positions:
        add_vcf_var(vcfvars, 'chr1', 'A', ('G',), pos, '')
    name_sites_by_genes(vcfvars, exons)
    return [vcfvars[('chr1', pos // VAR_BIN_SIZE * VAR_BIN_SIZE)][pos][2] for pos in positions]


def test_a_site_is_of_the_genes_with_an_exon_over_it(tmp_path):
    # exons chr1:101-200 and 301-400, 1-based; another gene's exon 381-500
    bed = _isoform_bed(tmp_path / 'isoforms.bed', ('T1', 'G1', 100, [(100, 200), (300, 400)]),
                       ('T2', 'G2', 380, [(380, 500)]))
    assert _site_genes(site_exons(bed, None), 100, 101, 200, 201, 300, 301, 390, 450, 501) == \
        ['', 'G1', 'G1', '', '', 'G1', 'G1,G2', 'G2', '']


def test_an_exon_reaching_across_a_bin_names_the_sites_of_both():
    b = VAR_BIN_SIZE
    exons = {'chr1': [(b - 10, b + 10, 'G1')]}
    assert _site_genes(exons, b - 20, b - 5, b, b + 5, b + 11) == ['', 'G1', 'G1', 'G1', '']


def test_a_gtfs_exons_are_by_gene_id(tmp_path):
    gtf = tmp_path / 'genes.gtf'
    rows = [('chr1', 'gene', 101, 400, 'gene_id "G1";'),
            ('chr1', 'transcript', 101, 400, 'gene_id "G1"; transcript_id "T1";'),
            ('chr1', 'exon', 101, 200, 'gene_id "G1"; transcript_id "T1"; exon_number 1;'),
            ('chr1', 'exon', 301, 400, 'gene_id "G1"; transcript_id "T1"; exon_number 2;'),
            ('chr1', 'exon', 101, 200, 'gene_id "G1"; transcript_id "T2";'),
            ('chr2', 'exon', 51, 60, 'transcript_id "T3"; gene_id "G2";')]
    gtf.write_text('#!genome-build GRCh38\n' + ''.join(f'{c}\tsrc\t{f}\t{s}\t{e}\t.\t+\t.\t{a}\n'
                                                       for c, f, s, e, a in rows))
    assert exons_of_gtf(str(gtf)) == {'chr1': [(100, 200, 'G1'), (300, 400, 'G1')], 'chr2': [(50, 60, 'G2')]}


def test_supplementary_alignments_are_counted_but_not_secondary(tmp_path):
    vcfvars = add_vcf_var({}, 'chr1', 'A', ['G'], 5, '')
    reads = [_read('p', 'AAAAGAAAAA'), _read('s', 'AAAAGAAAAA'), _read('x', 'AAAAGAAAAA')]
    reads[1].is_supplementary = True
    reads[2].is_secondary = True
    header = {'SQ': [{'SN': 'chr1', 'LN': 3 * VAR_BIN_SIZE}]}
    with pysam.AlignmentFile(str(tmp_path / 's1.bam'), 'wb', header=header) as bam:
        for read in reads:
            bam.write(read)
    tempdir = tmp_path / 'tmp'
    tempdir.mkdir()
    path, = parse_all_bam_files([['s1', str(tmp_path / 's1.bam')]], str(tempdir) + '/', vcfvars, 1)
    assert [line.split('\t')[1] for line in open(path)] == ['0__p', '0__s']


def _gtf(path, *exons):
    "a GTF of exons, (chrom, start, end, strand, gene) 1-based"
    path.write_text(''.join(f'{c}\tsrc\texon\t{s}\t{e}\t.\t{strand}\t.\tgene_id "{g}"; transcript_id "{g}.1";\n'
                            for c, s, e, strand, g in exons))
    return str(path)


def _bed12(path, *isoforms):
    "a BED12 of isoforms, (name, start, strand, [(exon start, end), ...]) on chr1, 0-based"
    with open(path, 'w') as out:
        for name, start, strand, exons in isoforms:
            sizes = ''.join(f'{end - s},' for s, end in exons)
            starts = ''.join(f'{s - start},' for s, _ in exons)
            fields = ('chr1', start, exons[-1][1], name, 0, strand, start, start, '0,0,0', len(exons), sizes, starts)
            out.write('\t'.join(str(f) for f in fields) + '\n')
    return str(path)


def test_a_bed12s_isoforms_are_of_the_gtf_genes_they_overlap_most(tmp_path):
    # the GTF's G1 on + has exons 101-200 and 301-400, G2 on - 151-450; the bed's isoform on +
    # has a novel exon 501-600, past them, and T2 overlaps no gene
    gtf = _gtf(tmp_path / 'genes.gtf', ('chr1', 101, 200, '+', 'G1'), ('chr1', 301, 400, '+', 'G1'),
               ('chr1', 151, 450, '-', 'G2'))
    bed = _bed12(tmp_path / 'isoforms.bed', ('T1', 100, '+', [(100, 200), (300, 400), (500, 600)]),
                 ('T2', 2000, '+', [(2000, 2100)]))
    exons = site_exons(bed, gtf)
    assert exons == {'chr1': [(100, 200, 'G1'), (300, 400, 'G1'), (500, 600, 'G1')]}
    # the novel exon's site is G1's, which the GTF alone misses, and G2 names no site of the bed
    assert _site_genes(exons, 160, 550, 2050) == ['G1', 'G1', '']
    assert _site_genes(site_exons(None, gtf), 160, 550) == ['G1,G2', '']


def test_a_flair_beds_genes_are_its_own_and_a_bed12_alone_names_no_sites(tmp_path, caplog):
    gtf = _gtf(tmp_path / 'genes.gtf', ('chr1', 101, 200, '+', 'G1'))
    flair_bed = _isoform_bed(tmp_path / 'isoforms.bed', ('T1', 'FLG1', 100, [(100, 200)]))
    assert site_exons(flair_bed, gtf) == {'chr1': [(100, 200, 'FLG1')]}
    bed = _bed12(tmp_path / 'plain.bed', ('T1', 100, '+', [(100, 200)]))
    assert site_exons(bed, None) == {}
    assert '1 of the 1 isoforms' in caplog.text


def test_sites_whose_ref_is_not_the_genomes_are_left_out(tmp_path, caplog):
    fasta = tmp_path / 'genome.fa'
    fasta.write_text('>chr1\nACGTACGTAC\n')
    pysam.faidx(str(fasta))
    # chr1:1 is A, as is its record; 2, C, recorded G>A, its complement, as of the transcript's
    # strand; 3, G, recorded A>G; 12 is past chr1's end; chr2 is not in the genome
    vcf = _vcf(tmp_path / 'sites.vcf', ('chr1', '1', '.', 'A', 'G'), ('chr1', '2', '.', 'G', 'A'),
               ('chr1', '3', '.', 'A', 'G'), ('chr1', '12', '.', 'A', 'G'), ('chr2', '5', '.', 'A', 'G'))
    vcfvars = read_vcf_sites([vcf])
    check_sites_against_genome(vcfvars, str(fasta))
    assert vcfvars == {('chr1', 0): {1: ('A', ('G',), '')}}
    assert "3 vcf sites whose REF is not the genome's base left out, 1 of them its complement" in caplog.text
    assert 'chr1:12 A>G at None' in caplog.text
    assert 'chr1:2 G>A at C' in caplog.text
    assert '1 vcf sites on sequences not in the genome left out' in caplog.text


def test_trusting_strand_counts_each_strand_of_reads_in_a_row_of_its_own(tmp_path):
    # chr1:5 A>G: + reads g and ref, - reads g2, ref2 and ref3, and low, of mapping quality 0
    reads = [_read('g', 'AAAAGAAAAA'), _read('ref', 'AAAAAAAAAA'), _read('g2', 'AAAAGAAAAA'),
             _read('ref2', 'AAAAAAAAAA'), _read('ref3', 'AAAAAAAAAA'), _read('low', 'AAAAGAAAAA')]
    for read in reads:
        read.mapping_quality = 60
    for read in reads[2:]:
        read.is_reverse = True
    reads[5].mapping_quality = 0
    header = {'SQ': [{'SN': 'chr1', 'LN': 3 * VAR_BIN_SIZE}]}
    with pysam.AlignmentFile(str(tmp_path / 's1.bam'), 'wb', header=header) as bam:
        for read in reads:
            bam.write(read)
    vcf = _vcf(tmp_path / 'sites.vcf', ('chr1', '5', '.', 'A', 'G'))

    def counts(**options):
        quantvarpos(bam=str(tmp_path / 's1.bam'), vcf=vcf, output=str(tmp_path / 'out'), min_coverage=1,
                    output_all=False, keep_intermediate=False, **options)
        return sorted(open(tmp_path / 'out.var.counts.tsv').read().splitlines())

    assert counts() == ['chr1:5:A>G\tno\tchr1:5\tA\tG\t.\t.\t.\t\t3', 'chr1:5:A>G\tyes\tchr1:5\tA\tG\t.\t.\t.\t\t3', 'var_name\tis_var\tvarpos\tref\talt\tstrand\tphase_set\thaplotype\tgene\ts1']
    assert counts(quality=1) == ['chr1:5:A>G\tno\tchr1:5\tA\tG\t.\t.\t.\t\t3', 'chr1:5:A>G\tyes\tchr1:5\tA\tG\t.\t.\t.\t\t2', 'var_name\tis_var\tvarpos\tref\talt\tstrand\tphase_set\thaplotype\tgene\ts1']
    assert counts(quality=1, trust_strand=True) == ['chr1:5:A>G(+)\tno\tchr1:5\tA\tG\t+\t.\t.\t\t1',
                                                    'chr1:5:A>G(+)\tyes\tchr1:5\tA\tG\t+\t.\t.\t\t1',
                                                    'chr1:5:A>G(-)\tno\tchr1:5\tA\tG\t-\t.\t.\t\t2',
                                                    'chr1:5:A>G(-)\tyes\tchr1:5\tA\tG\t-\t.\t.\t\t1',
                                                    'var_name\tis_var\tvarpos\tref\talt\tstrand\tphase_set\thaplotype\tgene\ts1']


def _tagged_bam(path, reads):
    "a bam of reads, (name, seq, HP, PS, reverse), HP None for an untagged read"
    header = {'SQ': [{'SN': 'chr1', 'LN': 3 * VAR_BIN_SIZE}]}
    with pysam.AlignmentFile(str(path), 'wb', header=header) as bam:
        for name, seq, hp, ps, reverse in reads:
            read = _read(name, seq)
            read.mapping_quality, read.is_reverse = 60, reverse
            if hp is not None:
                read.set_tag('HP', hp)
                read.set_tag('PS', ps)
            bam.write(read)
    return str(path)


def test_haplotypes_count_the_reads_of_each_phase_set_and_haplotype_in_a_row_of_their_own(tmp_path):
    # chr1:5 A>G: haplotype 1 of phase set 100 all A, haplotype 2 G but one, two reads untagged
    bam = _tagged_bam(tmp_path / 's1.bam', [('a1', 'AAAAAAAAAA', 1, 100, False), ('a2', 'AAAAAAAAAA', 1, 100, True),
                                            ('g1', 'AAAAGAAAAA', 2, 100, False), ('g2', 'AAAAGAAAAA', 2, 100, False),
                                            ('a3', 'AAAAAAAAAA', 2, 100, True), ('u1', 'AAAAGAAAAA', None, None, False),
                                            ('u2', 'AAAAAAAAAA', None, None, True)])
    vcf = _vcf(tmp_path / 'sites.vcf', ('chr1', '5', '.', 'A', 'G'))

    def counts(**options):
        # an allele's rows are of every group, haplotype 1, of no read showing G, too
        quantvarpos(bam=bam, vcf=vcf, output=str(tmp_path / 'out'), min_coverage=1, output_all=False,
                    keep_intermediate=False, **options)
        return sorted(open(tmp_path / 'out.var.counts.tsv').read().splitlines())

    assert counts(haplotypes=True) == ['chr1:5:A>G\tno\tchr1:5\tA\tG\t.\t.\t.\t\t1',
                                       'chr1:5:A>G\tyes\tchr1:5\tA\tG\t.\t.\t.\t\t1',
                                       'chr1:5:A>G|100|1\tno\tchr1:5\tA\tG\t.\t100\t1\t\t2',
                                       'chr1:5:A>G|100|1\tyes\tchr1:5\tA\tG\t.\t100\t1\t\t0',
                                       'chr1:5:A>G|100|2\tno\tchr1:5\tA\tG\t.\t100\t2\t\t1',
                                       'chr1:5:A>G|100|2\tyes\tchr1:5\tA\tG\t.\t100\t2\t\t2',
                                       'var_name\tis_var\tvarpos\tref\talt\tstrand\tphase_set\thaplotype\tgene\ts1']
    # with strand, each strand of each haplotype
    assert counts(haplotypes=True, trust_strand=True)[:3] == ['chr1:5:A>G(+)\tno\tchr1:5\tA\tG\t+\t.\t.\t\t0', 'chr1:5:A>G(+)\tyes\tchr1:5\tA\tG\t+\t.\t.\t\t1', 'chr1:5:A>G(+)|100|1\tno\tchr1:5\tA\tG\t+\t100\t1\t\t1']


def test_samples_sharing_no_phase_set_are_warned_of(tmp_path, caplog):
    bam1 = _tagged_bam(tmp_path / 's1.bam', [('a', 'AAAAGAAAAA', 1, 100, False)])
    bam2 = _tagged_bam(tmp_path / 's2.bam', [('b', 'AAAAGAAAAA', 1, 200, False)])
    bam3 = _tagged_bam(tmp_path / 's3.bam', [('c', 'AAAAGAAAAA', 2, 100, False)])
    vcf = _vcf(tmp_path / 'sites.vcf', ('chr1', '5', '.', 'A', 'G'))
    manifest = _manifest(tmp_path, f's1\t{bam1}', f's3\t{bam3}')
    quantvarpos(manifest=manifest, vcf=vcf, output=str(tmp_path / 'out'), min_coverage=1, output_all=False,
                keep_intermediate=False, haplotypes=True)
    assert 'share no phase set' not in caplog.text
    manifest = _manifest(tmp_path, f's1\t{bam1}', f's2\t{bam2}')
    quantvarpos(manifest=manifest, vcf=vcf, output=str(tmp_path / 'out'), min_coverage=1, output_all=False,
                keep_intermediate=False, haplotypes=True)
    assert 'samples s1 and s2 share no phase set' in caplog.text
