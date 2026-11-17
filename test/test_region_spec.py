import pytest
from flair import FlairInputDataError, SeqRange
from flair.region_spec import parse_region, parse_regions, write_regions_bed


def test_one_based_range_is_closed():
    assert parse_region('chr1:1-1000') == SeqRange('chr1', 0, 1000)
    assert parse_region('chr1:1000-1000') == SeqRange('chr1', 999, 1000)


def test_zero_based_range_is_half_open():
    assert parse_region('chr1 0 1000') == SeqRange('chr1', 0, 1000)
    assert parse_region('chr1\t999\t2000') == SeqRange('chr1', 999, 2000)


def test_coordinates_may_have_commas_and_space_around_them():
    assert parse_region(' chr1:1,000-2,000 ') == SeqRange('chr1', 999, 2000)
    assert parse_region('chr1 1,000 2,000') == SeqRange('chr1', 1000, 2000)


def test_sequence_name_may_contain_colons():
    assert parse_region('HLA-A*01:01:01:01:5-10') == SeqRange('HLA-A*01:01:01:01', 4, 10)


@pytest.mark.parametrize('spec, problem', [('chr1:0-1000', 'the start must be at least 1'),
                                           ('chr1 -1 10', 'not a sequence name with a range'),
                                           ('chr1:1000-999', 'the end must be after the start'),
                                           ('chr1 5 5', 'the end must be after the start'),
                                           ('chr1-1000', 'not a sequence name with a range'),
                                           ('chr1:1-1000-2000', 'not a sequence name with a range')])
def test_invalid_region_is_an_error(spec, problem):
    with pytest.raises(FlairInputDataError) as exc:
        parse_region(spec)
    assert problem in str(exc.value)


def test_no_regions_given_is_none():
    assert parse_regions([], None) is None


def test_regions_of_both_options_are_merged():
    assert parse_regions(['chr1:1-100', 'chr1 50 200'], 'input/partition.region.bed') == \
        [SeqRange('chr1', 0, 200), SeqRange('chr19', 368000, 368100), SeqRange('chr6', 74285000, 74285100)]


def test_empty_region_bed_is_an_error(tmp_path):
    empty = tmp_path / 'empty.bed'
    empty.write_text('')
    with pytest.raises(FlairInputDataError) as exc:
        parse_regions([], str(empty))
    assert 'it must have at least one BED record' in str(exc.value)


def test_regions_are_written_as_bed3(tmp_path):
    bed = str(tmp_path / 'regions.bed')
    write_regions_bed([SeqRange('chr1', 0, 100), SeqRange('chr2', 5, 6)], bed)
    assert open(bed).read() == 'chr1\t0\t100\nchr2\t5\t6\n'
