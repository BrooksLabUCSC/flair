"""
Tests for the sorted BED output of the io_utils module.
"""
import gzip
from flair.io_utils import open_sorted_bed, sort_bed_file
from flair.pycbio.hgdata.bed import Bed, BedReader

# out of order on chrom, on start, and two rows at one position
ROWS = (Bed('chr2', 100, 200, 'b'),
        Bed('chr10', 5, 9, 'd'),
        Bed('chr2', 9, 20, 'a'),
        Bed('chr2', 9, 30, 'c'))

# C collation, so chr10 comes before chr2, and the two chr2:9 rows keep the order
# they were written in
SORTED_NAMES = ['d', 'a', 'c', 'b']

def write_rows(path):
    with open_sorted_bed(path) as fh:
        for bed in ROWS:
            bed.write(fh)

def read_names(path):
    return [bed.name for bed in BedReader(path)]

def test_sorted_as_written(tmp_path):
    bed_file = str(tmp_path / 'x.bed')
    write_rows(bed_file)
    assert read_names(bed_file) == SORTED_NAMES

def test_sorted_compressed(tmp_path):
    "the suffix chooses the compression, the order is the same"
    bed_file = str(tmp_path / 'x.bed.gz')
    write_rows(bed_file)
    assert open(bed_file, 'rb').read(2) == b'\x1f\x8b'
    assert read_names(bed_file) == SORTED_NAMES

def test_sort_written_file(tmp_path):
    bed_file = str(tmp_path / 'x.bed')
    with open(bed_file, 'w') as fh:
        for bed in ROWS:
            bed.write(fh)
    sort_bed_file(bed_file)
    assert read_names(bed_file) == SORTED_NAMES

def test_sort_written_file_compressed(tmp_path):
    bed_file = str(tmp_path / 'x.bed.gz')
    with gzip.open(bed_file, 'wt') as fh:
        for bed in ROWS:
            bed.write(fh)
    sort_bed_file(bed_file)
    assert open(bed_file, 'rb').read(2) == b'\x1f\x8b'
    assert read_names(bed_file) == SORTED_NAMES
