"""Random access by isoform to a read map (read_map_tsv), as flair quantify
--generate_map writes it, so that only some isoforms' reads are read.

A read map has a row per isoform and can name tens of millions of reads, too many
to hold at once.  ReadMapIndex reads the file through once, keeping only the offset
of each isoform's row, and then reads the rows of the isoforms asked for.

The file may be plain, bgzip compressed, which is read by its virtual offsets, or
plain gzip, which cannot be read from an offset and is decompressed once into a
temporary directory first.  flair quantify writes bgzip under --compress.
"""
import os
import gzip
import shutil
from flair import FlairInputDataError
from flair.read_map_tsv import COLUMNS

_GZIP_MAGIC = b'\x1f\x8b'
_HEADER = '\t'.join(COLUMNS)


def is_bgzip(path):
    """whether a file is bgzip compressed: gzip with the extra field of a BGZF block,
    subfield BC"""
    with open(path, 'rb') as fh:
        head = fh.read(18)
    return len(head) == 18 and head[:2] == _GZIP_MAGIC and head[3] & 4 and head[12:14] == b'BC'


def is_gzip(path):
    with open(path, 'rb') as fh:
        return fh.read(2) == _GZIP_MAGIC


class _PlainLines:
    "lines of a plain file by byte offset"
    def __init__(self, path):
        self.fh = open(path, 'rb')

    def tell(self):
        return self.fh.tell()

    def seek(self, offset):
        self.fh.seek(offset)

    def readline(self):
        return self.fh.readline().rstrip(b'\r\n')

    def close(self):
        self.fh.close()


class _BgzipLines:
    "lines of a bgzip file by virtual offset"
    def __init__(self, path):
        from pysam.libcbgzf import BGZFile
        self.fh = BGZFile(path, 'rb')

    def tell(self):
        return self.fh.tell()

    def seek(self, offset):
        self.fh.seek(offset)

    def readline(self):
        return self.fh.readline().rstrip(b'\r\n')

    def close(self):
        self.fh.close()


def _decompressed(path, temp_dir):
    "a plain copy of a plain-gzip file, in temp_dir"
    plain = os.path.join(temp_dir, os.path.basename(path)[:-len('.gz')] if path.endswith('.gz')
                         else os.path.basename(path) + '.plain')
    with gzip.open(path, 'rb') as src, open(plain, 'wb') as dst:
        shutil.copyfileobj(src, dst)
    return plain


class ReadMapIndex:
    """The offset of each isoform's row in a read map.  Pickles small, so it can be
    handed to worker processes, each of which opens the file itself (open_reader)."""

    def __init__(self, read_map_tsv, temp_dir):
        self.read_map_tsv = read_map_tsv
        if is_bgzip(read_map_tsv):
            self.path, self.bgzip = read_map_tsv, True
        elif is_gzip(read_map_tsv):
            self.path, self.bgzip = _decompressed(read_map_tsv, temp_dir), False
        else:
            self.path, self.bgzip = read_map_tsv, False
        self.offsets = self._index()

    def _open(self):
        return _BgzipLines(self.path) if self.bgzip else _PlainLines(self.path)

    def _index(self):
        offsets = {}
        lines = self._open()
        try:
            first = True
            while True:
                offset = lines.tell()
                line = lines.readline()
                if not line:
                    break
                text = line.decode()
                if first and text == _HEADER:
                    first = False
                    continue
                first = False
                name, sep, _ = text.partition('\t')
                if not sep:
                    raise FlairInputDataError(f'{self.read_map_tsv}: expected a row of an isoform and its reads, '
                                              f'tab separated, found: {text[:100]}')
                if name in offsets:
                    raise FlairInputDataError(f'{self.read_map_tsv}: isoform {name} has two rows')
                offsets[name] = offset
        finally:
            lines.close()
        return offsets

    def open_reader(self):
        return ReadMapReader(self)


class ReadMapReader:
    "reads the rows of chosen isoforms through a ReadMapIndex; one per process"

    def __init__(self, index):
        self.index = index
        self.lines = index._open()

    def assigned(self, isoforms):
        "{read name: isoform} for the reads of the isoforms that have a row"
        out = {}
        for iso in isoforms:
            offset = self.index.offsets.get(iso)
            if offset is None:
                continue
            self.lines.seek(offset)
            _, _, reads = self.lines.readline().decode().partition('\t')
            for read in reads.split(','):
                if read:
                    out[read] = iso
        return out

    def close(self):
        self.lines.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
