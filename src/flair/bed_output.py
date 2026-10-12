"""Final BED output, sorted by position.

A BED that flair hands to the user is sorted by chrom and start, so that a run is
comparable with another and so that bedToBigBed and genome browsers get the order
they expect.

Rows that are all in memory are sorted with Bed.genome_sort_key.  A stream of
rows goes through open_sorted_bed, which pipes them through sort in the same
order.  A file joined from parts that are each sorted is put in order afterwards
by sort_bed_file.
"""
import os
import pipettor
from flair.pycbio.sys import fileOps

# chrom and start, which is Bed.genome_sort_key, kept stable so that rows at one
# position stay in the order they were written
SORT_CMD = ('sort', '-s', '-k1,1', '-k2,2n')

# C collation and numeric parsing, so the order does not depend on the locale of
# whoever runs flair.  pipettor takes no environment, so this is set here, in the
# environment its children inherit.
os.environ['LC_COLLATE'] = 'C'
os.environ['LC_NUMERIC'] = 'C'

def bed_row_sort_key(row):
    "position order of a BED row of strings, which is what Bed.genome_sort_key gives a record"
    return row[0], int(row[1])

def _write_cmds(path):
    cmds = [SORT_CMD]
    if fileOps.isCompressed(path):
        cmds.append(fileOps.compressCmd(path))
    return cmds

def open_sorted_bed(path):
    """Writable file for BED rows, sorting them as they are written.  Compression
    follows the file name, as it does in opengz."""
    return pipettor.Popen(_write_cmds(path), mode='w', stdout=path)

def sort_bed_file(path):
    "Sort a BED that is already written, leaving it in the same place."
    cmds = _write_cmds(path)
    if fileOps.isCompressed(path):
        cmds.insert(0, fileOps.decompressCmd(path) + [path])
        stdin = None
    else:
        stdin = path
    with fileOps.AtomicFileCreate(path) as tmp_path:
        pipettor.run(cmds, stdin=stdin, stdout=tmp_path)
