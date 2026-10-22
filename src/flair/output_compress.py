"""Whether this run writes its text outputs compressed.

A flair run compresses all of its text outputs or none of them, so this is one
setting for the run rather than an argument threaded through every function that
writes a file.  flair --compress sets it once, before any subcommand starts.

What it covers:

- the text files flair itself writes, which gain a .gz suffix: the tables, the
  isoform BED, the GTF and the FASTA
- FASTA, and flair quantify's read maps, are bgzip rather than plain gzip, so that
  samtools faidx and pysam can index what flair wrote, and flair isoalleles read a
  read map's rows by offset
- not BAM, which is compressed already, nor bigBed, nor the PDFs the R scripts
  write, none of which flair can gzip
- not the intermediates under <output>.intermediate, which are joined and
  measured as plain bytes, and are deleted at the end of the run anyway

Reading needs no option at all: every text input flair reads may be compressed,
whether or not this run writes compressed output.
"""
import os.path as osp
from flair import FlairInputDataError
from flair.pycbio.sys import fileOps

COMPRESS_SUFFIX = '.gz'

# set once from the command line, read wherever a final output path is built
_compress_output = False

def set_output_compression(compress):
    "called once, by the flair command line, before any subcommand runs"
    global _compress_output
    _compress_output = bool(compress)

def output_compression():
    "whether this run writes its text outputs compressed"
    return _compress_output

def compressed_path(path):
    """The name to write a text output under: the path itself, or the path plus
    the compression suffix when the run asked for compression.  Pass only final
    outputs; intermediates stay plain."""
    return path + COMPRESS_SUFFIX if _compress_output else path

def open_output(path):
    "writable text file under the name compressed_path gives it"
    return fileOps.opengz(compressed_path(path), 'w')

def open_seekable_output(path):
    """Writable text file under the name compressed_path gives it, bgzip compressed
    rather than gzip when the run asked for compression, so it can be read from an
    offset: a FASTA by samtools faidx, a read map by read_map_index."""
    return fileOps.opengz(compressed_path(path), 'w', bgzip=True)

def open_fasta_output(path):
    "Writable FASTA, bgzip compressed when the run asked for compression (open_seekable_output)"
    return open_seekable_output(path)

def resolved_input_path(path):
    """The file a derived input path names, when flair builds the name rather than
    being given it: the path itself, or its compressed form when only that exists.
    A prefix from a compressed run has .gz on every text file it wrote."""
    compressed = path + COMPRESS_SUFFIX
    have, have_compressed = osp.exists(path), osp.exists(compressed)
    if have and have_compressed:
        raise FlairInputDataError(f"both {path} and {compressed} exist; remove the one that is "
                                  "stale, or name the file to use rather than its prefix")
    return compressed if have_compressed else path
