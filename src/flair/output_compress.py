"""Whether this run writes its text outputs compressed.

A flair run compresses all of its text outputs or none of them, so this is one
setting for the run rather than an argument threaded through every function that
writes a file.  flair --compress sets it once, before any subcommand starts.

What it covers is narrow on purpose:

- the text files flair itself writes, which gain a .gz suffix
- not BAM, which is compressed already, nor bigBed, nor the PDFs the R scripts
  write, none of which flair can gzip
- not the intermediates under <output>.intermediate, which are joined and
  measured as plain bytes, and are deleted at the end of the run anyway

Reading needs no option at all: every text input flair reads may be compressed,
whether or not this run writes compressed output.
"""

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
