"""Filesystem / I/O helpers shared across FLAIR pipelines."""

import os
import shutil
from flair import FlairInputDataError
from flair.pycbio.sys.fileOps import isCompressed


def make_temp_dir(out_prefix):
    """The intermediate directory of a run, kept if it is already there, so a stage
    can reuse what an earlier run left.  A stage whose results depend on what the
    directory holds wants make_clean_temp_dir instead."""
    # FIXME: use TMPDIR unless directory explicitly specified
    temp_dir = out_prefix + ".intermediate"
    try:
        os.makedirs(temp_dir, exist_ok=True)
    except OSError as exc:
        raise OSError(f"Creation of the directory `{temp_dir}' failed") from exc
    return temp_dir + '/'


def make_clean_temp_dir(out_prefix):
    """The intermediate directory of a run, emptied first.  For a stage that reads
    back every file it finds there: a file left by an earlier run, which a crash or
    --keep_intermediate can leave behind, would otherwise be counted as this run's."""
    temp_dir = out_prefix + ".intermediate"
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)
    return make_temp_dir(out_prefix)


# a bgzip file is gzip with an extra field whose two id bytes are BC
GZIP_MAGIC = b'\x1f\x8b'
GZIP_FEXTRA = 0x04
BGZF_EXTRA_ID = b'BC'

def is_bgzip(fasta):
    "true when a compressed file is bgzip rather than plain gzip"
    with open(fasta, 'rb') as fh:
        header = fh.read(14)
    return ((len(header) == 14) and (header[:2] == GZIP_MAGIC)
            and bool(header[3] & GZIP_FEXTRA) and (header[12:14] == BGZF_EXTRA_ID))

def check_genome_fasta(genome):
    """Check a genome FASTA before a run starts, rather than after the setup that
    precedes the first sequence fetch.  A compressed genome must be bgzip, since the
    sequence reader cannot index plain gzip, and must carry its index, since the
    reader can only build one where it may write beside the genome."""
    if not os.path.exists(genome):
        raise FlairInputDataError(f'genome FASTA does not exist: {genome}')
    if not isCompressed(genome):
        return
    if not (genome.endswith('.gz') or genome.endswith('.bgz')):
        raise FlairInputDataError(
            f'genome FASTA must be plain or bgzip compressed, not {genome}; '
            'recompress it with bgzip and index it with samtools faidx')
    if not is_bgzip(genome):
        raise FlairInputDataError(
            f'genome FASTA {genome} is plain gzip, which cannot be indexed for random '
            'access; recompress it with bgzip and index it with samtools faidx')
    missing = [genome + suffix for suffix in ('.fai', '.gzi')
               if not os.path.exists(genome + suffix)]
    if len(missing) > 0:
        raise FlairInputDataError(
            f"genome FASTA {genome} is missing {', '.join(missing)}; "
            f'index it with: samtools faidx {genome}')
