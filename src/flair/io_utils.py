"""Filesystem / I/O helpers shared across FLAIR pipelines."""

import os
import tempfile
import pipettor
import pysam
from flair import FlairInputDataError
from flair.output_compress import compressed_path
from flair.pycbio.sys import fileOps

# chrom and start, which is Bed.genome_sort_key, kept stable so that rows at one
# position stay in the order they were written
BED_SORT_CMD = ('sort', '-s', '-k1,1', '-k2,2n')

# C collation and numeric parsing, so the order does not depend on the locale of
# whoever runs flair.  pipettor takes no environment, so this is set here, in the
# environment its children inherit.
os.environ['LC_COLLATE'] = 'C'
os.environ['LC_NUMERIC'] = 'C'


def make_temp_dir(out_prefix):
    """The intermediate directory of a run, kept if it is already there, so a stage
    can reuse what an earlier run left.  A stage whose results depend on what the
    directory holds wants make_run_temp_dir instead, a new directory each run."""
    # FIXME: use TMPDIR unless directory explicitly specified
    temp_dir = out_prefix + ".intermediate"
    try:
        os.makedirs(temp_dir, exist_ok=True)
    except OSError as exc:
        raise OSError(f"Creation of the directory `{temp_dir}' failed") from exc
    return temp_dir + '/'


def make_run_temp_dir(out_prefix, temp_dir=None):
    """Create a new directory for one run's temporary files, named for the output
    prefix, in temp_dir, or in the system temporary directory ($TMPDIR) if it is
    None.  Each run gets its own directory, so runs can share temp_dir."""
    base = temp_dir if temp_dir is not None else tempfile.gettempdir()
    try:
        os.makedirs(base, exist_ok=True)
        run_dir = tempfile.mkdtemp(prefix=os.path.basename(out_prefix) + '.', suffix='.intermediate', dir=base)
    except OSError as exc:
        raise OSError(f"Creation of a temporary directory in `{base}' failed") from exc
    return run_dir + '/'


def read_manifest_rows(path, field_counts, line_format, file_fields):
    """A manifest's rows, each a list of its tab-separated fields, for a line each but
    blank lines and '#' comments.  A line has one of field_counts fields, none of them
    empty, as line_format describes for errors; the fields at the indexes of
    file_fields, index -> what the file is, name files that must exist; and the first
    field, the sample's name, is unique."""
    rows = []
    # plain or compressed, as every flair text input may be
    with fileOps.opengz(path) as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.rstrip('\r\n')
            if line.strip() == '' or line.startswith('#'):
                continue
            fields = line.split('\t')
            where = f'{path}:{lineno}'
            if len(fields) not in field_counts or '' in fields:
                spaces = '; fields are separated by tabs, not spaces' if len(fields) == 1 and ' ' in line.strip() else ''
                raise FlairInputDataError(f'{where}: a manifest line is {line_format}, tab-separated{spaces}: {line!r}')
            for index, what in file_fields.items():
                if index < len(fields) and not os.path.exists(fields[index]):
                    raise FlairInputDataError(f'{where}: {what} {fields[index]} does not exist')
            rows.append(fields)
    if len(rows) == 0:
        raise FlairInputDataError(f'{path}: manifest has no samples')
    names = [row[0] for row in rows]
    duplicated = sorted({name for name in names if names.count(name) > 1})
    if duplicated:
        raise FlairInputDataError(f'{path}: samples named more than once: {", ".join(duplicated)}')
    return rows


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
    if not fileOps.isCompressed(genome):
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


##
# Final BED output, sorted by position.
#
# A BED that flair hands to the user is sorted by chrom and start, so that a run is
# comparable with another and so that bedToBigBed and genome browsers get the order
# they expect.
#
# Rows that are all in memory are sorted with Bed.genome_sort_key.  A stream of rows
# goes through open_sorted_bed, which pipes them through sort in the same order.  A
# file joined from parts that are each sorted is put in order afterwards by
# sort_bed_file.
##

def bed_row_sort_key(row):
    "position order of a BED row of strings, which is what Bed.genome_sort_key gives a record"
    return row[0], int(row[1])

def _bed_sort_cmds(path):
    cmds = [BED_SORT_CMD]
    if fileOps.isCompressed(path):
        cmds.append(fileOps.compressCmd(path))
    return cmds

def open_sorted_bed(path):
    """Writable file for BED rows, sorting them as they are written.  Compression
    follows the file name, as it does in opengz."""
    return pipettor.Popen(_bed_sort_cmds(path), mode='w', stdout=path)

def sort_bed_file(path):
    "Sort a BED that is already written, leaving it in the same place."
    cmds = _bed_sort_cmds(path)
    if fileOps.isCompressed(path):
        cmds.insert(0, fileOps.decompressCmd(path) + [path])
        stdin = None
    else:
        stdin = path
    with fileOps.AtomicFileCreate(path) as tmp_path:
        pipettor.run(cmds, stdin=stdin, stdout=tmp_path)


def write_sorted_vcf(path, header, records):
    """Write VCF records in position order, under the name compressed_path gives, with
    a tabix index beside a compressed one.  The callers build their records per variant
    rather than per position, and tabix needs the file sorted, so it is sorted here."""
    # imported here rather than at module scope, so that a FLAIR install missing
    # this package still runs every other subcommand
    import vcfpy
    out_path = compressed_path(path)
    with vcfpy.Writer.from_path(out_path, header) as writer:
        for record in sorted(records, key=lambda rec: (rec.CHROM, rec.POS)):
            writer.write_record(record)
    if fileOps.isCompressed(out_path):
        pysam.tabix_index(out_path, preset='vcf', force=True)
    return out_path
