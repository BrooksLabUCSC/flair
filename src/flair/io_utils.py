"""The run temporary directory, manifest reading, genome FASTA checks, and
position sorted BED and VCF output."""

import logging
import os
import shutil
import tempfile
from contextlib import contextmanager
import pipettor
import pysam
from flair import FlairInputDataError, FlairToolError
from flair.output_compress import compressed_path
from flair.run_logging import debug_logging
from flair.pycbio.sys import fileOps

# chrom and start, which is Bed.genome_sort_key, kept stable so that rows at one
# position stay in the order they were written
BED_SORT_CMD = ('sort', '-s', '-k1,1', '-k2,2n')

# C collation and numeric parsing, so the order does not depend on the locale of
# whoever runs flair.  pipettor takes no environment, so this is set here, in the
# environment its children inherit.  LC_ALL, which POSIX gives precedence over
# every other LC_ variable, is set rather than the two it would override.
os.environ['LC_ALL'] = 'C'


def make_temp_dir(out_prefix):
    """The intermediate directory of a run, kept if it is already there, so a stage
    can reuse what an earlier run left.  A stage whose results depend on what the
    directory contains requires make_run_temp_dir instead, a new directory each run."""
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


@contextmanager
def keep_or_remove(run_dir, keep):
    """Remove the run's temporary directory at the end, or report where it is when
    it was asked for.  It goes on a failure too, so a run leaves nothing of its
    own."""
    try:
        yield run_dir
    finally:
        if keep:
            logging.info('intermediate files kept in %s', run_dir)
        else:
            shutil.rmtree(run_dir, ignore_errors=True)


@contextmanager
def remove_on_failure(path, remove):
    """Remove a directory this run created when the run fails before writing a
    result, so that the rerun after a failure is not rejected for the directory
    being there.  remove is False for a directory the run was given rather than
    made, whose contents are not this run's to delete."""
    try:
        yield
    except Exception:
        if remove:
            shutil.rmtree(path, ignore_errors=True)
        raise


def _parse_manifest_line(where, line, field_counts, line_format, file_fields):
    "the fields of one manifest line, checked"
    fields = line.split('\t')
    if len(fields) not in field_counts or '' in fields:
        spaces = '; fields are separated by tabs, not spaces' if len(fields) == 1 and ' ' in line.strip() else ''
        raise FlairInputDataError(f'{where}: a manifest line is {line_format}, tab-separated{spaces}: {line!r}')
    for index, what in file_fields.items():
        if index < len(fields) and not os.path.exists(fields[index]):
            raise FlairInputDataError(f'{where}: {what} {fields[index]} does not exist')
    return fields


def read_manifest_rows(path, field_counts, line_format, file_fields, *, named='samples'):
    """A manifest's rows, each a list of its tab-separated fields, for a line each but
    blank lines and '#' comments.  A line has one of field_counts fields, none of them
    empty, as line_format describes for errors; the fields at the indexes of
    file_fields, index -> what the file is, name files that must exist; and the first
    field is unique.  named is what that field names, for the messages about an empty
    manifest and a repeated one."""
    rows = []
    # plain or compressed, as every flair text input may be
    with fileOps.opengz(path) as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.rstrip('\r\n')
            if (line.strip() != '') and not line.startswith('#'):
                rows.append(_parse_manifest_line(f'{path}:{lineno}', line, field_counts, line_format, file_fields))
    if len(rows) == 0:
        raise FlairInputDataError(f'{path}: manifest has no {named}')
    names = [row[0] for row in rows]
    duplicated = sorted({name for name in names if names.count(name) > 1})
    if duplicated:
        raise FlairInputDataError(f'{path}: {named} named more than once: {", ".join(duplicated)}')
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
    sequence reader cannot index plain gzip, and must have its index, since the
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


def _program_names(cmds):
    "the programs of a pipettor command or pipeline, for a message"
    commands = [cmds] if isinstance(cmds[0], str) else cmds
    return ' | '.join(os.path.basename(command[0]) for command in commands)


def _program_failed_msg(cmds, exc):
    """A failed program's message, with its own diagnostics when pipettor captured
    them.  The cause does not survive a process pool worker, so the diagnostics must
    be in the message itself to be reported at all."""
    diagnostics = (exc.stderr or '').strip()
    where = f'\n{diagnostics}' if diagnostics != '' else ', its diagnostics are above'
    return f'{_program_names(cmds)} failed{where}'


def run_program(cmds, **kwargs):
    """Run a program, or a pipeline of them, reporting a failure as a tool error with
    the program's own diagnostics, so that a traceback of the flair code that started
    it adds nothing.  With debug logging their stderr is inherited rather than
    captured, so that what a flair program run this way logs is printed as it runs."""
    if debug_logging() and ('stderr' not in kwargs):
        kwargs['stderr'] = None
    try:
        pipettor.run(cmds, **kwargs)
    except pipettor.ProcessException as exc:
        raise FlairToolError(_program_failed_msg(cmds, exc)) from exc


class GenomeFasta(pysam.FastaFile):
    """The genome, reporting a request for a sequence it does not have as an input
    error.  The reader signals a KeyError, which reports nothing about where the name
    came from, and the cause is almost always the genome being named differently
    from the alignments and the isoforms, chr1 against 1."""

    def fetch(self, *args, **kwargs):
        try:
            return super().fetch(*args, **kwargs)
        except KeyError as exc:
            raise FlairInputDataError(self._missing_seq_msg(kwargs.get('reference', args[0] if args else None))) from exc

    def get_reference_length(self, reference):
        try:
            return super().get_reference_length(reference)
        except KeyError as exc:
            raise FlairInputDataError(self._missing_seq_msg(reference)) from exc

    def _missing_seq_msg(self, seq_name):
        return (f"sequence `{seq_name}' is not in the genome FASTA {os.fsdecode(self.filename)}; the genome, "
                'the alignments and the isoforms must name their sequences the same way')


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

@contextmanager
def open_sorted_bed(path):
    """Writable file for BED rows, sorting them as they are written.  Compression
    follows the file name, as it does in opengz.  The rows are written under a
    temporary name in the same directory and installed under path when they are all
    written, as sort_bed_file and the TSV writers do, so a failed run leaves no part
    of a BED under the name of a whole one, and the directory of path is created if
    it does not exist."""
    with fileOps.AtomicFileCreate(path) as tmp_path:
        with pipettor.Popen(_bed_sort_cmds(path), mode='w', stdout=tmp_path) as fh:
            yield fh

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
    # in the order the header declares the contigs, which tools reading the file
    # check against it, rather than by name, which puts chr10 before chr2
    order = {line.id: num for num, line in enumerate(header.get_lines('contig'))}
    with vcfpy.Writer.from_path(out_path, header) as writer:
        for record in sorted(records, key=lambda rec: (order.get(rec.CHROM, len(order)), rec.CHROM, rec.POS)):
            writer.write_record(record)
    if fileOps.isCompressed(out_path):
        pysam.tabix_index(out_path, preset='vcf', force=True)
    return out_path
