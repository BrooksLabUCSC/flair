"""The checks a subcommand makes before it starts work, and what it reports about
the result when it ends.  A mistake in the command line costs a second here rather
than the run, and a run that produced nothing says so."""

import os
import shutil
import logging
import pysam
from flair import FlairInputDataError
from flair.io_utils import check_genome_fasta
from flair.pycbio.sys import fileOps

# the range a numeric option takes: its bounds, how the message says them, and what
# to do about a value outside them
_OPTION_RANGES = {
    'fraction': (0.0, 1.0, 'a fraction of 0 to 1', 'a percentage goes in as its fraction, 5 percent as 0.05'),
    'percent': (0, 100, 'a percentage of 0 to 100', None),
    'count': (0, None, '0 or more', None),
    'positive': (1, None, '1 or more', None),
}


def check_option_ranges(args, **kinds):
    """The numeric options of a subcommand, in range, as fraction=, percent=, count=
    or positive=, each a tuple of argparse dests.  A value outside the range makes
    the run report nothing or fail somewhere else, so it is an error here."""
    for kind, dests in kinds.items():
        _check_option_kind(args, kind, dests)


def _check_option_kind(args, kind, dests):
    low, high, what, hint = _OPTION_RANGES[kind]
    for dest in dests:
        value = getattr(args, dest)
        if (value is not None) and ((value < low) or ((high is not None) and (value > high))):
            raise FlairInputDataError(f"--{dest} is {value}, which is not {what}"
                                      + (f'; {hint}' if hint is not None else ''))


def check_input_files(files):
    """every input file of a run exists; files maps what a file is, for the message,
    to its path, or to several paths of that kind.  A None path is an option that was
    not given"""
    for what, paths in files.items():
        for path in ((paths,) if isinstance(paths, str) else (paths or ())):
            if not os.path.exists(path):
                raise FlairInputDataError(f'{what} file does not exist: {path}')


def check_programs(programs):
    """every program a subcommand runs is installed, rather than the run stopping at
    the first one it executes, which may be after the alignment"""
    missing = [program for program in programs if shutil.which(program) is None]
    if missing:
        raise FlairInputDataError(f"not on PATH: {', '.join(missing)}; flair runs these programs, which "
                                  'docs/source/installing.rst says how to install')


def _bed_seq_names(bed_file):
    "the sequence names of a BED"
    with fileOps.opengz(bed_file) as fh:
        return {line.split('\t', 1)[0] for line in fh if line.strip() != ''}


def _bam_seq_names(bam_file):
    "the sequence names of a BAM or CRAM header"
    with pysam.AlignmentFile(bam_file) as bam:
        return set(bam.references)


def _vcf_seq_names(vcf_file):
    """the sequence names a VCF header declares, none when it declares no contigs,
    which the format does not require"""
    with pysam.VariantFile(vcf_file) as vcf:
        return set(vcf.header.contigs)


def _check_seq_names(what, path, names, genome, genome_names):
    """the sequences of an input are the genome's, a naming difference such as chr1
    against 1 being the common cause; sharing none of them is the error, so a BAM or
    BED covering part of the genome passes"""
    if names and genome_names.isdisjoint(names):
        raise FlairInputDataError(f'{what} {path} has none of the sequences of the genome {genome}: it names '
                                  f"{', '.join(sorted(names)[:3])} and the genome names "
                                  f"{', '.join(sorted(genome_names)[:3])}; the alignments, the isoforms and "
                                  'the genome must name their sequences the same way')


def check_seq_names(genome, alignments=(), isoform_beds=(), vcfs=()):
    "the sequence names of the alignments, the isoform BEDs and the VCFs are the genome's"
    with pysam.FastaFile(genome) as genome_fa:
        genome_names = set(genome_fa.references)
    for what, paths, seq_names in (('alignments', alignments, _bam_seq_names),
                                   ('isoforms', isoform_beds, _bed_seq_names),
                                   ('variants', vcfs, _vcf_seq_names)):
        for path in paths:
            _check_seq_names(what, path, seq_names(path), genome, genome_names)


def precheck(*, files=None, programs=(), genome=None, alignments=(), isoform_beds=(), vcfs=()):
    """The checks a subcommand makes before it starts work: files, what a file is
    mapped to its path, exist; programs, the external programs it runs, are
    installed; genome, the FASTA, is compressed and indexed as the sequence reader
    needs; and alignments, BAMs, isoform_beds and vcfs name the genome's sequences.
    A None path is an option that was not given, and is skipped."""
    check_input_files(dict(files or {}, **({'genome FASTA': genome} if genome is not None else {})))
    check_programs(programs)
    if genome is not None:
        check_genome_fasta(genome)
        check_seq_names(genome, _given(alignments), _given(isoform_beds), _given(vcfs))


def _given(paths):
    "the paths of an option that was given, the others being None"
    return [path for path in paths if path is not None]


def report_result_count(what, count, *, remedy):
    """What a finished run produced, warning when it is nothing: an empty result is
    otherwise the same to look at as a run with nothing to report.  remedy says what
    to do about an empty one."""
    if count == 0:
        logging.warning('no %s; %s', what, remedy)
    else:
        logging.info('%d %s', count, what)


def count_records(path):
    "the lines of a file a run wrote, for report_result_count"
    with fileOps.opengz(path) as fh:
        return sum(1 for line in fh if line.strip() != '')
