"""The checks a subcommand makes before it starts work, and what it reports about
the result when it ends.  A mistake in the command line costs a second here rather
than the run, and a run that produced nothing reports that."""

import os
import shutil
import logging
import pysam
from flair import FlairInputDataError
from flair.io_utils import check_genome_fasta
from flair.run_logging import plural, plural_noun
from flair.pycbio.sys import fileOps

# the range a numeric option takes: its bounds, how the message states them, and what
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


# what a sample name may not contain, and what it would break: the name is written
# into a comma separated list, into tab separated tables, and into the names of the
# files of that sample
_SAMPLE_NAME_REJECTED = (
    (',', 'a comma, which separates the samples of an isoform BED and of a GTF attribute'),
    ('\t', 'a tab, which separates the fields of every table flair writes'),
    ('/', 'a slash, which would name a directory where a file is named for the sample'))


def check_sample_name(sample_name, where='--sample_name'):
    """A sample name that survives the files it is written into.  A name with a comma
    was read back as two samples a run later, by whatever read the list it is written
    into, and the error named the wrong cause."""
    for char, what in _SAMPLE_NAME_REJECTED:
        if char in sample_name:
            raise FlairInputDataError(f'{where} {sample_name!r} contains {what}')


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
                                  'the installing page of the FLAIR documentation states how to install: '
                                  'https://flair.readthedocs.io/en/latest/installing.html')


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


def check_bam_indexes(alignments):
    """every BAM a subcommand reads by region has its index.  Without it the fetch
    fails in a worker process, where the program's own message is all there is to go
    on, and after the setup that precedes the first fetch"""
    for path in alignments:
        with pysam.AlignmentFile(path) as bam:
            if not bam.has_index():
                raise FlairInputDataError(f'BAM {path} has no index, which is needed to read it by region; '
                                          f'make one with: samtools index {path}')


def _check_seq_names(what, path, names, genome, genome_names):
    """the sequences of an input are the genome's, a naming difference such as chr1
    against 1 being the common cause; sharing none of them is the error, so a BAM or
    BED covering part of the genome passes"""
    if names and genome_names.isdisjoint(names):
        raise FlairInputDataError(f'{what} {path} has none of the sequences of the genome {genome}: it names '
                                  f"{', '.join(sorted(names)[:3])} and the genome names "
                                  f"{', '.join(sorted(genome_names)[:3])}; the genome, the annotation, the "
                                  'alignments and the isoforms must name their sequences the same way')


def check_seq_names(genome, alignments=(), isoform_beds=(), vcfs=()):
    "the sequence names of the alignments, the isoform BEDs and the VCFs are the genome's"
    with pysam.FastaFile(genome) as genome_fa:
        genome_names = set(genome_fa.references)
    for what, paths, seq_names in (('alignments', alignments, _bam_seq_names),
                                   ('isoforms', isoform_beds, _bed_seq_names),
                                   ('variants', vcfs, _vcf_seq_names)):
        for path in paths:
            _check_seq_names(what, path, seq_names(path), genome, genome_names)


def check_annotation(gtf_data, genome, *, region_filtered=False):
    """A parsed annotation has transcripts, on the genome's sequences.  An annotation
    that matches nothing is otherwise the same as no annotation at all: the run
    finishes, with fewer isoforms and none of them named.  With region_filtered the
    parse kept only the transcripts of the regions asked for, so having none of them is
    not an error."""
    chroms = set(gtf_data.get_chroms())
    if len(chroms) == 0:
        if region_filtered:
            return
        raise FlairInputDataError(f'annotation GTF {gtf_data.gtf_file} has no transcripts; it must have '
                                  'transcript and exon lines, with GTF attributes, gene_id "GENE"; '
                                  'transcript_id "TRANSCRIPT";')
    with pysam.FastaFile(genome) as genome_fa:
        genome_names = set(genome_fa.references)
    _check_seq_names('annotation', gtf_data.gtf_file, chroms, genome, genome_names)


def precheck(*, files=None, programs=(), genome=None, alignments=(), indexed_alignments=(),
             isoform_beds=(), vcfs=()):
    """The checks a subcommand makes before it starts work: files, what a file is
    mapped to its path, exist; programs, the external programs it runs, are
    installed; genome, the FASTA, is compressed and indexed as the sequence reader
    needs; alignments, BAMs, indexed_alignments, the BAMs read by region, isoform_beds
    and vcfs name the genome's sequences; and each of indexed_alignments has its index.
    A None path is an option that was not given, and is skipped."""
    check_input_files(dict(files or {}, **({'genome FASTA': genome} if genome is not None else {})))
    check_programs(programs)
    check_bam_indexes(_given(indexed_alignments))
    if genome is not None:
        check_genome_fasta(genome)
        check_seq_names(genome, _given(alignments) + _given(indexed_alignments),
                        _given(isoform_beds), _given(vcfs))


def _given(paths):
    "the paths of an option that was given, the others being None"
    return [path for path in paths if path is not None]


def report_result_count(what, count, *, remedy, plural_what=None):
    """What a finished run produced, warning when it is nothing: an empty result is
    otherwise the same to look at as a run with nothing to report.  what is the
    singular noun of what was produced, plural_what its plural where adding s is
    wrong, and remedy states what to do about an empty result."""
    if count == 0:
        logging.warning('no %s; %s', plural_noun(0, what, plural_what), remedy)
    else:
        logging.info('%s', plural(count, what, plural_what))


def count_records(path):
    "the lines of a file a run wrote, for report_result_count"
    with fileOps.opengz(path) as fh:
        return sum(1 for line in fh if line.strip() != '')
