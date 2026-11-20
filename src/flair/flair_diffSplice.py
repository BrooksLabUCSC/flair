#! /usr/bin/env python3

import os
import os.path as osp
import sys
import pipettor
import argparse
import logging
from flair import FlairInputDataError, FlairToolError
from flair.io_utils import keep_or_remove, make_run_temp_dir
from flair.run_checks import (check_option_ranges, count_records, precheck,
                              report_result_count)
from flair.conditions import select_condition_pair
from flair.counts_matrix_tsv import read_isoform_ids, describe_ids
from flair.drimseq_input import check_feature_ids_unique
from flair.event_quant_tsv import EventQuantReader
from flair.sample_info_tsv import read_sample_info, write_sample_info
from flair.pycbio.hgdata.bed import BedReader
from flair.pycbio.sys import fileOps
from flair.pycbio.tsv import TsvWriter
from flair.output_compress import compressed_path, output_compression

pkgdir = osp.dirname(osp.realpath(__file__))
diffSplice_drimSeq = osp.join(pkgdir, "diffSplice_drimSeq.R")
call_diffsplice_events = osp.join(pkgdir, "call_diffsplice_events.py")
es_as = osp.join(pkgdir, "es_as.py")
es_as_inc_excl_to_counts = osp.join(pkgdir, "es_as_inc_excl_to_counts.py")

BED_NAME_COLUMN = 3
FEATURE_ID_COLUMN = 0

# the alternative splicing events called, in the order they are reported
EVENT_TYPES = ('es', 'alt5', 'alt3', 'ir')

# the programs a run executes, checked before it starts
PROGRAMS = ('Rscript',)


def add_subparser(subparsers):
    desc = "Call alternative splicing events from isoforms and test them for differential usage"
    parser = subparsers.add_parser('diffsplice', help="Call and test alternative splicing events",
                                   description=desc)
    required = parser.add_argument_group('required named arguments')
    required.add_argument('--isoform_bed', type=str, required=True,
                          help='isoforms in bed format')
    required.add_argument('--counts_matrix', type=str, required=True,
                          help='tab-delimited isoform count matrix from flair quantify')
    required.add_argument('-o', '--output', type=str, required=True,
                          help='output directory for tables and plots')
    parser.add_argument('-t', '--threads', type=int, default=12,
                        help='number of threads for parallel DRIMSeq (default: %(default)s)')
    parser.add_argument('--test', action='store_true',
                        help='run DRIMSeq statistical testing')
    parser.add_argument('--min_samps_gene_expr', type=int, default=6,
                        help='DRIMSeq dmFilter min_samps_gene_expr: minimum number of samples that have '
                             'coverage over an AS event inclusion/exclusion; events with too few samples '
                             'are filtered out and not tested (default: %(default)s)')
    parser.add_argument('--min_samps_feature_expr', type=int, default=3,
                        help='DRIMSeq dmFilter min_samps_feature_expr: minimum number of samples expressing '
                             'the inclusion of an AS event; events with too few samples are filtered out '
                             'and not tested (default: %(default)s)')
    parser.add_argument('--min_gene_expr', type=int, default=15,
                        help='DRIMSeq dmFilter min_gene_expr: minimum number of reads covering an AS event '
                             'inclusion/exclusion; events with too few reads are filtered out and not '
                             'tested (default: %(default)s)')
    parser.add_argument('--min_feature_expr', type=int, default=5,
                        help='DRIMSeq dmFilter min_feature_expr: minimum number of reads covering an AS '
                             'event inclusion; events with too few reads are filtered out and not tested '
                             '(default: %(default)s)')
    parser.add_argument('--batch', action='store_true',
                        help='implies --test. DRIMSeq corrects for the batch of each sample, '
                             'as the sample info TSV beside the counts matrix gives it')
    parser.add_argument('--condition_a', default='',
                        help='implies --test. The reference condition, as named in the counts matrix '
                             'columns; the comparison is --condition_b against this. With neither this '
                             'nor --condition_b given, the two conditions are taken in sorted order')
    parser.add_argument('--condition_b', default='',
                        help='the condition compared against --condition_a')
    parser.add_argument('--temp_dir',
                        help='directory for temporary files; each run makes its own directory, named for the output, in it (default: $TMPDIR or the system temporary directory)')
    parser.add_argument('--keep_intermediate', action='store_true',
                        help="keep the run's own directory, which holds the event tables DRIMSeq was given "
                             "and its diagnostics, and whose path the run reports")
    parser.add_argument('--overwrite_output', action='store_true', help=argparse.SUPPRESS)
    parser.set_defaults(entry=diffsplice_cmd)

def diffsplice_cmd(args):
    if args.overwrite_output:
        raise FlairInputDataError('--overwrite_output has been removed: -o names a file name base, as it does in every other subcommand, and the outputs are replaced as theirs are')
    if bool(args.condition_a) != bool(args.condition_b):
        # checked here as well as in select_condition_pair, which only runs when
        # testing is asked for; --condition_b alone would otherwise be ignored.
        # raise, not return 1: flair_cli discards the return value, so
        # flair diffsplice --condition_a X printed a line and exited 0 having
        # done no testing at all
        raise FlairInputDataError('--condition_a and --condition_b must both be given, '
                                  'or both left out to take the two conditions in sorted order')
    check_option_ranges(args, count=('min_samps_gene_expr', 'min_samps_feature_expr', 'min_gene_expr',
                                     'min_feature_expr'),
                        positive=('threads',))
    precheck(files={'isoform models BED': args.isoform_bed, 'counts matrix': args.counts_matrix},
             programs=PROGRAMS)
    diff_splice(isoform_bed=args.isoform_bed, counts_matrix=args.counts_matrix,
                output=args.output, threads=args.threads, test=args.test,
                min_samps_gene_expr=args.min_samps_gene_expr,
                min_samps_feature_expr=args.min_samps_feature_expr,
                min_gene_expr=args.min_gene_expr, min_feature_expr=args.min_feature_expr,
                batch=args.batch, condition_a=args.condition_a, condition_b=args.condition_b,
                temp_dir=args.temp_dir, keep_intermediate=args.keep_intermediate)

def check_isoform_ids(isoform_bed, counts_matrix):
    """Every counted isoform must have a BED record, since the events are called from
    the BED and the counts attached to them by id.  Extra BED records are allowed; they
    contribute no counts."""
    bed_names = frozenset(bed.name for bed in BedReader(isoform_bed, fixScores=True))
    missing = [iso_id for iso_id in read_isoform_ids(counts_matrix)
               if iso_id not in bed_names]
    if len(missing) > 0:
        raise FlairInputDataError(
            f"{len(missing)} isoform ids in counts matrix {counts_matrix} have no record in "
            f"{isoform_bed}: {describe_ids(missing)}; pass the isoform BED that was quantified "
            "to produce this counts matrix")

def check_batch_modelled(sample_infos, counts_matrix):
    """--batch adds the batch of each sample to the DRIMSeq model, which can be fitted
    only when the samples of the comparison are in more than one batch and the batches
    are not the conditions under other names.  Without this check R stops in
    model.matrix with "contrasts can be applied only to factors with 2 or more
    levels", which names neither the option nor the file to correct."""
    conditions = [si.condition for si in sample_infos]
    batches = [si.batch for si in sample_infos]
    combos = {(conditions.index(condition), batches.index(batch))
              for condition, batch in zip(conditions, batches)}
    if len(set(batches)) < 2:
        raise FlairInputDataError(f"--batch needs the samples of {conditions[0]} and {conditions[-1]} in more "
                                  f"than one batch; the sample info beside {counts_matrix} puts them all in "
                                  f"{batches[0]}.  Drop --batch, or give the batches in the flair quantify "
                                  "manifest and quantify again")
    if len(combos) == len(set(conditions)):
        raise FlairInputDataError(f"--batch needs a batch containing samples of more than one condition; in the "
                                  f"sample info beside {counts_matrix} each condition has a batch of its own, "
                                  "so the batch and the condition cannot be told apart.  Drop --batch")

def diff_splice(*, isoform_bed, counts_matrix, output, threads, test, min_samps_gene_expr,  # noqa: C901 - FIXME: reduce complexity
                min_samps_feature_expr, min_gene_expr, min_feature_expr, batch,
                condition_a, condition_b, temp_dir=None, keep_intermediate=False):
    if not os.path.exists(counts_matrix):
        raise FlairInputDataError(f'counts matrix file does not exist: {counts_matrix}')
    if not os.path.exists(isoform_bed):
        raise FlairInputDataError(f'isoform bed file does not exist: {isoform_bed}')
    check_isoform_ids(isoform_bed, counts_matrix)

    # the event tables DRIMSeq is given and its diagnostics, in a directory of this
    # run's own, removed at the end unless they were asked for
    workdir = make_run_temp_dir(output, temp_dir)

    with keep_or_remove(workdir, keep_intermediate):
        call_and_test_events(isoform_bed=isoform_bed, counts_matrix=counts_matrix, output=output,
                             workdir=workdir, threads=threads, test=test, batch=batch,
                             condition_a=condition_a, condition_b=condition_b,
                             min_samps_gene_expr=min_samps_gene_expr,
                             min_samps_feature_expr=min_samps_feature_expr,
                             min_gene_expr=min_gene_expr, min_feature_expr=min_feature_expr)


def report_event_counts(filebase):
    """The events of each type the run called, which are otherwise four tables to
    open."""
    counts = {event: max(count_records(compressed_path(f'{filebase}.{event}.events.quant.tsv')) - 1, 0)
              for event in EVENT_TYPES}
    logging.info('events called: %s', ', '.join(f'{event} {count}' for event, count in counts.items()))
    report_result_count('splicing event', sum(counts.values()),
                        remedy='check that the isoform BED is the one the counts matrix was made from, '
                        'and that its genes have more than one isoform')


def report_tested_events(output, condition_a, condition_b):
    "what DRIMSeq reported for each event type, which an empty table does not state"
    for event in EVENT_TYPES:
        path = drimseq_result_path(output, event, condition_a, condition_b)
        report_result_count(f'{event} event with differential usage', max(count_records(path) - 1, 0),
                            plural_what=f'{event} events with differential usage',
                            remedy='check that the conditions compared are the ones meant, or lower '
                            '--min_samps_gene_expr and --min_samps_feature_expr')


def call_and_test_events(*, isoform_bed, counts_matrix, output, workdir, threads, test, batch,
                         condition_a, condition_b, min_samps_gene_expr, min_samps_feature_expr,
                         min_gene_expr, min_feature_expr):
    "call the events of the isoforms and their counts, and test them when asked to"
    # the event tables are written beside the output prefix, as every other
    # subcommand writes its outputs
    filebase = output
    # the event callers are separate processes, so each is told whether to compress;
    # the exon skipping one writes to stdout, so this process opens that file
    call_args = [isoform_bed, filebase, counts_matrix]
    if output_compression():
        call_args.append('--compress')
    call_events(call_diffsplice_events, call_args)
    es_events_tsv = filebase + '.es.events.tsv'
    with open(es_events_tsv, 'w') as es_fh:
        call_events(es_as, [isoform_bed], stdout=es_fh)
    with fileOps.opengz(compressed_path(filebase + '.es.events.quant.tsv'), 'w') as quant_fh:
        call_events(es_as_inc_excl_to_counts, [counts_matrix, es_events_tsv], stdout=quant_fh)
    os.unlink(es_events_tsv)

    report_event_counts(filebase)

    # --condition_a and --batch imply --test, as each is only used in the testing
    # below: given without it, they were read and nothing was done with them
    if test or condition_a or batch:
        logging.info('DRIMSeq testing for each AS event type')
        # resolved here, not in the R script, so that diffexp and diffsplice choose
        # the same way and R is always told both names
        sample_infos = read_sample_info(counts_matrix)
        condition_a, condition_b = select_condition_pair([si.condition for si in sample_infos],
                                                         condition_a, condition_b, counts_matrix)
        if batch:
            check_batch_modelled([si for si in sample_infos if si.condition in (condition_a, condition_b)],
                                 counts_matrix)
        # the condition and batch of each sample, so that the R script does not have
        # to take them back out of the matrix column names
        formula_tsv = os.path.join(workdir, 'formula_matrix.tsv')
        write_sample_info(formula_tsv, sample_infos)
        ds_command = ['Rscript', diffSplice_drimSeq, '--threads', threads,
                      '--out_prefix', output, '--workdir', workdir,
                      '--condition_a', condition_a, '--condition_b', condition_b,
                      '--formula', formula_tsv,
                      '--min_samps_gene_expr', min_samps_gene_expr,
                      '--min_samps_feature_expr', min_samps_feature_expr,
                      '--min_gene_expr', min_gene_expr,
                      '--min_feature_expr', min_feature_expr]
        if output_compression():
            ds_command.append('--compress')
        if batch:
            ds_command += ['--batch']
        with open(workdir + '/ds.stderr.txt', 'w') as ds_stderr:
            for event in EVENT_TYPES:
                run_drimseq_event(ds_command, event, filebase, workdir, ds_stderr,
                                  output, condition_a, condition_b, sample_infos)
        report_tested_events(output, condition_a, condition_b)
    else:
        logging.info('no statistical test was run: give --condition_a and --condition_b, or --test, '
                     'to test the events with DRIMSeq')

def call_events(program, args, stdout=None):
    """Run one of the event callers.  Its own message is already on stderr, so an
    unhandled ProcessException would only add a traceback of the flair code that
    started it."""
    try:
        pipettor.run([sys.executable, program] + args, stdout=stdout)
    except pipettor.ProcessException as exc:
        raise FlairToolError(f"{osp.basename(program)} failed, see the error above") from exc

def drimseq_result_path(output, event, condition_a, condition_b):
    "the result file for one event type, as diffSplice_drimSeq.R names it"
    return compressed_path(f'{output}.drimseq_{event}_{condition_a}_v_{condition_b}.tsv')

def write_empty_drimseq_result(output, event, condition_a, condition_b, sample_infos):
    """The file DRIMSeq would have written had the event type called anything, so that
    a run always writes the four results it declares.  The columns are the ones
    diffSplice_drimSeq.R ends with: the two ids, the samples of this comparison in
    formula order, then lr and adj_pvalue."""
    samples = [si.sample_id for si in sample_infos
               if si.condition in (condition_a, condition_b)]
    columns = ['feature_id', 'gene_id'] + samples + ['lr', 'adj_pvalue']
    with TsvWriter(drimseq_result_path(output, event, condition_a, condition_b),
                   columns=columns, defaultColType=str):
        pass

def run_drimseq_event(ds_command, event, filebase, workdir, ds_stderr,
                      output, condition_a, condition_b, sample_infos):
    matrixfile = compressed_path(f'{filebase}.{event}.events.quant.tsv')
    if empty_matrix(matrixfile):
        logging.info(f'{event} event matrix file empty, not running DRIMSeq')
        write_empty_drimseq_result(output, event, condition_a, condition_b, sample_infos)
    else:
        check_feature_ids_unique(matrixfile)
        cur_command = ds_command + ['--matrix', matrixfile, '--prefix', event]
        try:
            pipettor.run(cur_command, stderr=ds_stderr)
        except pipettor.ProcessException as exc:
            raise FlairToolError(f"DRIMSeq failed on `{event}' event. "
                                 f'Check {workdir}/ds.stderr.txt for details') from exc

def empty_matrix(event_quant_tsv):
    "true when an event type called no events, leaving a file of only its header"
    with EventQuantReader(event_quant_tsv) as reader:
        return next(iter(reader), None) is None
