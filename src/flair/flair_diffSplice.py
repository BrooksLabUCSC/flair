#! /usr/bin/env python3

import os
import os.path as osp
import sys
import pipettor
import logging
from flair import FlairError, FlairInputDataError
from flair.conditions import select_condition_pair
from flair.counts_matrix_tsv import read_isoform_ids, describe_ids
from flair.event_quant_tsv import EventQuantReader
from flair.sample_info_tsv import read_sample_info, write_sample_info
from flair.pycbio.hgdata.bed import BedReader
from flair.pycbio.sys import fileOps
from flair.output_compress import compressed_path, output_compression

pkgdir = osp.dirname(osp.realpath(__file__))
diffSplice_drimSeq = osp.join(pkgdir, "diffSplice_drimSeq.R")
call_diffsplice_events = osp.join(pkgdir, "call_diffsplice_events.py")
es_as = osp.join(pkgdir, "es_as.py")
es_as_inc_excl_to_counts = osp.join(pkgdir, "es_as_inc_excl_to_counts.py")

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
    parser.add_argument('-t', '--threads', type=int, default=4,
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
                        help='with --test, DRIMSeq will perform batch correction')
    parser.add_argument('--condition_a', default='',
                        help='implies --test. The reference condition, as named in the counts matrix '
                             'columns; the comparison is --condition_b against this. With neither this '
                             'nor --condition_b given, the two conditions are taken in sorted order')
    parser.add_argument('--condition_b', default='',
                        help='the condition compared against --condition_a')
    parser.add_argument('--overwrite_output', action='store_true',
                        help='overwrite files in an existing output directory')
    parser.set_defaults(entry=diffsplice_cmd)

def diffsplice_cmd(args):
    if bool(args.condition_a) != bool(args.condition_b):
        # checked here as well as in select_condition_pair, which only runs when
        # testing is asked for; --condition_b alone would otherwise be ignored.
        # raise, not return 1: flair_cli discards the return value, so
        # flair diffsplice --condition_a X printed a line and exited 0 having
        # done no testing at all
        raise FlairInputDataError('--condition_a and --condition_b must both be given, '
                                  'or both left out to take the two conditions in sorted order')
    diff_splice(isoform_bed=args.isoform_bed, counts_matrix=args.counts_matrix,
                output=args.output, threads=args.threads, test=args.test,
                min_samps_gene_expr=args.min_samps_gene_expr,
                min_samps_feature_expr=args.min_samps_feature_expr,
                min_gene_expr=args.min_gene_expr, min_feature_expr=args.min_feature_expr,
                batch=args.batch, condition_a=args.condition_a, condition_b=args.condition_b,
                overwrite_output=args.overwrite_output)

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

def diff_splice(*, isoform_bed, counts_matrix, output, threads, test, min_samps_gene_expr,  # noqa: C901 - FIXME: reduce complexity
                min_samps_feature_expr, min_gene_expr, min_feature_expr, batch,
                condition_a, condition_b, overwrite_output):
    if not os.path.exists(counts_matrix):
        raise FlairInputDataError('Counts matrix file path does not exist')
    if not os.path.exists(isoform_bed):
        raise FlairInputDataError('Isoform bed file path does not exist')
    check_isoform_ids(isoform_bed, counts_matrix)

    # Create output directory including a working directory for intermediate files.
    workdir = os.path.join(output, 'workdir')
    if overwrite_output:
        if not os.path.exists(workdir):
            os.makedirs(workdir)
        pass
    elif not os.path.exists(output):
        try:
            os.makedirs(workdir, 0o700)
        except OSError as ex:
            raise OSError("** ERROR cannot create directory %s" % (workdir)) from ex
    else:
        raise FlairInputDataError(f'** Error. Name {output} already exists. Choose another name for out_dir')

    filebase = os.path.join(output, 'diffsplice')
    # the event callers are separate processes, so each is told whether to compress;
    # the exon skipping one writes to stdout, so this process opens that file
    call_cmd = [sys.executable, call_diffsplice_events, isoform_bed, filebase, counts_matrix]
    if output_compression():
        call_cmd.append('--compress')
    pipettor.run(call_cmd)
    es_events_tsv = filebase + '.es.events.tsv'
    with open(es_events_tsv, 'w') as es_fh:
        pipettor.run([sys.executable, es_as, isoform_bed], stdout=es_fh)
    with fileOps.opengz(compressed_path(filebase + '.es.events.quant.tsv'), 'w') as quant_fh:
        pipettor.run([sys.executable, es_as_inc_excl_to_counts, counts_matrix, es_events_tsv],
                     stdout=quant_fh)
    os.unlink(es_events_tsv)

    if test or condition_a:
        logging.info('DRIMSeq testing for each AS event type')
        # resolved here, not in the R script, so that diffexp and diffsplice choose
        # the same way and R is always told both names
        sample_infos = read_sample_info(counts_matrix)
        condition_a, condition_b = select_condition_pair([si.condition for si in sample_infos],
                                                         condition_a, condition_b, counts_matrix)
        # the condition and batch of each sample, so that the R script does not have
        # to take them back out of the matrix column names
        formula_tsv = os.path.join(workdir, 'formula_matrix.tsv')
        write_sample_info(formula_tsv, sample_infos)
        ds_command = ['Rscript', diffSplice_drimSeq, '--threads', threads, '--out_dir', output,
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
            for event in ('es', 'alt5', 'alt3', 'ir'):
                run_drimseq_event(ds_command, event, filebase, workdir, ds_stderr)
    else:
        # workdir only necessary for drimseq output
        os.rmdir(workdir)

def run_drimseq_event(ds_command, event, filebase, workdir, ds_stderr):
    matrixfile = compressed_path(f'{filebase}.{event}.events.quant.tsv')
    if empty_matrix(matrixfile):
        logging.info(f'{event} event matrix file empty, not running DRIMSeq')
    else:
        cur_command = ds_command + ['--matrix', matrixfile, '--prefix', event]
        try:
            pipettor.run(cur_command, stderr=ds_stderr)
        except pipettor.ProcessException as exc:
            raise FlairError(f"DRIMSeq failed on `{event}' event. "
                             f'Check {workdir}/ds.stderr.txt for details') from exc

def empty_matrix(event_quant_tsv):
    "true when an event type called no events, leaving a file of only its header"
    with EventQuantReader(event_quant_tsv) as reader:
        return next(iter(reader), None) is None
