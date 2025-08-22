#! /usr/bin/env python3

import sys
import argparse
import os
import os.path as osp
import pipettor
import logging
from flair import FlairError, set_unix_path, FlairInputDataError

pkgdir = osp.dirname(osp.realpath(__file__))
diffSplice_drimSeq = osp.join(pkgdir, "diffSplice_drimSeq.R")

# FIXME: restructure, odd the say argument parsing is done bases on function
# arguments

def diffSplice(isoforms='', counts_matrix=''):
    set_unix_path()
    parser = argparse.ArgumentParser()
    required = parser.add_argument_group('required named arguments')
    if not isoforms:
        required.add_argument('-i', '--isoforms', action='store', required=True,
                type=str, help='isoforms in bed format')
        required.add_argument('-q', '--counts_matrix', action='store',
                type=str, required=True, help='tab-delimited isoform count matrix from flair quantify module')
    required.add_argument('-o', '--out_dir', action='store',
            type=str, required=True, help='Output directory for tables and plots.')
    parser.add_argument('-t', '--threads', action='store',
            type=int, required=False, default=4, help='Number of threads for parallel DRIMSeq (4)')
    parser.add_argument('--test', action='store_true', dest='test',
            required=False, default=False, help='Run DRIMSeq statistical testing')
    parser.add_argument('--drim1', action='store', dest='drim1', type=int, required=False, default=6,
            help='''The minimum number of samples that have coverage over an AS event inclusion/exclusion
            for DRIMSeq testing; events with too few samples are filtered out and not tested (6)''')
    parser.add_argument('--drim2', action='store', dest='drim2', type=int, required=False, default=3,
            help='''The minimum number of samples expressing the inclusion of an AS event;
             events with too few samples are filtered out and not tested (3)''')
    parser.add_argument('--drim3', action='store', dest='drim3', type=int, required=False, default=15,
            help='''The minimum number of reads covering an AS event inclusion/exclusion for DRIMSeq testing,
             events with too few samples are filtered out and not tested (15)''')
    parser.add_argument('--drim4', action='store', dest='drim4', type=int, required=False, default=5,
            help='''The minimum number of reads covering an AS event inclusion for DRIMSeq testing,
             events with too few samples are filtered out and not tested (5)''')
    parser.add_argument('--batch', action='store_true', dest='batch', required=False, default=False,
            help='''If specified with --test, DRIMSeq will perform batch correction''')
    parser.add_argument('--conditionA', action='store', dest='conditionA', required=False, default='',
            help='''Implies --test. Specify one condition corresponding to samples in the counts_matrix to be compared against
            condition2; by default, the first two unique conditions are used''')
    parser.add_argument('--conditionB', action='store', dest='conditionB', required=False, default='',
            help='''Specify another condition corresponding to samples in the counts_matrix to be compared against
            conditionA''')
    parser.add_argument('-of', '--out_dir_force', action='store_true',
            required=False, help='''Specify this argument to force overwriting of files in
            an existing output directory''')
    args = parser.parse_args()

    if isoforms:
        args.isoforms = isoforms
        args.counts_matrix = counts_matrix

    if not os.path.exists(args.counts_matrix):
        raise FlairInputDataError('Counts matrix file path does not exist')
    if not os.path.exists(args.isoforms):
        raise FlairInputDataError('Isoform bed file path does not exist')

    # Create output directory including a working directory for intermediate files.
    workdir = os.path.join(args.out_dir, 'workdir')
    if args.out_dir_force:
        if not os.path.exists(workdir):
            os.makedirs(workdir)
        pass
    elif not os.path.exists(args.out_dir):
        try:
            os.makedirs(workdir, 0o700)
        except OSError as ex:
            raise OSError("** ERROR cannot create directory %s" % (workdir)) from ex
    else:
        raise FlairInputDataError(f'** Error. Name {args.out_dir} already exists. Choose another name for out_dir')
    if args.isoforms.endswith('psl'):
        raise FlairInputDataError('** Error. Flair no longer accepts PSL input. Please use psl_to_bed first.')

    filebase = os.path.join(args.out_dir, 'diffsplice')
    pipettor.run(['call_diffsplice_events.py', args.isoforms, filebase, args.counts_matrix])
    pipettor.run(['es_as.py', args.isoforms], stdout=open(filebase+'.es.events.tsv', 'w'))
    pipettor.run(['es_as_inc_excl_to_counts.py', args.counts_matrix, filebase+'.es.events.tsv'],
                 stdout=open(filebase+'.es.events.quant.tsv', 'w'))
    os.unlink(filebase+'.es.events.tsv')

    if args.test or args.conditionA:
        logging.info('DRIMSeq testing for each AS event type')
        ds_command = ['Rscript', diffSplice_drimSeq, '--threads', args.threads, '--outDir', args.out_dir,
                      '--drim1', args.drim1, '--drim2', args.drim2, '--drim3', args.drim3, '--drim4', args.drim4]
        if args.batch:
            ds_command += ['--batch']
        if args.conditionA:
            if not args.conditionB:
                logging.info('Both conditionA and conditionB must be specified, or both left unspecified')
                return 1
            ds_command += ['--conditionA', args.conditionA, '--conditionB', args.conditionB]

        with open(workdir+'/ds.stderr.txt', 'w') as ds_stderr:
            for event in ['es', 'alt5', 'alt3', 'ir']:
                matrixfile = f'{filebase}.{event}.events.quant.tsv'
                if emptyMatrix(matrixfile):
                    logging.info(f'{event} event matrix file empty, not running DRIMSeq\n')
                    continue
                cur_command = ds_command + ['--matrix', matrixfile, '--prefix', event]
                try:
                    pipettor.run(cur_command, stderr=ds_stderr)
                except pipettor.ProcessException as exc:
                    raise FlairError(f"DRIMSeq failed on `{event}' event"
                                     f'Check {workdir}/ds.stderr.txt for details') from exc
    else:
        # workdir only necessary for drimseq output
        os.rmdir(workdir)
    return 0

def emptyMatrix(infile):
    '''Returns true if file has only a header line'''
    with open(infile, 'r') as inf:
        if len(inf.readlines()) <= 1:
            return True
    return False


if __name__ == '__main__':
    exit(diffSplice())
