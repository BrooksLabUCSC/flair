#! /usr/bin/env python3

import sys
import argparse
import subprocess
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'

def diffSplice(isoforms='', counts_matrix=''):
    parser = argparse.ArgumentParser(description='flair-diffSplice parse options',
            usage='flair diffSplice -i isoforms.bed -q counts_matrix.tsv [options]')
    required = parser.add_argument_group('required named arguments')
    if not isoforms:
        required.add_argument('-i', '--isoforms', action='store', dest='i', required=True,
                type=str, help='isoforms in bed format')
        required.add_argument('-q', '--counts_matrix', action='store', dest='q',
                type=str, required=True, help='tab-delimited isoform count matrix from flair quantify module')
    required.add_argument('-o', '--out_dir', action='store', dest='o',
            type=str, required=True, help='Output directory for tables and plots.')
    parser.add_argument('-t', '--threads', action='store', dest='t',
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
    parser.add_argument('-of', '--out_dir_force', action='store_true', dest='of',
            required=False, help='''Specify this argument to force overwriting of files in
            an existing output directory''')
    args = parser.parse_args()

    if isoforms:
        args.i = isoforms
        args.q = counts_matrix

    if not os.path.exists(args.q):
        sys.stderr.write('Counts matrix file path does not exist\n')
        return 1
    if not os.path.exists(args.i):
        sys.stderr.write('Isoform bed file path does not exist\n')
        return 1

    # Create output directory including a working directory for intermediate files.
    workdir = os.path.join(args.o, 'workdir')
    if args.of:
        if not os.path.exists(workdir):
            os.makedirs(workdir)
        pass
    elif not os.path.exists(args.o):
        try:
            os.makedirs(workdir, 0o700)
        except OSError as ex:
            raise OSError("** ERROR cannot create directory %s" % (workdir)) from ex
    else:
        sys.stderr.write(f'** Error. Name {args.o} already exists. Choose another name for out_dir\n')
        return 1
    if args.i.endswith('psl'):
        sys.stderr.write('** Error. Flair no longer accepts PSL input. Please use psl_to_bed first.\n')
        return 1

    filebase = os.path.join(args.o, 'diffsplice')
    subprocess.check_call(['call_diffsplice_events.py', args.i, filebase, args.q])
    subprocess.check_call(['es_as.py', args.i], stdout=open(filebase+'.es.events.tsv', 'w'))
    subprocess.check_call(['es_as_inc_excl_to_counts.py', args.q, filebase+'.es.events.tsv'],
            stdout=open(filebase+'.es.events.quant.tsv', 'w'))
    subprocess.check_call(['rm', filebase+'.es.events.tsv'])

    if args.test or args.conditionA:
        sys.stderr.write('DRIMSeq testing for each AS event type\n')
        drim1, drim2, drim3, drim4 = [str(x) for x in [args.drim1, args.drim2, args.drim3, args.drim4]]
        ds_command = ['runDS.py', '--threads', str(args.t), '--outDir', args.o,
                '--drim1', drim1, '--drim2', drim2, '--drim3', drim3, '--drim4', drim4]
        if args.batch:
            ds_command += ['--batch']
        if args.conditionA:
            if not args.conditionB:
                sys.stderr.write('Both conditionA and conditionB must be specified, or both left unspecified\n')
                return 1
            ds_command += ['--conditionA', args.conditionA, '--conditionB', args.conditionB]

        with open(workdir+'/ds.stderr.txt', 'w') as ds_stderr:
            for event in ['es', 'alt5', 'alt3', 'ir']:
                matrixfile = f'{filebase}.{event}.events.quant.tsv'
                if emptyMatrix(matrixfile):
                    sys.stderr.write(f'{event} event matrix file empty, not running DRIMSeq\n')
                    continue
                cur_command = ds_command + ['--matrix', matrixfile, '--prefix', event]
                if subprocess.call(cur_command, stderr=ds_stderr):
                    print(f'\nDRIMSeq failed on {event} event')
                    print(f'Check {workdir}/ds.stderr.txt for details.\nCommand was:')
                    print(' '.join([x for x in cur_command]))
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
