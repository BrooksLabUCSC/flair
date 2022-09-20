flair diffSplice
=======

usage: python flair.py diffSplice -i isoforms.bed|isoforms.psl -q counts_matrix.tsv [options]

flair-diffSplice parse options

positional arguments:
  diffSplice

optional arguments:
  -h, --help            show this help message and exit
  -o O, --output O      output file name base for FLAIR isoforms (default:
                        flair.diffsplice)
  --test                Run DRIMSeq statistical testing
  -t T, --threads T     Number of threads for parallel DRIMSeq (1)
  --drim1 DRIM1         The minimum number of samples that have coverage over
                        an AS event inclusion/exclusion for DRIMSeq testing;
                        events with too few samples are filtered out and not
                        tested (6)
  --drim2 DRIM2         The minimum number of samples expressing the inclusion
                        of an AS event; events with too few samples are
                        filtered out and not tested (3)
  --drim3 DRIM3         The minimum number of reads covering an AS event
                        inclusion/exclusion for DRIMSeq testing, events with
                        too few samples are filtered out and not tested (15)
  --drim4 DRIM4         The minimum number of reads covering an AS event
                        inclusion for DRIMSeq testing, events with too few
                        samples are filtered out and not tested (5)
  --batch               If specified with --test, DRIMSeq will perform batch
                        correction
  --conditionA CONDITIONA
                        Specify one condition corresponding to samples in the
                        counts_matrix to be compared against condition2; by
                        default, the first two unique conditions are used
  --conditionB CONDITIONB
                        Specify another condition corresponding to samples in
                        the counts_matrix to be compared against conditionA

required named arguments:
  -i I, --isoforms I    isoforms in bed or psl format
  -q Q, --counts_matrix Q
                        tab-delimited isoform count matrix from flair quantify
                        module
