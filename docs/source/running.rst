.. _running-label:

Running FLAIR
^^^^^^^^^^^^^

``flair`` is one program and each analysis is a sub-command, as in
``flair transcriptome -o sample1 ...``.  The options below are global and come
before the sub-command.  Giving one after it is an error that reports it.

.. include:: cli/flair.rst

.. _logging-label:

What a run prints
=================

Every message goes to standard error, so the output files of a run and its log are
separate.  A record is stamped with the time it was logged, and a warning with its
level.

At the default level a run reports its version and command line, each stage as it
starts, each partition or batch as it finishes, and what it produced::

   2026-11-19 15:28:57 FLAIR 3.0.0: flair transcriptome -b reads.bam -g genome.fa ...
   2026-11-19 15:28:57 loading annotation GTF
   2026-11-19 15:28:57 running partitions
   2026-11-19 15:29:41 partition 1 of 3: chr20:35542020-35557634, 34 items, 43.9 s
   2026-11-19 15:30:02 7 isoforms
   2026-11-19 15:30:02 flair transcriptome took 0:01:05

The partition lines are the progress of a long run: a partition is a cluster of
reads and annotated transcripts, and no read crosses one, so the count is of the
work and not of the genome.

``--log-level`` takes the level names of the Python logging module.  ``WARNING``
prints only the warnings and the error that ends a run, ``DEBUG`` prints every
program FLAIR runs with its arguments, every read dropped and the reason for it,
and the diagnostics of the FLAIR programs a run starts.  ``--log-debug`` is a
short form of ``--log-level DEBUG``, and its output is what to send with a bug
report::

   flair --log-debug transcriptome ... 2> run.log

An error is printed whatever level was given.

.. _compression-label:

Compression
===========

Any text file FLAIR reads may be gzipped, in any sub-command, decided by the file
name extension.  The genome needs random access, so a compressed genome must be
bgzip compressed and indexed::

   bgzip genome.fa
   samtools faidx genome.fa.gz

With ``--compress``, or ``-z``, the files FLAIR writes gain a ``.gz`` suffix: the
TSVs, the isoform BED, the GTF, the FASTA and the VCF.  The FASTA is bgzip, so
``samtools faidx`` can index it, and the VCF is bgzip with a ``.tbi`` beside it.

Reading and writing are independent: a run may read gzipped input and write plain
output, or the reverse.

The additional programs name their output file rather than a prefix, so they take no
option: end the name in ``.gz`` and the output is gzipped, as with
``gtf_to_bed ANNOTATION.gtf ANNOTATION.bed.gz``.  ``bed_to_gtf`` is the exception.  It
writes to standard output, so redirect it through ``gzip``::

   bed_to_gtf ISOFORMS.bed | gzip > ISOFORMS.gtf.gz
