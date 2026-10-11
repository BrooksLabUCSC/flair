.. _running-label:

Running FLAIR
^^^^^^^^^^^^^

``flair`` is one program and each analysis is a sub-command, as in
``flair transcriptome -o sample1 ...``.  The options below are global and come
before the sub-command.  Giving one after it is an error that says so.

.. include:: cli/flair.rst

.. _compression-label:

Compression
===========

Any text file FLAIR reads may be gzipped, in any sub-command.  FLAIR decides by
the file name, so ``reads.fa.gz`` and ``annotation.gtf.gz`` are read as they are.

The genome is the exception.  FLAIR needs random access to it, so a compressed
genome must be bgzip compressed and indexed::

   bgzip genome.fa
   samtools faidx genome.fa.gz

A genome that is plain gzip, or bgzip without its index, is reported before the
run starts.

With ``--compress``, or ``-z``, the tables FLAIR writes gain a ``.gz`` suffix.
BAM, bigBed and the PDF plots stay plain, as do the isoform BED, GTF and FASTA,
which genome browsers and ``bedToBigBed`` read and which do not take gzip.

Reading and writing are independent: a run may read gzipped input and write
plain output, or the reverse.

The additional programs, such as ``diff_iso_usage`` and ``bed_to_gtf``, name
their output file rather than a prefix, so they take no option: end the name in
``.gz`` and the output is gzipped.
