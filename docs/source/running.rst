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

The additional programs, such as ``diff_iso_usage`` and ``bed_to_gtf``, name their
output file rather than a prefix, so they take no option: end the name in ``.gz``
and the output is gzipped.
