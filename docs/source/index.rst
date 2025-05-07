.. Flair documentation master file, created by
   sphinx-quickstart on Fri Jul  8 08:42:03 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to FLAIR's Documentation!
=================================

Staying in touch and getting help
=================================

**Please subscribe to the FLAIR announce mailing list at:**

FLAIR Announce Group

   https://groups.google.com/a/ucsc.edu/g/flair-announce-group

This is a read-only, low volume list that will only have announcements of new
FLAIR releases, publications, and other FLAIR-related user information.

If you have questions, please ask in FLAIR GitHub Discussions:

   https://github.com/BrooksLabUCSC/flair/discussions

Report bugs to FLAIR GitHub Issues:

   https://github.com/BrooksLabUCSC/flair/issues
   

Installing FLAIR
================

See  :ref:`installing-label` for instructions on the various approache to install FLAIR.


Using FLAIR
===========

FLAIR can be run optionally with short-read data to help increase splice
site accuracy of the long read splice junctions. FLAIR uses multiple
alignment steps and splice site filters to increase confidence in the
set of isoforms defined from noisy data. FLAIR was designed to be able
to sense subtle splicing changes in nanopore data from `Tang et
al. (2020) <https://www.nature.com/articles/s41467-020-15171-6>`__.
Please read for more description of the methods.

.. figure:: img/flair_workflow_compartmentalized.png


If you have multiple samples and want to compare them on a single
transcriptome, you have two options:

Run flair correct and collapse individually on each sample, then combine your
transcriptomes using ``flair combine``. This method
will be faster and easier, but you may miss some low-expression transcripts.

Your other option is to merge your samples before running FLAIR. If using
PacBio reads, be careful doing this, as the reads may not have unique IDs. You
may need to label each read with its sample ID to keep the read IDs
unique. You can either merge the FASTA/FASTQ files before running FLAIR
(simplest, recommended), or merge the bed files after running FLAIR correct,
making sure to run FLAIR collapse with a list of all of your read files.
            
It is recommended to combine all samples together prior to running
``flair collapse`` for isoform assembly by concatenating corrected read ``BED``
files together. Following the creation of an isoform reference from
``flair collapse``, consequent steps will assign reads from each sample
individually to isoforms of the combined assembly for downstream analyses.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installing.rst
   modules.rst
   scripts.rst
   flair2_functions.rst
   other_ways.rst
   testrun.rst
   example_files.rst
   faqs.rst
   cite.rst



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
