.. Flair documentation master file, created by
   sphinx-quickstart on Fri Jul  8 08:42:03 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to FLAIR's documentation!
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

FLAIR can be conda installed using 

.. code:: sh

   conda create -n flair -c conda-forge -c bioconda flair
   conda activate flair

On Apple Silicon Mac systems (ARM64 CPUs: M1, M2, ...) you must use

.. code:: sh

   CONDA_SUBDIR=osx-64 conda create -n flair
   conda activate flair
   conda config --env --set subdir osx-64
   conda config --add channels bioconda
   conda config --add channels conda-forge
   conda install flair

Note that mamba currently fails to install FLAIR on Mac ARM64.

FLAIR can be run optionally with short-read data to help increase splice
site accuracy of the long read splice junctions. FLAIR uses multiple
alignment steps and splice site filters to increase confidence in the
set of isoforms defined from noisy data. FLAIR was designed to be able
to sense subtle splicing changes in nanopore data from `Tang et
al. (2020) <https://www.nature.com/articles/s41467-020-15171-6>`__.
Please read for more description of the methods.

.. figure:: img/flair_workflow_compartmentalized.png


It is recommended to combine all samples together prior to running
flair-collapse for isoform assembly by concatenating corrected read
``bed`` files together. Following the creation of an isoform
reference from flair-collapse, consequent steps will assign reads from
each sample individually to isoforms of the combined assembly for
downstream analyses.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   requirements.rst
   modules.rst
   scripts.rst
   flair2_functions.rst
   other_ways.rst
   other_environments.rst
   testrun.rst
   example_files.rst
   faqs.rst
   cite.rst



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
