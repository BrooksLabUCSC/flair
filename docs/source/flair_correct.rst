flair correct
=============

.. code:: sh

usage: python flair.py correct -q query.bed12 [-f annotation.gtf]|[-j introns.tab] -g genome.fa [options]


This module corrects misaligned splice sites using genome annotations and/or short-read splice junctions. It creates a ``<prefix>.corrected.bed`` for use in subsequent steps and a ``<prefix>.inconsistent.bed`` with rejected alignments.

Make sure that the genome annotation and genome sequences are compatible (if the genome sequence contains the 'chr' prefix, the annotations must too).

Please do use GTF instead of GFF; annotations should not split single exons into multiple entries. 


Options
-------

Required arguments
~~~~~~~~~~~~~~~~~~
``--query`` Uncorrected bed12 file, e.g. output of ``flair align``.

``--genome`` Reference genome in fasta format.

At least one of the following arguments is required:

``--shortread`` Bed format splice junctions from short-read sequencing. You can generate these from ``sam`` format files using the ``junctions_from_sam`` program that comes with Flair.

``--gtf`` GTF annotation file.

Optional arguments
~~~~~~~~~~~~~~~~~~
``--help`` Show all options 

``--chromsizes`` Tab separated file of chromosome sizes. `Here <https://raw.githubusercontent.com/igvteam/igv/master/genomes/sizes/hg38.chrom.sizes>`__ is one for GRCh38/hg38. If you use this parameter, a ``<prefix>_corrected.psl`` file will be created in addition to `<prefix)_corrected.bed``

``--nvrna`` Specify this flag to make the strand of a read consistent with the input annotation during correction.

``--threads`` Number of processors to use (default 4).

``--ss_window`` Window size for correcting splice sites (default 15).

``-output`` Name base for output files (default: ``flair``). You can supply an output directory (e.g. ``output/flair``)
but it has to exist; Flair will not create it. If you run the same command twice, Flair will overwrite the files without warning.

``--print_check``         Print err.txt with step checking.


