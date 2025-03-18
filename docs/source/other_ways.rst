Other ways to run FLAIR modules
===============================

A beta version of the collapse module, called collapse-range, has been
developed. The corrected reads are divided into many independent
regions, which are then subject to isoform calling separately and
parallelized over the number of threads specified. This dramatically
decreases the memory footprint of intermediate files and increases the
speed in which the module runs without altering the final isoforms. This
version can be invoked by specifying collapse-range as the module (or
3.5 if using numbers). An additional program,
`bedPartition <http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/>`__,
needs to be in your $PATH.

.. code:: sh

   flair collapse-range -r reads.bam -q query.bed -g genome.fa -f annotation.gtf -o flair.output --temp_dir temp_flair [optional arguments]

If you would prefer not to use python’s multiprocessing module, a
bash script has also been provided
(``run_flair_collapse_ranges.sh``) that runs collapse-range for the
user that parallelizes using GNU parallel, which you can alter as
they see fit for their system.

