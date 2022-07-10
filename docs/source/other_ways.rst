Other ways to run FLAIR modules
===============================

For convenience, multiple FLAIR modules can be run in the same command.
In place of a single module name, multiple module numbers can be
specified (module numbers: align=1, correct=2, collapse=3,
collapse-range=3.5, quantify=4, diffExp=5, diffSplice=6). All arguments
for the modules that will be run must be provided. For example, to run
the align, correct, and collapse modules, the command might look like:

.. code:: sh

   python flair.py 123 -r reads.fa -g genome.fa -f annotation.gtf -o flair.output --temp_dir temp_flair [optional arguments]

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

   python flair.py collapse-range -r reads.bam -q query.bed -g genome.fa -f annotation.gtf -o flair.output --temp_dir temp_flair [optional arguments]

If the user would prefer not to use python’s multiprocessing module, a
bash script has also been provided
(``bin/run_flair_collapse_ranges.sh``) that runs collapse-range for the
user that parallelizes using GNU parallel, which the user can alter as
they see fit for their system.

