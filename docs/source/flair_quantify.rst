flair quantify
==============

.. code:: sh

    usage: python flair.py quantify -r reads_manifest.tsv -i isoforms.fa [options]

Output isoform-by-sample counts file that can be used in the diffExp and diffSplice modules.

Options
-------

Required arguments
~~~~~~~~~~~~~~~~~~
``--isoforms`` Fasta of Flair collapsed isoforms

``--reads_manifest`` Tab delimited file containing sample id, condition, batch, reads.fq, where ``reads.fq`` is the path to the sample fastq file. Example:

.. code:: sh

   sample1      condition1      batch1  mydata/sample1.fq
   sample2      condition1      batch1  mydata/sample2.fq
   sample3      condition1      batch1  mydata/sample3.fq
   sample4      condition2      batch1  mydata/sample4.fq
   sample5      condition2      batch1  mydata/sample5.fq
   sample6      condition2      batch1  mydata/sample6.fq

Note: Do **not** use underscores in the first three fields, see below for details.


Optional arguments
~~~~~~~~~~~~~~~~~~
``--help`` Show all options

``--threads`` Number of processors to use (default 4).

``--output`` Name base for output files (default: ``flair.quantify``). You can supply an output directory (e.g. ``output/flair_quantify``)

``--temp_dir`` directory to put temporary files. use "./" to indicate current directory (default: python tempfile directory)

``--salmon`` Path to salmon executable, specify if salmon quantification is desired. Please note that salmon is not installed with Flair's conda or docker installations.

``--quality`` Minimum MAPQ of read assignment to an isoform (default 1). If using salmon, all alignments are used.

``--trust_ends`` Specify if reads are generated from a long read method with minimal fragmentation

``--sample_id_only`` Only use sample id in output header instead of a concatenation of id, condition, and batch

``--generate_map`` Create read-to-isoform assignment files for each sample

``--tpm`` Convert counts matrix to transcripts per million and output as a separate file named <output>.tpm.tsv

``--isoform_bed`` isoform .bed or .psl file, must be specified if --stringent or --check-splice is specified

``--stringent`` Supporting reads must cover 80% of their isoform and extend at least 25 nt into the first and last exons. If those exons are themselves shorter than 25 nt, the requirement becomes 'must start within 4 nt from the start' or 'end within 4 nt from the end'

``--check_splice`` Enforces coverage of 4 out of 6 bp around each splice site and no insertions greater than 3 bp at the splice site


Obsolete arguments (will be removed soon)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you ``conda install`` Flair or use its docker container, minimap2 and samtools are in its $PATH

``--minimap2`` Path to minimap2 if not in $PATH.

``--samtools`` Path to samtools if not in $PATH (only needed when --quality is specified)


Unless ``--sample_id_only`` is specified, the output counts file concatenates id, condition and batch info for each sample. Flair diffExp and diffSplice expect this information.

.. code:: sh

id   sample1_condition1_batch1  sample2_condition1_batch1  sample3_condition1_batch1  sample4_condition2_batch1  sample5_condition2_batch1  sample6_condition2_batch1
ENST00000225792.10_ENSG00000108654.15   21.0    12.0    10.0    10.0    14.0    13.0
ENST00000256078.9_ENSG00000133703.12    7.0     6.0     7.0     15.0    12.0    7.0


