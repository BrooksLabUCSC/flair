flair align
===========

.. code:: sh

    usage: python flair.py align -g genome.fa -r <reads.fq>|<reads.fa> [options]


This module aligns reads to the genome using `minimap2 <https://github.com/lh3/minimap2>`__, 
and converts the ``sam`` output to `BED12 <https://genome.ucsc.edu/FAQ/FAQformat.html#format14>`__ and
optionally `PSL <https://genome.ucsc.edu/FAQ/FAQformat.html#format2>`__ formats.  
Aligned reads in ``psl`` format can be visualized in `IGV <https://igv.org/>`__ or the 
`UCSC Genome browser <https://genome.ucsc.edu/cgi-bin/hgGateway>`__. 

If you're using human sequences, the best reference genome is 
`GCA_000001405.15_GRCh38_no_alt_analysis_set <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz>`__ as described in this 
`helpful blog post by Heng Li <https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use>`__

If your input sequences are Oxford nanopore reads, please use `Pychopper <https://github.com/epi2me-labs/pychopper>`__ before running Flair.

If your reads are already aligned, you can convert the sorted ``bam`` output to ``bed12`` using
``bam2Bed12`` to supply for flair-correct. This step smoothes gaps in the alignment.


Options
-------

Required arguments
~~~~~~~~~~~~~~~~~~
``-reads`` Raw reads in fasta or fastq format. This argument accepts multiple files if they are separated by spaces.

At least one of the following arguments is required:

``-genome`` Reference genome in fasta format. Flair will minimap index this file unless there already is a ``.mmi`` file in the same location.

``-mm_index`` If there already is a ``.mmi`` index for the genome it can be supplied directly using this option. 


Optional arguments
~~~~~~~~~~~~~~~~~~

``-help`` Show all options.

``-output`` Name base for output files (default: ``flair.aligned``). You can supply an output directory (e.g. ``output/flair_aligned``) 
but it has to exist; Flair will not create it. If you run the same command twice, Flair will overwrite the files without warning.

``--threads`` Number of processors to use (default 4).

``--junction_bed`` Annotated isoforms/junctions bed file for splice site-guided minimap2 genomic alignment.

``--nvrna`` Use native-RNA specific alignment parameters for minimap2 (``-u f -k 14``, see `minimap2's manual <https://lh3.github.io/minimap2/minimap2.html>`__ for details).

``--psl`` Output PSL format file in addition to bed. PSL can be used as input to other Flair modules as well.

``--chromsizes`` Tab separated file of chromosome sizes, needed to make the ``psl`` file genome browser compatible. `Here <https://raw.githubusercontent.com/igvteam/igv/master/genomes/sizes/hg38.chrom.sizes>`__ is one for GRCh38/hg38.

``--quality`` Minimum `MAPQ score <http://www.acgt.me/blog/2014/12/16/understanding-mapq-scores-in-sam-files-does-37-42>`__ of read alignment to the genome. The default is 1, which is the lowest possible score.

``-N`` Retain at most INT secondary alignments from minimap2 (default 0). Please proceed with caution, changing this setting is only useful if you know there are closely related homologs elsewhere in the genome. It will likely decrease the quality of Flair's final results.

``--quiet`` Dont print progress statements.

Obsolete arguments (will be removed soon)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you ``conda install`` Flair or use its docker container, minimap2 and samtools are in $PATH.

``--minimap2`` Path to minimap2 if not in $PATH.

``--samtools`` Path to samtools if not in $PATH

``--pychopper`` `Pychopper <https://github.com/epi2me-labs/pychopper>`__ is a preprocessing tool for Oxford Nanopore reads. It's easy to conda install and run locally; please do so before running Flair.


