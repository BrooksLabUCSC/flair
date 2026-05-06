Modules
^^^^^^^

``flair`` is a program which implement the follow module sub-commands.

If you want to compare multiple samples, there are two primary ways of doing this:
 - Combine the fastq or fasta reads of all samples and run FLAIR transcriptome on all samples together (will generate the most comprehensive transcriptome)
 - Run FLAIR transcriptome on each sample separately (better for large sets of samples) and
   use FLAIR combine to merge results

.. figure:: img/flair_workflow_compartmentalized.svg


.. figure:: img/flair_workflow_compartmentalized.svg


.. _transcriptome-label:

flair transcriptome
===================

.. code:: text

    usage: usage: flair transcriptome -b reads.genomealigned.bam [options]


This module generates a transcriptome of high confidence isoforms (bed, gtf, and fasta files) directly from a bam file of aligned reads.
To get aligned reads, you can use FLAIR align or just run the following command to generate the bam file to use as input.
minimap2 -ax splice -s 80 -G 200k -t 20 --secondary=no genome.fa sample.fastq | samtools view -hb - | samtools sort - > sample.genomealigned.bam; samtools index sample.genomealigned.bam
If you want to run downstream fusion detection with FLAIR fusion, run flair align with --filtertype separate to generate a separate file of chimeric alignments.


**Outputs**

 - ``flair.isoforms.bed``
 - ``flair.isoforms.gtf``
 - ``flair.isoforms.fa``
 - ``flair.read.map.txt``

Options
-------

Required arguments
~~~~~~~~~~~~~~~~~~

.. code:: text

    -b --genomealignedbam    Sorted and indexed bam file aligned to the genome
    -g --genome    Reference genome in fasta format


Optional arguments
~~~~~~~~~~~~~~~~~~

.. code:: text

 -o --output
                        output file name base for FLAIR isoforms (default: flair.collapse)
  -t --threads
                        minimap2 number of threads (4)
  -f --gtf              [HIGHLY RECOMMENDED] GTF annotation file, used for renaming FLAIR isoforms 
                        to annotated isoforms and adjusting TSS/TESs
  -j --shortread
                        [HIGHLY RECOMMENDED] bed format splice junctions from short-read sequencing. 
                        NO NOVEL SPLICE SITES WILL BE DETECTED WITHOUT ORTHOGONAL SHORT READS
  --ss_window
                        window size for correcting splice sites (15)
  -s  --support
                        minimum number of supporting reads for an isoform (3)
  --stringent           [HIGHLY RECOMMENDED] specify if all supporting reads need to be full-length 
                        (spanning 25 bp of the first and last exons)
  --check_splice        [HIGHLY RECOMMENDED] enforce coverage of 4 out of 6 bp around each splice site 
                        and no insertions greater than 3 bp at the splice site. DON'T USE WITH DATA WITH ERROR RATES of 5% or more, such
                        as older direct-RNA (anything before using the RNA-specific flow-cell with the Dorado basecaller).
  -w --end_window
                        window size for comparing TSS/TES (100)
  --noaligntoannot      related to old annotation_reliant, now specify if you don't want an initial alignment 
                        to the annotated sequences and only want transcript detection from the
                        genomic alignment. Will be slightly faster but less accurate if the annotation is good
  -n --no_redundant 
                        For each unique splice junction chain, report options include: none--best TSSs/TESs chosen for each unique set of splice junctions; longest--single TSS/TES
                        chosen to maximize length; best_only--single most supported TSS/TES used in conjunction chosen (none)
  --max_ends            maximum number of TSS/TES picked per isoform (2)
  --filter              Report options include: 
                            default--subset isoforms are removed based on support;
                            nosubset--any isoforms that are a proper set of another isoform are removed;
                            comprehensive--default set + all subset isoforms; 
                            ginormous--comprehensive set + single exon subset isoforms
  --splittoregion       force running on each region of non-overlapping reads, no matter the file size 
                        default: parallelize by chromosome if file is <1G, otherwise parallelize on all regions of non-overlapping reads
  --predictCDS          specify if you want to predict the CDS of the final isoforms. 
                        Will be output in the final bed file but not the gtf file. 
                        Productivity annotation is also added in the name field, 
                        which is detailed further in the predictProductivity documentation



.. _align-label:

flair align
===========

.. code:: text

    usage: flair align -g genome.fa -r <reads.fq>|<reads.fa> [options]

Use of this modules is deprecated, as ``FLAIR transcriptome`` operates on a BAM, and other output of ``flair align`` is no longer used.

This module aligns reads to the genome using `minimap2 <https://github.com/lh3/minimap2>`__, 
and converts the `SAM <https://en.wikipedia.org/wiki/SAM_(file_format)>`__ output to `BED12 <https://genome.ucsc.edu/FAQ/FAQformat.html#format14>`__.
Aligned reads in BED12 format can be visualized in `IGV <https://igv.org/>`__ or the 
`UCSC Genome browser <https://genome.ucsc.edu/cgi-bin/hgGateway>`__. 


**Outputs**

 - ``flair.aligned.bam``
 - ``flair.aligned.bam.bai``
 - ``flair.aligned.bed``

Options
-------

Required arguments
~~~~~~~~~~~~~~~~~~

.. code:: text

    --reads     Raw reads in fasta or fastq format. This argument accepts multiple 
                (comma/space separated) files.

    At least one of the following arguments is required:
    --genome    Reference genome in fasta format. Flair will minimap index this file 
                unless there already is a .mmi file in the same location.
    --mm_index  If there already is a .mmi index for the genome it can be supplied 
                directly using this option. 


Optional arguments
~~~~~~~~~~~~~~~~~~

.. code:: text

  -o OUTPUT, --output OUTPUT
                        output file name base (default: flair.aligned)
  -t THREADS, --threads THREADS
                        minimap2 number of threads (4)
  --junction_tab JUNCTION_TAB
                        short-read junctions in SJ.out.tab format. Use this option
                        if you aligned your short-reads with STAR, STAR will
                        automatically output this file
  --junction_bed JUNCTION_BED
                        annotated isoforms/junctions bed file for splice site-guided minimap2 genomic alignment
  --nvrna               specify this flag to use native-RNA specific alignment parameters for minimap2
  --quality QUALITY     minimum MAPQ of read alignment to the genome (0)
  --minfragmentsize MINFRAGMENTSIZE
                        minimum size of alignment kept, used in minimap -s. More important when doing downstream fusion detection
  --maxintronlen MAXINTRONLEN
                        maximum intron length in genomic alignment. Longer can help recover more novel isoforms with long introns
  --filtertype FILTERTYPE
                        method of filtering chimeric alignments (potential fusion reads). Options: removesup (default), separate (required for downstream work with fusions), keepsup
                        (keeps supplementary alignments for isoform detection, does not allow gene fusion detection)
  --quiet               Suppress minimap progress statements from being printed
  --remove_internal_priming
                        specify if want to remove reads with internal priming
  -f GTF, --gtf GTF     reference annotation, only used if --remove_internal_priming is specified, recommended if so
  --intprimingthreshold INTPRIMINGTHRESHOLD
                        number of bases that are at leas 75% As required to call read as internal priming
  --intprimingfracAs INTPRIMINGFRACAS
                        number of bases that are at least 75% As required to call read as internal priming
  --remove_singleexon   specify if want to remove unspliced reads
    

Notes
-----
If you're using human sequences, the best reference genome is 
`GCA_000001405.15_GRCh38_no_alt_analysis_set <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz>`__ as described in this 
`helpful blog post by Heng Li <https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use>`__

If your input sequences are Oxford nanopore reads, please use `Pychopper <https://github.com/epi2me-labs/pychopper>`__ before running Flair.

nvrna settings: See `minimap2's manual <https://lh3.github.io/minimap2/minimap2.html>`__ for details.

quality: `More info on MAPQ scores <http://www.acgt.me/blog/2014/12/16/understanding-mapq-scores-in-sam-files-does-37-42>`__ 

.. _fusions-label:

flair fusion
============

.. code:: text

    usage: flair fusion -g genome.fa -r sample.fastq -b sample.genomealigned_chimeric.bam -f annot.gtf [-o OUTPUT_PREFIX]

This identifies gene fusions and generates a fusion transcriptome. 
To incorporate this fusion transcriptome in downstream analysis, 
use flair combine to merge it with normal isoforms. 

**Output**

sample.fusions.isoforms.bed
    Bed file of fusion transcriptome (each fusion has a line for each locus in the fusion, 
    and position in the fusion is specified by the fusiongeneX prefix in the name field
sample.fusions.isoforms.fa
    Fasta file of fusion transcriptome
sample.syntheticAligned.isoform.read.map
    read map of reads to fusion isoforms

Required Options
----------------

.. code:: text

  -g --genome
                        FastA of reference genome
  -r READS [READS ...], --reads READS [READS ...]
                        FastA/FastQ files of raw reads, can specify multiple files
  -b --genomechimbam
                        bam file of chimeric reads from genomic alignment from flair align run with --filtertype separate
  -f --gtf              GTF annotation file

Other Options
-------------

.. code:: text

  --transcriptchimbam TRANSCRIPTCHIMBAM
                        Optional: bam file of chimeric reads from transcriptomic alignment. 
                        If not provided, this will be made for you
  -o OUTPUT, --output OUTPUT
                        output file name base for FLAIR isoforms
  -t --threads
                        minimap2 number of threads (4)
  --minfragmentsize 
                        minimum size of alignment kept, used in minimap -s (40)
  -s --support
                        minimum number of supporting reads for a fusion (3)
  --maxloci             max loci detected in fusion. Set higher for detection of 3-gene+ fusions

  --keep_intermediate   keep intermediate and temporary files debugging purposes

.. _combine-label:

flair combine
=============
.. code:: sh

    usage: flair_combine [-h] -m MANIFEST [-o OUTPUT_PREFIX] [-w ENDWINDOW]
                         [-p MINPERCENTUSAGE] [-c] [-s] [-f FILTER]

    options:
      -h, --help            show this help message and exit
      -m MANIFEST, --manifest MANIFEST
                            path to manifest files that points to transcriptomes to combine.
                            Each line of file should be tab separated with sample name, sample
                            type (isoform or fusionisoform), path/to/isoforms.bed,
                            path/to/isoforms.fa, path/to/isoform.read.map.txt. fa and
                            read.map.txt files are not required, although if .fa files are not
                            provided for each sample a .fa output will not be generated
      -o OUTPUT_PREFIX, --output_prefix OUTPUT_PREFIX
                            path to collapsed_output.bed file. default: 'collapsed_flairomes'
      -w ENDWINDOW, --endwindow ENDWINDOW
                            window for comparing ends of isoforms with the same intron chain.
                            Default:200bp
      -p MINPERCENTUSAGE, --minpercentusage MINPERCENTUSAGE
                            minimum percent usage required in one sample to keep isoform in
                            combined transcriptome. Default:10
      -c, --convert_gtf     [optional] whether to convert the combined transcriptome bed file
                            to gtf
      -s, --include_se      whether to include single exon isoforms. Default: dont include
      -f FILTER, --filter FILTER
                            type of filtering. Options: usageandlongest(default), usageonly,
                            none, or a number for the total count of reads required to call an
                            isoform

    Combines FLAIR transcriptomes with other FLAIR transcriptomes or annotation transcriptomes to generate accurate combined transcriptome. Only the manifest file is required. Manifest file is in the following format. If the transcriptome is from FLAIR collapse or transcriptome, but isoform in the second column, if it is from FLAIR fusion, put fusionisoform in the second column:

Manifest example (we suggest using absolute file paths to point to your files though):

.. code:: text

    sample1	isoform	sample1.FLAIR.isoforms.bed	sample1.FLAIR.isoforms.fa	sample1.read.map.txt
    sample2	isoform	sample2.FLAIR.isoforms.bed	sample2.FLAIR.isoforms.fa	sample2.read.map.txt
    sample1	fusionisoform	sample1.fusion.isoforms.bed	sample1.fusion.isoforms.fa	sample1.fusion.isoform.read.map.txt
    sample2	fusionisoform	sample2.fusion.isoforms.bed	sample2.fusion.isoforms.fa	sample2.fusion.isoform.read.map.txt

For each line, the sample name and bed path is required. The fasta and
read.map.txt file is optional. Without these files there is less ability to
filter and more isoforms will be included. If a sample is a FLAIR run, we
highly recommend including the read.map.txt file. If you want to combine FLAIR
transcriptomes with annotated transcripts, you can convert an annotation gtf
file to a bed file using gtf_to_bed (see Additional Programs)

Flair combine will generate a counts file, but for the most accurate quantification, we recommend 
running FLAIR quantify using all samples against the combined transcriptome



.. _quantify-label:


flair quantify
==============

.. code:: text

    usage: flair quantify -r reads_manifest.tsv -i isoforms.fa [options]

**Output**

Default: identifes the best isoform assignment based on alignment quality, fraction of read aligned, and fraction of transcript aligned

check_splice: adds check for read matching reference transcript at all splice sites

stringent: adds requirement for read to cover at least 25bp of the first and last exons

If you need your reads to match your isoforms well, use --check_splice and --stringent, while if you need more reads assigned to isoforms for better statistical comparison, use the default.

--quality 0 is also reccommended, as this allows slightly better recall as FLAIR can disambiguate some similar isoform alignments.

Options
-------

Required arguments
~~~~~~~~~~~~~~~~~~

.. code:: text

    --isoforms          Fasta of Flair collapsed or combined isoforms
    --reads_manifest    Tab delimited file containing sample id, condition, batch, 
                        reads.fq, where reads.fq is the path to the sample fastq file. 

Manifest example (we suggest using absolute file paths to point to your files though):

.. code:: text

   sample1      condition1      batch1  mydata/sample1.fq
   sample2      condition1      batch1  mydata/sample2.fq
   sample3      condition1      batch1  mydata/sample3.fq
   sample4      condition2      batch1  mydata/sample4.fq
   sample5      condition2      batch1  mydata/sample5.fq
   sample6      condition2      batch1  mydata/sample6.fq

Note: Do **not** use underscores in the first three fields, see below for details.


Optional arguments
~~~~~~~~~~~~~~~~~~

.. code:: text

    --help	        Show all options
    --output	        Name base for output files (default: flair.quantify). You 
                        can supply an output directory (e.g. output/flair_quantify).
    --threads	        Number of processors to use (default 4).
    --temp_dir	        Directory to put temporary files. use ./ to indicate current 
                        directory (default: python tempfile directory).
    --sample_id_only	Only use sample id in output header instead of a concatenation 
                        of id, condition, and batch.
    --quality	        Minimum MAPQ of read assignment to an isoform (default 0). 
    --trust_ends	Specify if reads are generated from a long read method with 
                        minimal fragmentation.
    --generate_map	Create read-to-isoform assignment files for each sample.
    --isoform_bed	isoform .bed file, must be specified if --stringent or 
                        --check-splice is specified.
    --stringent	        Supporting reads must cover 80% of their isoform and extend 
                        at least 25 nt into the first and last exons. If those exons 
                        are themselves shorter than 25 nt, the requirement becomes 
                        'must start within 4 nt from the start' or 'end within 4 nt 
                        from the end'.
    --check_splice	Enforces coverage of 4 out of 6 bp around each splice site 
                        and no insertions greater than 3 bp at the splice site.
    --output_bam	If selected, forces output of each reads file aligned to the 
                        FLAIR transcriptome. This will be a bam with no secondary alignments

Other info
----------
Unless ``--sample_id_only`` is specified, the output counts file concatenates id, condition and batch info for each sample. The `flair diffexp` and `flair diffsplice` modules expect this information.

.. code:: text

   id   sample1_condition1_batch1  sample2_condition1_batch1  sample3_condition1_batch1  sample4_condition2_batch1  sample5_condition2_batch1  sample6_condition2_batch1
   ENST00000225792.10_ENSG00000108654.15   21.0    12.0    10.0    10.0    14.0    13.0
   ENST00000256078.9_ENSG00000133703.12    7.0     6.0     7.0     15.0    12.0    7.0



.. _variants-label:


flair variants
==============

.. code:: text

    usage: flair variants -m manifest.tsv -i isoforms.fa -b isoforms.bed -g genome.fa -f annot.gtf [-o OUTPUT_PREFIX]

This does not call variants, it integrates already called variants with 
isoforms to understand allele-specific isoform expression and allele bias.
Before running this module, you need to run a variant caller on each of your
samples individually. We recommend longshot with the following command:
longshot --force_overwrite --bam sample.genomealigned.bam --ref genome.fa --out sample.genomealigned.longshot.vcf --min_cov 3 --min_alt_count 3 --strand_bias_pvalue_cutoff 0.000001
You can use other variant calling tools or even variants called from WGS though.
You will also need to have run FLAIR quantify with the --output_bam option
so you have files of each sample aligned to the transcriptome.

**Output**

sample.isoforms.productivity.bed
    This is your isoforms with CDS annotation. Does not account for impact of variants.
sample.isovars.genomicpos.bed
    Genomic position of final set of variants
sample.isoswithvars.fa
    Sequences of variant-aware isoforms
sample.isoswithvars.counts.tsv
    Counts of variant-aware isoforms for each sample (large set, hard to do stats)
sample.aaseq.counts.tsv
    Counts of amino acid sequences for each sample (compact set, great for stats)
sample.aaseq.key.tsv
    Key of actual amino acid sequence associated with isoform/aaseq ID

Options
-------

.. code:: text

  -m --manifest
                        path to manifest files that points to sample files (see below). Each line of file
                        should be tab separated.
  -o --output_prefix
                        path to collapsed_output.bed file. default: 'flair'
  -i --isoforms
                        path to transcriptome fasta file
  -b --bedisoforms
                        path to transcriptome bed file
  -g --genome
                        FastA of reference genome
  -f --gtf              GTF annotation file

Manifest example:

Make sure bam files are from FLAIR quantify with --output_bam, 
not aligned to the genome

.. code:: text

   sample1      sample1.flair.aligned.bam      sample1.genomealigned.variants.vcf
   sample2      sample2.flair.aligned.bam      sample2.genomealigned.variants.vcf
   sample3      sample3.flair.aligned.bam      sample3.genomealigned.variants.vcf


.. _diffexp-label:

flair diffexp
=============


The standard `conda` environment no long installed `R` and the required packages.
These maybe added do the environment as describe in :ref:`installing-label` 

.. code:: text

   usage: flair diffexp -q counts_matrix.tsv --out_dir out_dir [options]


This module performs differential *expression* and differential *usage* analyses between **exactly two** conditions with 
3 or more replicates. Please have your control condition name (from the flair quantify manifest file) be alphabetically lower than your test condition for best results (eg ctl and test = good, untreated and treated = less good). It does so by running these R packages:

 - `DESeq2 <https://bioconductor.org/packages/release/bioc/html/DESeq2.html>`__ on genes and isoforms. This tests for differential expression.
 - `DRIMSeq <http://bioconductor.org/packages/release/bioc/html/DRIMSeq.html>`__ is used on isoforms only and tests for differential usage. This is done by testing if the ratio of isoforms changes between conditions.

If you do not have replicates you can use the `diff_iso_usage <#diffisoscript>`__ standalone script.

If you have more than two sample condtions, either split your counts matrix ahead of time or run DESeq2 and DRIMSeq yourself. 

**Outputs**

After the run, the output directory (``--out_dir``) contains the following, where COND1 and COND2 are the names of the sample groups.

 - ``genes_deseq2_MCF7_v_A549.tsv`` Filtered differential gene expression table.
 - ``genes_deseq2_QCplots_MCF7_v_A549.pdf`` QC plots, see the `DESeq2 manual <https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html>`__ for details.
 - ``isoforms_deseq2_MCF7_v_A549.tsv`` Filtered differential isoform expression table.
 - ``isoforms_deseq2_QCplots_MCF7_v_A549.pdf`` QC plots
 - ``isoforms_drimseq_MCF7_v_A549.tsv`` Filtered differential isoform usage table
 - ``workdir`` Temporary files including unfiltered output files.


Options
-------

Required arguments
~~~~~~~~~~~~~~~~~~

.. code:: text
    
    --counts_matrix	Tab-delimited isoform count matrix from flair quantify
    --out_dir	        Output directory for tables and plots.
    
Optional arguments
~~~~~~~~~~~~~~~~~~

.. code:: text
    
    --help	        Show this help message and exit
    --threads	        Number of threads for parallel DRIMSeq.
    --exp_thresh	Read count expression threshold. Isoforms in which both 
                        conditions contain fewer than E reads are filtered out (Default E=10) 
                        (This option requires that all replicates in either condition have > exp_thresh reads)
    --out_dir_force	Specify this argument to force overwriting of files in 
                        an existing output directory


Notes
-----

DESeq2 and DRIMSeq are optimized for short read experiments and expect many reads for each expressed gene. Lower coverage (as expected when using long reads) will tend to result in false positives.

For instance, look at this counts table with two groups (s and v) of three samples each:

.. code:: text

    gene   s1    s2      s3      v1      v2      v3
       A    1     0       2       0       4       2
       B  100    99     101     100     104     102

Gene A has an average expression of 1 in group s, and 2 in group v but the total variation in read count is 0-4. The same variation is true for gene B, but it will not be considered differentially expressed.

Flair does not remove low count genes as long as they are expressed in all samples of at least one group so please be careful when interpreting results.

Results tables are filtered and reordered by p-value so that only p<0.05 differential genes/isoforms remain. Unfiltered tables can be found in ``workdir``

.. _diffsplice-label:

flair diffsplice
================

The standard `conda` environment no long installed `R` and the required packages.
These maybe added do the environment as describe in :ref:`installing-label` 

.. code:: text

   usage: flair diffsplice -i isoforms.bed -q counts_matrix.tsv [options]

This module calls alternative splicing (AS) events from isoforms. Currently supports
the following AS events: 

 - intron retention (ir)
 - alternative 3’ splicing (alt3)
 - alternative 5’ splicing (alt5)
 - cassette exons (es)

If there are 3 or more samples per condition, then you can run with
``--test`` and DRIMSeq will be used to calculate differential usage of
the alternative splicing events between two conditions. See below for
more DRIMSeq-specific arguments. 

If conditions were sequenced without replicates, then the diffSplice output files can
be input to the `diffsplice_fishers_exact <#diffsplice_fishers>`__
script for statistical testing instead.

**Outputs**

After the run, the output directory (``--out_dir``) contains the following tab separated files:

 - ``diffsplice.alt3.events.quant.tsv``
 - ``diffsplice.alt5.events.quant.tsv``
 - ``diffsplice.es.events.quant.tsv``
 - ``diffsplice.ir.events.quant.tsv``

If DRIMSeq was run (where ``A`` and ``B`` are conditionA and conditionB, see below):

 - ``drimseq_alt3_A_v_B.tsv``
 - ``drimseq_alt5_A_v_B.tsv``
 - ``drimseq_es_A_v_B.tsv``
 - ``drimseq_ir_A_v_B.tsv``
 - ``workdir`` Temporary files including unfiltered output files.

Options
-------

Required arguments
~~~~~~~~~~~~~~~~~~

.. code:: text

    --isoforms	        Isoforms in bed format from Flair collapse.
    --counts_matrix	Tab-delimited isoform count matrix from Flair quantify.
    --out_dir	        Output directory for tables and plots.
    
Optional arguments
~~~~~~~~~~~~~~~~~~

.. code:: text
    
    --help	        Show all options.
    --threads	        Number of processors to use (default 4).
    --test	        Run DRIMSeq statistical testing.
    --drim1	        The minimum number of samples that have coverage over an 
                        AS event inclusion/exclusion for DRIMSeq testing; events 
                        with too few samples are filtered out and not tested (6).
    --drim2	        The minimum number of samples expressing the inclusion of 
                        an AS event; events with too few samples are filtered out 
                        and not tested (3).
    --drim3	        The minimum number of reads covering an AS event 
                        inclusion/exclusion for DRIMSeq testing, events with too 
                        few samples are filtered out and not tested (15).
    --drim4	        The minimum number of reads covering an AS event inclusion 
                        for DRIMSeq testing, events with too few samples are 
                        filtered out and not tested (5).
    --batch	        If specified with --test, DRIMSeq will perform batch correction.
    --conditionA	Specify one condition corresponding to samples in the 
                        counts_matrix to be compared against condition2; by default, 
                        the first two unique conditions are used. This implies --test.
    --conditionB	Specify another condition corresponding to samples in the 
                        counts_matrix to be compared against conditionA.
    --out_dir_force	Specify this argument to force overwriting of files in an 
                        existing output directory

Notes
-----

Results tables are filtered and reordered by p-value so that only p<0.05 differential genes/isoforms remain. Unfiltered tables can be found in ``workdir``

For a complex splicing example, please note the 2 alternative 3’ SS, 3
intron retention, and 4 exon skipping events in the following set of
isoforms that ``flair diffSplice`` would call and the isoforms that are
considered to include or exclude the each event:

.. figure:: img/toy_isoforms_coord.svg

.. code::

   a3ss_feature_id     coordinate                  sample1 sample2 ... isoform_ids
   inclusion_chr1:80   chr1:80-400_chr1:80-450     75.0    35.0    ... a,e
   exclusion_chr1:80   chr1:80-400_chr1:80-450     3.0     13.0    ... c
   inclusion_chr1:500  chr1:500-650_chr1:500-700   4.0     18.0    ... d
   exclusion_chr1:500  chr1:500-650_chr1:500-700   70.0    17.0    ... e

.. code::

   ir_feature_id           coordinate      sample1 sample2 ... isoform_ids
   inclusion_chr1:500-650  chr1:500-650    46.0    13.0    ... g
   exclusion_chr1:500-650  chr1:500-650    4.0     18.0    ... d
   inclusion_chr1:500-700  chr1:500-700    46.0    13.0    ... g
   exclusion_chr1:500-700  chr1:500-700    70.0    17.0    ... e
   inclusion_chr1:250-450  chr1:250-450    50.0    31.0    ... d,g
   exclusion_chr1:250-450  chr1:250-450    80.0    17.0    ... b

.. code::

   es_feature_id           coordinate      sample1 sample2 ... isoform_ids
   inclusion_chr1:450-500  chr1:450-500    83.0    30.0    ... b,c
   exclusion_chr1:450-500  chr1:450-500    56.0    15.0    ... f
   inclusion_chr1:200-250  chr1:200-250    80.0    17.0    ... b
   exclusion_chr1:200-250  chr1:200-250    3.0     13.0    ... c
   inclusion_chr1:200-500  chr1:200-500    4.0     18.0    ... d
   exclusion_chr1:200-500  chr1:200-500    22.0    15.0    ... h
   inclusion_chr1:400-500  chr1:400-500    75.0    35.0    ... e,a
   exclusion_chr1:400-500  chr1:400-500    56.0    15.0    ... f

