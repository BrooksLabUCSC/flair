flair diffExp
=============

.. code:: sh

   usage: python flair.py diffExp -q counts_matrix.tsv --out_dir out_dir [options]


This module performs differential *expression* and differential *usage* analyses between **exactly two** conditions with 
3 or more replicates. It does so by running these R packages:

 - `DESeq2 <https://bioconductor.org/packages/release/bioc/html/DESeq2.html>`__ on genes and isoforms. This tests for differential expression.
 - `DRIMSeq <http://bioconductor.org/packages/release/bioc/html/DRIMSeq.html>`__ is used on isoforms only and tests for differential usage. This is done by testing if the ratio of isoforms changes between conditions.

If you do not have replicates you can use the `diff_iso_usage <#diffisoscript>`__ standalone script.

If you have more than two sample condtions, either split your counts matrix ahead of time or run DESeq2 and DRIMSeq yourself. 

Outputs

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
``--counts_matrix`` Tab-delimited isoform count matrix from flair quantify module.

``--out_dir`` Output directory for tables and plots.

Optional arguments
~~~~~~~~~~~~~~~~~~
``--help`` show this help message and exit

``--threads`` Number of threads for parallel DRIMSeq.

``--exp_thresh`` Read count expression threshold. Isoforms in which both conditions contain fewer than E reads are filtered out (Default E=10)

``--out_dir_force`` Specify this argument to force overwriting of files in an existing output directory



Notes
-----

DESeq2 and DRIMSeq are optimized for short read experiments and expect many reads for each expressed gene. Lower coverage (as expected when using long reads) will tend to result in false positives.

For instance, look at this counts table with two groups (s and v) of three samples each:

.. code:: sh

    gene   s1    s2      s3      v1      v2      v3
       A    1     0       2       0       4       2
       B  100    99     101     100     104     102

Gene A has an average expression of 1 in group s, and 2 in group v but the total variation in read count is 0-4. The same variation is true for gene B, but it will not be considered differentially expressed.

Flair does not remove low count genes as long as they are expressed in all samples of at least one group so please be careful when interpreting results.

Results tables are filtered and reordered by p-value so that only p<0.05 differential genes/isoforms remain. Unfiltered tables can be found in ``workdir``

Code requirements
~~~~~~~~~~~~~~~~~
This module requires python modules and R packages that are not necessary for other Flair modules (except diffSplice).  

**If you are not using the docker container or the conda installed version of Flair** you may have to install these separately:

1. python modules: pandas, numpy, rpy2
2. `DESeq2 <https://bioconductor.org/packages/release/bioc/html/DESeq2.html>`__
3. `ggplot2 <https://ggplot2.tidyverse.org>`__
4. `qqman <https://cran.r-project.org/web/packages/qqman/index.html>`__
5. `DRIMSeq <http://bioconductor.org/packages/release/bioc/html/DRIMSeq.html>`__
6. `stageR <http://bioconductor.org/packages/release/bioc/html/stageR.html>`__


