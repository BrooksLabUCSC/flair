Example Files 
=============

We have provided the following example files `here <https://people.ucsc.edu/~atang14/flair/example_files/>`__:

* ``star.firstpass.gm12878.junctions.3.tab``, a file of splice junctions observed from short read sequencing of GM18278 that can be used in the correction step with ``-j``. Junctions with fewer than 3 uniquely mapping reads have been filtered out.
* ``promoter.gencode.v27.20.bed``, promoter regions determined from `ENCODE promoter chromatin states for GM12878 <http://hgdownload.cse.ucsc.edu/goldenPath/hg18/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmmGm12878HMM.bed.gz>`_ and 20 bp around annotated TSS in GENCODE v27. Can be supplied to flair-collapse with ``-p`` to build the initial firstpass set with only reads with start positions falling within these regions

Other downloads:
 
* `Native RNA Pass reads <https://github.com/nanopore-wgs-consortium/NA12878/blob/master/RNA.md>`__ Running these 10 million nanopore reads from ``fastq`` through flair align, correct, and collapse modules to assembled isoforms with 8 threads requires ~3.5 hours (includes ~2.5 hours of minimap2 alignment) 
* `NanoSim_Wrapper.py <https://github.com/BrooksLabUCSC/labtools/blob/master/NanoSim_Wrapper.py>`__, a wrapper script written for simulating nanopore transcriptome data using `Nanosim <https://github.com/bcgsc/NanoSim>`__



