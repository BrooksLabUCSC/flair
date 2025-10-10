# Test case development

One has to use great care in selecting test cases so to avoid create
large data files checked into git and slow tests.

Data files names should be in the form:
`input/<testset_name>.<something_useful>.<ext>`.

For a given test data set, the results from a previous modules in the
`expected/` directory are used as input to following steps.  While this makes
tests interdependent, it ensure that downstream tests are using the latest
output from preceding modules.

## basic - basic test cases
These are LRGASP WTC-11 cDNA short reads pretending to be long 
reads.  
WTC11 PacBio cDNA from LRGASP 
Sample: ENCFF245IPA
For the three genes: KRAS, ERGIC3, DDX5
Took a random subset of one third of the reads aligned to each locus, this is the core set for collapse

For quantify, from the third subset of IPA used for collapse:
took 3 random subsets at 1/10th of that file: reps 1-3
took 3 random subsets at 1/20th of that file: reps 4-6

Leads to a fast test that the code paths just work.

## set seg1 - single exon genes
LRGASP WTC11 ONT single-exon gene tests on GRCh38

Cases:
* chr20:32188017-32190331 ENSG00000293164 single-exon protein-coding gene

Files:
* seg1.cdna-ont.fastq -  LRGASP WTC11 cDNA ONT reads
* seg1.gencodeV47.gtf - from GENCODE V47, as some test cases added from LRGASP
* seg1.promoter-regions.bed - from encodeCcreCombined

## tardigrade-regress
A test case derived from tardigrade ONT RNA data that caused flair correct to
generate invalid BAMs.  This was done on a pre-release genome so only the
alignment SAM and BED, and the introns for the region are check in.  When the
genome is release, this can be turned into a full test case.

## Others
* tiny.sam - A few records from test-align.bam output
