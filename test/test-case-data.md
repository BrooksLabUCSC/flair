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
reads.  Leads to a fast test that the code paths just work.

## set seg1 - single exon genes
ONT single-exon gene tests on GRCh38

Cases:
* chr20:32188017-32190331 ENSG00000293164 single-exon protein-coding gene

Files:
* seg1.cdna-ont.fastq -  LRGASP WTC11 cDNA ONT reads
* seg1.gencodeV47.gtf - from GENCODE V47, as some test cases added from LRGASP
* seg1.promoter-regions.bed - from encodeCcreCombined
