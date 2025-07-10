library(argparse)
library(methods)
library(DRIMSeq)
library(BiocParallel)
options(error = function() traceback(2))

# CommandLine
parse_args <- function() {
  parser <- ArgumentParser(description='run DRIMSeq for flair_diffExp')
  
  parser$add_argument("--group1", required=TRUE, help='Sample group 1.')
  parser$add_argument("--group2", required=TRUE, help='Sample group 2.')
  parser$add_argument("--batch", required=FALSE, default=NULL, help='Secondary sample attribute (used in design matrix).')
  parser$add_argument("--matrix", required=TRUE, help='Input count files.')
  parser$add_argument("--outDir", required=TRUE, help='Write to specified output directory.')
  parser$add_argument("--prefix", required=TRUE, help='Specify file prefix.')
  parser$add_argument("--formula", required=TRUE, help='Formula design matrix.')
  parser$add_argument("--threads", type="integer", default=4, help='Number of threads for running DRIM-Seq. BBPARAM')
  
  args <- parser$parse_args()
  return(args)
}

main <- function() {
  args <- parse_args()
  
  outdir <- args$outDir
  group1 <- args$group1
  group2 <- args$group2
  matrix <- args$matrix
  prefix <- args$prefix
  formula <- args$formula
  threads <- args$threads
  
  rundrimseq(outdir, group1, group2, matrix, prefix, formula, threads)
}

rundrimseq <- function(outdir, group1, group2, matrix, prefix, formula, threads) {
  cat(sprintf('input file: %s\n', matrix), file=stderr())
  
  # create output working directory if it doesn't exist
  data_folder <- file.path(getwd(), outdir)
  workdir <- file.path(data_folder, 'workdir')
  if (!dir.exists(workdir)) {
    dir.create(workdir, recursive=TRUE)
  }
  
  resOut <- file.path(workdir, sprintf("%s_%s_v_%s_results.tsv", prefix, group1, group2))
  cleanOut <- file.path(data_folder, sprintf("%s_%s_v_%s.tsv", prefix, group1, group2))
  
  # Import data
  quantDF <- read.table(matrix, header=TRUE, sep='\t', row.names=1, check.names=FALSE)
  formulaDF <- read.table(formula, header=TRUE, sep="\t", check.names=FALSE)
  
  samples <- formulaDF
  data <- dmDSdata(counts = quantDF, samples = samples)
  
  # DRIMSEQ part
  if ("batch" %in% names(formulaDF)) {
    batch <- samples$batch
  }
  condition <- samples$condition
  
  
  filtered <- dmFilter(data, min_samps_gene_expr = 6, min_samps_feature_expr = 3, min_gene_expr = 15, min_feature_expr = 5)
  
  if ("batch" %in% names(formulaDF)) {
    design_full <- model.matrix(~ condition + batch, data = samples(filtered))
  } else {
    design_full <- model.matrix(~ condition, data = samples(filtered))
  }
  
  set.seed(123)

  d <- dmPrecision(filtered, design = design_full, BPPARAM=MulticoreParam(threads))
  d <- dmFit(d, design = design_full, verbose = 1, BPPARAM=MulticoreParam(threads))
  
  contrast <- colnames(design_full)[2]
  
  d <- dmTest(d, coef = contrast, verbose = 1, BPPARAM=MulticoreParam(threads))
  res <- merge(proportions(d), results(d, level="feature"), by=c("feature_id","gene_id"))
  
  # Write raw output
  write.table(res, file=resOut, sep="\t", quote=FALSE, row.names=FALSE)

  # Order by adjusted p-value and keep only significant
  res <- res[order(res[,"adj_pvalue"]),]
  res <- subset(res, adj_pvalue <= 0.05)
  
  # Final data preparation for output
  fina <- ncol(res) - 4
  outdf <- data.frame(res[1:2], round(res[3:fina],3), lr = round(res[["lr"]],2), adj_pvalue = signif(res[["adj_pvalue"]],3))
  
  write.table(outdf, file=cleanOut, row.names=FALSE, quote=FALSE, sep="\t")
}

main()
