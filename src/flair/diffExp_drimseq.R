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
  parser$add_argument("--matrix", required=TRUE, help='Input count files.')
  parser$add_argument("--out_dir", required=TRUE, help='Write to specified output directory.')
  parser$add_argument("--prefix", required=TRUE, help='Specify file prefix.')
  parser$add_argument("--formula", required=TRUE, help='Formula design matrix.')
  parser$add_argument("--threads", type="integer", default=4, help='Number of threads for running DRIM-Seq. BBPARAM')
  # the same four dmFilter thresholds diffSplice_drimSeq.R exposes; they were inline
  # here, so the same analysis was tunable in one module and fixed in the other
  parser$add_argument('--min_samps_gene_expr', type="integer", default=6, help='Minimum number of samples expressing a gene (6).')
  parser$add_argument('--min_samps_feature_expr', type="integer", default=3, help='Minimum number of samples expressing a feature (3).')
  parser$add_argument('--min_gene_expr', type="integer", default=15, help='Minimum number of reads covering a gene (15).')
  parser$add_argument('--min_feature_expr', type="integer", default=5, help='Minimum number of reads covering a feature (5).')
  
  args <- parser$parse_args()
  return(args)
}

main <- function() {
  args <- parse_args()
  
  outdir <- args$out_dir
  group1 <- args$group1
  group2 <- args$group2
  matrix <- args$matrix
  prefix <- args$prefix
  formula <- args$formula
  threads <- args$threads
  
  rundrimseq(outdir, group1, group2, matrix, prefix, formula, threads,
             args$min_samps_gene_expr, args$min_samps_feature_expr, args$min_gene_expr, args$min_feature_expr)
}

rundrimseq <- function(outdir, group1, group2, matrix, prefix, formula, threads,
                       min_samps_gene_expr, min_samps_feature_expr, min_gene_expr, min_feature_expr) {
  cat(sprintf('input file: %s\n', matrix), file=stderr())
  
  # create output working directory if it doesn't exist
  data_folder <- normalizePath(outdir, mustWork = FALSE)
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
  
  
  filtered <- dmFilter(data, min_samps_gene_expr = min_samps_gene_expr, min_samps_feature_expr = min_samps_feature_expr,
                       min_gene_expr = min_gene_expr, min_feature_expr = min_feature_expr)
  
  # group1 is the reference, so the reported fold change has the direction the output
  # file name states.  Without this, condition is a character column and model.matrix
  # orders its levels alphabetically
  filtered_samples <- samples(filtered)
  filtered_samples$condition <- relevel(factor(filtered_samples$condition), ref=group1)

  if ("batch" %in% names(formulaDF)) {
    design_full <- model.matrix(~ condition + batch, data = filtered_samples)
  } else {
    design_full <- model.matrix(~ condition, data = filtered_samples)
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
