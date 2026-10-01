#!/usr/bin/env Rscript

library(argparse)
library(DRIMSeq)
library(data.table)
suppressWarnings(library(BiocParallel))
options(error = function() traceback(2))

# Function to parse command-line arguments
parse_arguments <- function() {
  parser <- ArgumentParser(description = 'run dirmSeq')
  parser$add_argument("--matrix", required=TRUE, help="Input DRIM-Seq formatted count files.")
  # required, not default='': file.path("", "workdir") is "/workdir", at the
  # filesystem root, so the default could never work
  parser$add_argument("--out_dir", required=TRUE, help="Write to specified output directory.")
  parser$add_argument("--prefix", required=TRUE, help="Specify file prefix.")
  parser$add_argument('--min_samps_gene_expr', type="integer", default=6, help="Minimum number of samples expressing event inclusion/exclusion (6).")
  parser$add_argument('--min_samps_feature_expr', type="integer", default=3, help="Minimum number of samples expressing the inclusion of an event (3).")
  parser$add_argument('--min_gene_expr', type="integer", default=15, help="Minimum number of reads covering an event inclusion/exclusion (15).")
  parser$add_argument('--min_feature_expr', type="integer", default=5, help="Minimum number of reads covering an event inclusion (5).")
  parser$add_argument("--threads", type="integer", default=4, help="Number of threads for running DRIM-Seq.")
  parser$add_argument('--batch', action='store_true', default=FALSE, help="If specified, batch correction will be performed.")
  parser$add_argument('--condition_a', required=TRUE, help="Reference condition; the comparison is condition_b against this.")
  parser$add_argument('--condition_b', required=TRUE, help="Condition compared against condition_a.")
  parser$add_argument("--formula", required=TRUE, help="TSV of sample_id, condition and batch, one row per matrix sample column.")
  return(parser$parse_args())
}

# Function to run DRIMSeq
run_DRIMSeq <- function(args) {
  cat('Input file:', args$matrix, "\n", file=stderr())

  # Create output working directory if it doesn't exist
  workdir <- file.path(args$out_dir, 'workdir')
  if (!dir.exists(workdir)) {
    dir.create(workdir, recursive=TRUE)
  }

  # The condition and batch of each sample come from the formula file rather than
  # from the matrix column names, which carry only the sample id.
  formulaDF <- fread(args$formula, colClasses="character")
  formulaDF <- formulaDF[condition %in% c(args$condition_a, args$condition_b)]

  if (nrow(formulaDF) == 0) {
    cat(sprintf('\n**ERROR** Could not find %s and/or %s in input file, exiting\n\n', args$condition_a, args$condition_b), file=stderr())
    stop(sprintf('Could not find %s and/or %s in input file', args$condition_a, args$condition_b))
  }

  # Read counts and prepare the quantification data frame
  quantDF <- fread(args$matrix)
  # formulaDF$sample_id, not samples: the matrix used to keep every sample column
  # while the pseudocount loop below covered only the selected ones, so the counts
  # handed to dmDSdata were inflated in some columns and not others
  quantDF <- quantDF[, c("feature_id", "coordinate", formulaDF$sample_id, "isoform_ids"), with=FALSE]
  setnames(quantDF, "coordinate", "gene_id")

  # Add pseudocount
  for (sample in formulaDF$sample_id) {
    quantDF[, eval(sample) := .SD[[sample]] + 1, .SDcols = sample]
  }

  # Convert data.table to data.frame for DRIMSeq compatibility
  sample_df <- as.data.frame(formulaDF)
  count_df <- as.data.frame(quantDF)

  # Initialize and filter data
  data <- dmDSdata(counts = count_df, samples = sample_df)
  filtered <- dmFilter(data, min_samps_gene_expr = args$min_samps_gene_expr, min_samps_feature_expr = args$min_samps_feature_expr, min_gene_expr = args$min_gene_expr, min_feature_expr = args$min_feature_expr)

  # condition_a is the reference, so the reported fold change has the direction the
  # output file name states.  Without this, condition is a character column and
  # model.matrix orders its levels alphabetically
  filtered_samples <- samples(filtered)
  filtered_samples$condition <- relevel(factor(filtered_samples$condition), ref=args$condition_a)

  # Design matrix and fitting
  design_full <- if (args$batch) {
    model.matrix(~ condition + batch, data = filtered_samples)
  } else {
    model.matrix(~ condition, data = filtered_samples)
  }

  set.seed(123)
  BPPARAM <- MulticoreParam(workers = args$threads)
  d <- dmPrecision(filtered, design = design_full, BPPARAM = BPPARAM)
  d <- dmFit(d, design = design_full, verbose = 1, BPPARAM = BPPARAM)

  # Perform the differential test
  contrast <- grep("condition", colnames(design_full), value=TRUE)
  d <- dmTest(d, coef = contrast, verbose = 1, BPPARAM = BPPARAM)

  # Capture results
  res <- merge(proportions(d), results(d, level = "feature"), by=c("feature_id", "gene_id"))
  res_out_path <- file.path(workdir, sprintf("%s_%s_v_%s_drimseq_results.tsv", args$prefix, args$condition_a, args$condition_b))
  write.table(res, file=res_out_path, quote=FALSE, sep='\t', row.names=FALSE)

  # Finalize and write cleaned results
  clean_out_path <- file.path(args$out_dir, sprintf("drimseq_%s_%s_v_%s.tsv", args$prefix, args$condition_a, args$condition_b))
  res <- res[order(res$adj_pvalue),]
  res <- subset(res, adj_pvalue <= 0.05)
  final_idx <- ncol(res) - 4
  outdf <- data.frame(res[,1:2], round(res[,3:final_idx], 3), lr = round(res[["lr"]], 2), adj_pvalue = signif(res[["adj_pvalue"]], 3))
  write.table(outdf, file=clean_out_path, quote=FALSE, sep='\t', row.names=FALSE)

  return(clean_out_path)
}

# Main function to orchestrate tasks
main <- function() {
  args <- parse_arguments()
  outfile <- run_DRIMSeq(args)
}

main()
