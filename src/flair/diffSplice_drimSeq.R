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
  parser$add_argument("--outDir", default='', help="Write to specified output directory.")
  parser$add_argument("--prefix", required=TRUE, help="Specify file prefix.")
  parser$add_argument('--drim1', type="integer", default=6, help="Minimum number of samples expressing event inclusion/exclusion (6).")
  parser$add_argument('--drim2', type="integer", default=3, help="Minimum number of samples expressing the inclusion of an event (3).")
  parser$add_argument('--drim3', type="integer", default=15, help="Minimum number of reads covering an event inclusion/exclusion (15).")
  parser$add_argument('--drim4', type="integer", default=5, help="Minimum number of reads covering an event inclusion (5).")
  parser$add_argument("--threads", type="integer", default=4, help="Number of threads for running DRIM-Seq.")
  parser$add_argument('--batch', action='store_true', default=FALSE, help="If specified, batch correction will be performed.")
  parser$add_argument('--conditionA', default='', help="Specify one condition to compare against conditionB.")
  parser$add_argument('--conditionB', default='', help="Specify one condition to compare against conditionA.")
  return(parser$parse_args())
}

# Function to run DRIMSeq
run_DRIMSeq <- function(args) {
  cat('Input file:', args$matrix, "\n", file=stderr())

  # Create output working directory if it doesn't exist
  workdir <- file.path(args$outDir, 'workdir')
  if (!dir.exists(workdir)) {
    dir.create(workdir, recursive=TRUE)
  }

  # Read sample info from the matrix
  sample_info <- fread(args$matrix, nrows=1)
  samples <- colnames(sample_info)[-c(1:2, length(sample_info))]
  groups <- sapply(strsplit(samples, "_"), `[`, 2)
  batches <- sapply(strsplit(samples, "_"), `[`, 3)

  # Determine conditionA and conditionB if not provided
  if (args$conditionA == '') {
    args$conditionA <- groups[1]
    args$conditionB <- groups[which(groups != args$conditionA)[1]]
  }

  # Fill the formula data.frame
  formulaDF <- data.table(sample_id = character(), condition = character(), batch = character())
  for (i in seq_along(samples)) {
    if (args$conditionA %in% groups[i] || args$conditionB %in% groups[i]) {
      formulaDF <- rbind(formulaDF, list(samples[i], groups[i], batches[i]))
    }
  }

  if (nrow(formulaDF) == 0) {
    cat(sprintf('\n**ERROR** Could not find %s and/or %s in input file, exiting\n\n', args$conditionA, args$conditionB), file=stderr())
    stop(sprintf('Could not find %s and/or %s in input file', args$conditionA, args$conditionB))
  }

  # Read counts and prepare the quantification data frame
  quantDF <- fread(args$matrix)
  quantDF <- quantDF[, c("feature_id", "coordinate", samples, "isoform_ids"), with=FALSE]
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
  filtered <- dmFilter(data, min_samps_gene_expr = args$drim1, min_samps_feature_expr = args$drim2, min_gene_expr = args$drim3, min_feature_expr = args$drim4)

  # Design matrix and fitting
  design_full <- if (args$batch) {
    model.matrix(~ condition + batch, data = samples(filtered))
  } else {
    model.matrix(~ condition, data = samples(filtered))
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
  res_out_path <- file.path(workdir, sprintf("%s_%s_v_%s_drimseq_results.tsv", args$prefix, args$conditionA, args$conditionB))
  write.table(res, file=res_out_path, quote=FALSE, sep='\t', row.names=FALSE)

  # Finalize and write cleaned results
  clean_out_path <- file.path(args$outDir, sprintf("drimseq_%s_%s_v_%s.tsv", args$prefix, args$conditionA, args$conditionB))
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
