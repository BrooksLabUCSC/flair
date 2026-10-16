library(DESeq2)
library(ggplot2)
library(qqman)
library(argparse)
# traceback alone leaves the exit status at zero, so a failed run looks like a
# successful one to whatever called this script
options(error = function() {traceback(2); quit(save="no", status=1)})
source(file.path(dirname(sub("^--file=", "",
                            grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))),
                 "tsv_io.R"))

# Function for parsing command line arguments
parse_arguments <- function() {
  parser <- ArgumentParser(description = 'run DESeq2 for flair_diffExp')

  parser$add_argument("--condition_a", required = TRUE, help = "Reference condition; fold changes are relative to this.")
  parser$add_argument("--condition_b", required = TRUE, help = "Condition compared against condition_a.")
  parser$add_argument("--matrix", required = TRUE, help = "Input count files.")
  parser$add_argument("--out_dir", required = TRUE, help = "Write to specified output directory.")
  parser$add_argument("--prefix", required = TRUE, help = "Specify file prefix.")
  parser$add_argument("--formula", required = TRUE, help = "Formula design matrix.")
  parser$add_argument("--compress", action = "store_true", default = FALSE,
                      help = "Compress the output tables with gzip.")

  args <- parser$parse_args()
  return(args)
}

# Function to run DESeq2 and save results
run_deseq_analysis <- function(args) {
  outdir <- args$out_dir
  condition_a <- args$condition_a
  condition_b <- args$condition_b
  prefix <- args$prefix
  matrixFile <- args$matrix
  formulaFile <- args$formula
  
  data_folder <- normalizePath(outdir, mustWork = FALSE)
  workdir <- file.path(data_folder, 'workdir')
  if (!dir.exists(workdir)) {
    dir.create(workdir, recursive = TRUE)
  }
  
  lfcOut <- tsvCompressedPath(file.path(workdir, sprintf("%s_%s_v_%s_results_shrinkage.tsv", prefix, condition_a, condition_b)), args$compress)
  resOut <- tsvCompressedPath(file.path(workdir, sprintf("%s_%s_v_%s_results.tsv", prefix, condition_a, condition_b)), args$compress)
  cleanOut <- tsvCompressedPath(file.path(data_folder, sprintf("%s_%s_v_%s.tsv", prefix, condition_a, condition_b)), args$compress)
  
  countData <- read.table(matrixFile, header = TRUE, sep = "\t", row.names = 1)
  colData <- read.table(formulaFile, header = TRUE, sep = "\t", row.names = 1)
  
  design <- if ("batch" %in% colnames(colData)) {
    ~ condition + batch
  } else {
    ~ condition
  }
  
  dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = design)
  dds$condition <- relevel(dds$condition, ref=condition_a)
  dds <- DESeq(dds)
  name <- paste('condition_', condition_b, '_vs_', condition_a, sep='')

  res <- results(dds, name = name)
  resLFC <- lfcShrink(dds, coef = name)
  tsvWriteTable(as.data.frame(res), resOut)
  tsvWriteTable(as.data.frame(resLFC), lfcOut)
  
  resLFC <- na.omit(resLFC[order(resLFC$padj), ])
  resLFC <- as.data.frame(resLFC[resLFC$padj < 0.05, ])
  outdf <- data.frame(sample = rownames(resLFC), round(resLFC[, 1:3], 2), signif(resLFC[, 4:5], 3))
  tsvWriteTable(outdf, cleanOut, row.names = FALSE)
  
  return(dds)
}

# Function for plotting results
plot_results <- function(dds, args) {
  condition_a <- args$condition_a
  condition_b <- args$condition_b
  prefix <- args$prefix
  outdir <- args$out_dir
  matrixFile <- args$matrix

  data_folder <- normalizePath(outdir, mustWork = FALSE)
  qcOut <- file.path(data_folder, sprintf("%s_QCplots_%s_v_%s.pdf", prefix, condition_a, condition_b))
  
  pdf(qcOut)

  # the named coefficient, as the results table uses: bare results(dds) takes the last
  # coefficient in resultsNames, which is the batch term once batch is in the design
  name <- paste('condition_', condition_b, '_vs_', condition_a, sep='')
  plotMA(results(dds, name = name), ylim = c(-3, 3),
         main = sprintf("MA-plot: %s vs %s", condition_b, condition_a))
  plotDispEsts(dds, main = "Dispersion Estimates")

  nsub <- min(nrow(read.table(matrixFile, header = TRUE, sep = "\t", row.names = 1)), 1000)
  vsd <- tryCatch({
    vst(dds, nsub = nsub, blind = FALSE)
  }, error = function(e) {
    dev.off()
    stop('DESeq2 and other QC plots ran OK but the PCA plot failed, probably because the number of input genes is very low')
  })

  colData <- read.table(args$formula, header = TRUE, sep = "\t", row.names = 1)
  pcaData <- plotPCA(vsd, intgroup = if ("batch" %in% colnames(colData)) c("condition", "batch") else "condition", returnData = TRUE)
  percentVar <- attr(pcaData, "percentVar")
  
  x_label <- sprintf("PC1: %.1f%% variance", percentVar[1] * 100)
  y_label <- sprintf("PC2: %.1f%% variance", percentVar[2] * 100)
  
  plot <- ggplot(pcaData, aes(x = PC1, y = PC2, color = condition)) +
    geom_point(size = 3) +
    xlab(x_label) +
    ylab(y_label) +
    theme_classic()

  if ("batch" %in% colnames(colData)) {
    plot <- plot + aes(shape = batch)
  }
  
  print(plot)
  
  dev.off()
}

# Main function to handle the workflow
main <- function() {
  args <- parse_arguments()
  dds <- run_deseq_analysis(args)
  plot_results(dds, args)
}

main()
