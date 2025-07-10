library(DESeq2)
library(ggplot2)
library(qqman)
library(argparse)
options(error = function() traceback(2))

# Function for parsing command line arguments
parse_arguments <- function() {
  parser <- ArgumentParser(description = 'run DESeq2 for flair_diffExp')

  parser$add_argument("--group1", required = TRUE, help = "Sample group 1.")
  parser$add_argument("--group2", required = TRUE, help = "Sample group 2.")
  parser$add_argument("--batch", required = FALSE, default = NULL, help = "Secondary sample attribute (used in design matrix).")
  parser$add_argument("--matrix", required = TRUE, help = "Input count files.")
  parser$add_argument("--outDir", required = TRUE, help = "Write to specified output directory.")
  parser$add_argument("--prefix", required = TRUE, help = "Specify file prefix.")
  parser$add_argument("--formula", required = TRUE, help = "Formula design matrix.")

  args <- parser$parse_args()
  return(args)
}

# Function to run DESeq2 and save results
run_deseq_analysis <- function(args) {
  outdir <- args$outDir
  group1 <- args$group1
  group2 <- args$group2
  prefix <- args$prefix
  matrixFile <- args$matrix
  formulaFile <- args$formula
  
  data_folder <- normalizePath(outdir, mustWork = FALSE)
  workdir <- file.path(data_folder, 'workdir')
  if (!dir.exists(workdir)) {
    dir.create(workdir, recursive = TRUE)
  }
  
  lfcOut <- file.path(workdir, sprintf("%s_%s_v_%s_results_shrinkage.tsv", prefix, group1, group2))
  resOut <- file.path(workdir, sprintf("%s_%s_v_%s_results.tsv", prefix, group1, group2))
  cleanOut <- file.path(data_folder, sprintf("%s_%s_v_%s.tsv", prefix, group1, group2))
  
  countData <- read.table(matrixFile, header = TRUE, sep = "\t", row.names = 1)
  colData <- read.table(formulaFile, header = TRUE, sep = "\t", row.names = 1)
  
  design <- if ("batch" %in% colnames(colData)) {
    ~ condition + condition
  } else {
    ~ condition
  }
  
  dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = design)
  dds <- DESeq(dds)
  name <- grep("condition", resultsNames(dds), value=TRUE)
  
  res <- results(dds, name = name)
  resLFC <- lfcShrink(dds, coef = name)
  write.table(as.data.frame(res), file = resOut, quote = FALSE, sep = "\t")
  write.table(as.data.frame(resLFC), file = lfcOut, quote = FALSE, sep = "\t")
  
  resLFC <- na.omit(resLFC[order(resLFC$padj), ])
  resLFC <- as.data.frame(resLFC[resLFC$padj < 0.05, ])
  outdf <- data.frame(sample = rownames(resLFC), round(resLFC[, 1:3], 2), signif(resLFC[, 4:5], 3))
  write.table(outdf, file = cleanOut, row.names = FALSE, quote = FALSE, sep = "\t")
  
  return(dds)
}

# Function for plotting results
plot_results <- function(dds, args) {
  group1 <- args$group1
  group2 <- args$group2
  prefix <- args$prefix
  outdir <- args$outDir
  matrixFile <- args$matrix

  data_folder <- normalizePath(outdir, mustWork = FALSE)
  qcOut <- file.path(data_folder, sprintf("%s_QCplots_%s_v_%s.pdf", prefix, group1, group2))
  
  pdf(qcOut)

  plotMA(results(dds), ylim = c(-3, 3), main = "MA-plot results")
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
