#!/usr/bin/env Rscript
# One DESeq2 contrast.   Rscript deseq2_contrast.R <counts.tsv> <coldata.tsv> <mode> <out.tsv> [<lengths.tsv>]
#   mode "dose"     : design ~ dose (numeric rank 0/1/2); Wald on the slope; apeglm-shrunk slope; log2fc = 2 * shrunk slope
#                     (= top dose vs control); pvalue/padj from the slope test.
#   mode "pairwise" : design ~ condition (2 levels, control first); Wald on treated vs control; apeglm-shrunk log2fc.
#   lengths.tsv     : optional gene x well matrix of gene lengths (any positive unit; 1 for wells without a length term; de.py
#                     passes one length per gene per method, so within a method the offset cancels exactly).
#                     Stored as assay "avgTxLength", from which DESeq() builds per-gene normalisation factors exactly as it does
#                     for tximport input (DESeqDataSetFromTximport): nf = size factor x length / geometric mean of the gene's
#                     lengths across wells -- the length enters the model as an offset, the counts stay counts; apeglm takes
#                     the same offsets.
# Output columns: gene_id baseMean log2fc (shrunk effect) log2fc_raw lfcSE stat pvalue padj + log2(normalised count + 1) per well.
# The stdout line reports the size factors, or with lengths the per-well geometric mean of the normalisation factors.
suppressPackageStartupMessages({library(DESeq2)})
a <- commandArgs(trailingOnly = TRUE); cts <- as.matrix(read.delim(a[1], row.names = 1, check.names = FALSE)); cd <- read.delim(a[2], stringsAsFactors = FALSE)
rownames(cd) <- cd$well; cd <- cd[colnames(cts), ]; cts <- cts[rowSums(cts) >= 10, ]
len <- if (length(a) >= 5) as.matrix(read.delim(a[5], row.names = 1, check.names = FALSE))[rownames(cts), colnames(cts)] else NULL
if (a[3] == "dose") {
  dds <- DESeqDataSetFromMatrix(cts, cd, ~ dose); coef <- "dose"; mult <- 2
} else {
  cd$condition <- factor(cd$condition, levels = unique(cd$condition[order(cd$dose)]))
  dds <- DESeqDataSetFromMatrix(cts, cd, ~ condition); mult <- 1
}
if (!is.null(len)) { stopifnot(all(is.finite(len)), all(len > 0)); dimnames(len) <- dimnames(cts); assays(dds)[["avgTxLength"]] <- len }
dds <- DESeq(dds, quiet = TRUE); if (a[3] != "dose") coef <- resultsNames(dds)[2]
res <- results(dds, name = coef); shr <- lfcShrink(dds, coef = coef, type = "apeglm", quiet = TRUE)
out <- data.frame(gene_id = rownames(res), baseMean = res$baseMean, log2fc = mult * shr$log2FoldChange, log2fc_raw = mult * res$log2FoldChange,
                  lfcSE = res$lfcSE, stat = res$stat, pvalue = res$pvalue, padj = res$padj)
nc <- log2(counts(dds, normalized = TRUE) + 1); colnames(nc) <- paste0("log2norm_", colnames(nc))
write.table(cbind(out, nc[rownames(res), ]), a[4], sep = "\t", quote = FALSE, row.names = FALSE)
sf <- if (is.null(len)) sizeFactors(dds) else exp(colMeans(log(normalizationFactors(dds))))
cat(sprintf("%s: %d genes; coef %s;%s size factors %s\n", a[3], nrow(cts), coef, if (is.null(len)) "" else " length offsets (avgTxLength -> per-gene normalisation factors, reported as their geometric mean per well);", paste(sprintf("%.2f", sf), collapse = " ")))
