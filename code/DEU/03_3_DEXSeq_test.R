OK <- requireNamespace("devtools", quietly = TRUE)
if (!OK) {
  stop("devtools package required but is not installed (or can't be loaded)")
}
# source("/vast/projects/lab_chen/tam/tools/Subread_to_DEXSeq/load_SubreadOutput_modified.R")
source("../../tools/Subread_to_DEXSeq/load_SubreadOutput_modified.R")
library(dplyr)
library(DEXSeq)
library(BiocParallel)

# Getting parameters
args <- R.utils::commandArgs(asValues = TRUE)
print(args)

DIR <- as.character(args[['DIR']])
flat_GTF <- as.character(args[['flat_GTF']])
libs <- as.integer(args[['libs']])
workers <- as.integer(args[['workers']])
target <- data.frame(row.names = paste0("Rep", 1:(libs*2)), 
                     condition = rep(c("Group_1", "Group_2"), each=libs))

i <- "S1"

### 0. Check parameter
message("Check valid parameters")
print(DIR)
print(flat_GTF)
print(target)
print(libs)
print(workers)

message(paste0("Processing: ", i, " ..."))
count_f <- paste0(DIR, "STAR_output/featureCounts_subread_fit_RL_75/", i, "/featureCounts_output.txt")
print(count_f)
DEXSeq_output <- paste0(DIR, "STAR_output/DEXSeq_fit_RL_75_test/", i)
print(DEXSeq_output)
dir.create(DEXSeq_output, showWarnings = FALSE, recursive = TRUE)

### Create DEXSeq object
message("Create DEXSeq object ...")
dxd <- DEXSeqDataSetFromFeatureCounts(count_f, flattenedfile = flat_GTF, sampleData = target)

### DEU analysis
message("Estimate size factors and dispersion ...")
dxr1 <- DEXSeq(dxd, fullModel=design(dxd),
               reducedModel = ~ sample + exon,
               BPPARAM=MulticoreParam(workers=workers),
               fitExpToVar="condition",
               quiet=TRUE)
# dxd <- estimateSizeFactors(dxd)
# dxd <- estimateDispersions(dxd)
# message("Test for DEU ...")
# dxd  <- testForDEU(dxd)
# dxd <- estimateExonFoldChanges(dxd, fitExpToVar="condition")
# dxr1 <- DEXSeqResults(dxd)

message("Extract significant DEUs ...")
sigGenes.result <- as.data.frame(dxr1[which(dxr1$padj < 0.05), ])
sigGenes.result$transcripts <- sapply(sigGenes.result$transcripts, function(x) paste(x, collapse = ","))

### Remove ambiguous genes (genes have common exons)
message("Remove ambiguous genes (genes have common exons) ...")
amb.genes <- grep("\\+", sigGenes.result$groupID)
sigGenes.result.flt <- sigGenes.result[-amb.genes, ]
print(head(sigGenes.result.flt))

message("Save final DEU results ...")
write.table(sigGenes.result, file.path(DEXSeq_output, paste0("DEUs_detected.0.05.tsv")), sep="\t", quote=FALSE, row.names=FALSE)
write.table(sigGenes.result.flt, file.path(DEXSeq_output, paste0("DEUs_detected.0.05.flt.tsv")), sep="\t", quote=FALSE, row.names=FALSE)

