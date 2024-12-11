OK <- requireNamespace("devtools", quietly = TRUE)
if (!OK) {
  stop("devtools package required but is not installed (or can't be loaded)")
}
source("../../tools/Subread_to_DEXSeq/load_SubreadOutput_modified.R")
library(dplyr)
library(DEXSeq)
library(BiocParallel)

reformat_target <- function(targetp, pair) {
  target <- read.table(targetp, header=T)
  rownames(target) <- target$sampleID
  target <- target[-1]
  colnames(target)[1] <- "condition" 
  p <- strsplit(pair, "-")[[1]] 
  target$condition <- factor(target$condition, levels=c(p[2], p[1]))
  print(target$condition)
  return(target)
}

DEU_analysis <- function(count_f, flat_GTF, target, workers) {
  
  message("Create DEXSeq object ...")
  dxd <- DEXSeqDataSetFromFeatureCounts(count_f, 
                                        flattenedfile = flat_GTF, 
                                        sampleData = target)
  
  message("Estimate size factors and dispersion ...")
  dxd <- estimateSizeFactors(dxd)
  dxd <- estimateDispersions(dxd, 
                             BPPARAM=MulticoreParam(workers=workers))
  
  message("Test for DEU ...")
  dxd  <- testForDEU(dxd, 
                     BPPARAM=MulticoreParam(workers=workers))
  dxd <- estimateExonFoldChanges(dxd, 
                                 fitExpToVar="condition", 
                                 BPPARAM=MulticoreParam(workers=workers))
  dxr_exon <- DEXSeqResults(dxd)
  
  message("Calculate FDR at gene-level ...")
  dxr_gene <- perGeneQValue(dxr_exon)
  
  return(list("dxr_exon"=dxr_exon, "dxr_gene"=dxr_gene))
}

reformat_DEUres <- function(dxr_exon, dxr_gene) {
  
  message("Reformat results ...")
  exon_res  <- as.data.frame(dxr_exon)
  exon_res$transcripts <- sapply(exon_res$transcripts, 
                                 function(x) paste(x, collapse = ","))
  
  gene_res <- data.frame(groupID=names(dxr_gene), padj.gene=dxr_gene)
  rownames(gene_res) <- NULL
  
  message("Merge padj at gene level to final DEU result ...")
  print(head(exon_res))
  print(head(gene_res))
  res <- merge(exon_res, gene_res, by="groupID", all=TRUE)
  
  message("Filter aggregated genes ...")
  agg_genes <- grep("\\+", res$groupID)
  res_flt <- res[-agg_genes, ]
  
  return(list("res"=res, "res_flt"=res_flt, "gene_res"=gene_res))
}

# Getting parameters
args <- R.utils::commandArgs(asValues = TRUE)
print(args)
DIR <- as.character(args[['DIR']])
REF <- as.character(args[['REF']])
flat_GTF <- as.character(args[['flat_GTF']])
targetp <- as.character(args[['target']])
pair <- as.character(args[['pair']])
workers <- as.integer(args[['workers']])
i <- "S1"

message("Check parameters ...")
print(DIR)
print(REF)
print(flat_GTF)
print(targetp)
print(pair)

message("Reformat target.tsv ...")
target <- reformat_target(targetp, pair)
    
message(paste0("Processing: ", i, " ..."))
count_f <- paste0(DIR, "featureCounts_subread/", i, "/featureCounts_output.txt")
print(count_f)
DEXSeq_o <- paste0(DIR, "DEXSeq/", i)
print(DEXSeq_o)
dir.create(DEXSeq_o, showWarnings = FALSE, recursive = TRUE)

message("Run DEXSeq ...")
dxr <- DEU_analysis(count_f, flat_GTF, target, workers)

message("Post-process DEU results ...")
final_res <- reformat_DEUres(dxr$dxr_exon, dxr$dxr_gene)
