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

DEU_analysis <- function(count_f, flat_GTF, target) {
  
  message("Create DEXSeq object ...")
  dxd <- DEXSeqDataSetFromFeatureCounts(count_f, 
                                        flattenedfile = flat_GTF, 
                                        sampleData = target)
  
  message("Estimate size factors and dispersion ...")
  dxd <- estimateSizeFactors(dxd)
  dxd <- estimateDispersions(dxd, 
                             BPPARAM=MulticoreParam(workers=4))
  
  message("Test for DEU ...")
  dxd  <- testForDEU(dxd, 
                     BPPARAM=MulticoreParam(workers=4))
  dxd <- estimateExonFoldChanges(dxd, 
                                 fitExpToVar="condition", 
                                 BPPARAM=MulticoreParam(workers=4))
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

save_results <- function(res, res_flt, gene_res, fdr_cutoff) {
  
  message("Save DEXSeq result ...")
  write.table(res[!is.na(res$padj.gene), ], 
              "DEUs.allGenes.tsv", 
              sep="\t", quote=FALSE, row.names=FALSE)
  write.table(res_flt[!is.na(res_flt$padj.gene), ], 
              "DEUs.flt.allGenes.tsv", 
              sep="\t", quote=FALSE, row.names=FALSE)
  write.table(gene_res, 
              "DEUs.geneWise.allGenes.tsv", 
              sep="\t", quote=FALSE, row.names=FALSE)
  
  if (length(which(gene_res$padj.gene < fdr_cutoff)) > 0) {
    message("Save significant DEUs ...")
    sig_res <- res[which(res$padj.gene < fdr_cutoff), ]
    write.table(sig_res, 
                paste0("DEUs.sigGenes.", fdr_cutoff, ".tsv"), 
                sep="\t", quote=FALSE, row.names=FALSE)
    
    message("Save significant non-aggregated DEUs ...")    
    sig_res_flt <- res_flt[which(res_flt$padj.gene < fdr_cutoff), ]
    write.table(sig_res_flt, 
                paste0("DEUs.flt.sigGenes.", fdr_cutoff, ".tsv"), 
                sep="\t", quote=FALSE, row.names=FALSE)
  }
  
}

filter_gene <- function(REF, gene_res, fdr_cutoff, DEXSeq_o) {
  sigGene_res <- gene_res[gene_res$padj.gene < fdr_cutoff, ]

  message("Remove sig DEU genes without gene symbol ...")  
  gene_info_p <- paste0(REF, "/gene_info_GRCm39_M32_ensembl109.tsv")
  gene_info <- read.table(gene_info_p, header=T, na.strings="." )
  m <- match(sigGene_res$groupID, gene_info$GeneID)
  sigGene_res$symbol <- gene_info$Symbol[m]
  keep <- !is.na(sigGene_res$symbol)
  sigGene_res_flt <- sigGene_res[keep,]
  
  write.table(sigGene_res_flt, 
              paste0(DEXSeq_o, "DEUs.geneWise.sigGenes.", fdr_cutoff, ".withGeneSymbol.tsv"), 
              sep="\t", quote=FALSE, row.names=FALSE)   
}

# Getting parameters
args <- R.utils::commandArgs(asValues = TRUE)
print(args)
mode  <- as.character(args[['MODE']])
DIR <- as.character(args[['DIR']])
REF <- as.character(args[['REF']])
flat_GTF <- as.character(args[['flat_GTF']])
targetp <- as.character(args[['target']])
pair <- as.character(args[['pair']])
fdr_cutoff <- as.numeric(args[['fdr_cutoff']])
wd <- getwd()

message("Check parameters ...")
print(wd)
print(DIR)
print(REF)
print(flat_GTF)
print(targetp)
print(pair)
print(fdr_cutoff)

message("Reformat target.tsv ...")
target <- reformat_target(targetp, pair)

if (mode == "simulation") {
  
  seed <- as.integer(args[['seed']])
  noOfSim <- as.integer(args[['noOfSim']])
  workers <- as.integer(args[['workers']])
  simID <- paste0("S", 1:noOfSim)
  
  print(simID)
  print(seed)
  print(workers)
  
  ### 01. Register BPPARAM
  message('Register BPPARAM')
  BPPARAM <- MulticoreParam(workers = workers, progressbar = TRUE, RNGseed = seed)
  register(BPPARAM = BPPARAM)
  
  ### 02. Run DEXSeq
  message('Start running DEXSeq ...')
  bplapply(simID, FUN = function(i){
    
    message(paste0("Processing: ", i, " ..."))
    count_f <- paste0(DIR, "featureCounts_subread/", i, "/featureCounts_output.txt")
    print(count_f)
    DEXSeq_o <- paste0(DIR, "DEXSeq/", i)
    print(DEXSeq_o)
    dir.create(DEXSeq_o, showWarnings = FALSE, recursive = TRUE)
    
    message("Run DEXSeq ...")
    dxr <- DEU_analysis(count_f, flat_GTF, target)
    
    message("Post-process DEU results ...")
    final_res <- reformat_DEUres(dxr$dxr_exon, dxr$dxr_gene)
    
    setwd(DEXSeq_o)    
    message("Save DEXSeq result ...")
    res <- save_results(final_res$res, final_res$res_flt, final_res$gene_res, fdr_cutoff)
    setwd(wd)
    
  }, BPPARAM = BPPARAM)
  
} else if (mode == "case_study") {

  count_f <- paste0(DIR, "featureCounts_subread/", pair, "/featureCounts_output.txt") 
  print(count_f)
  DEXSeq_o <- paste0(DIR, "DEXSeq/", pair, "/")
  print(DEXSeq_o)
  dir.create(DEXSeq_o, showWarnings = FALSE, recursive = TRUE)
  
  message("Run DEXSeq ...")
  dxr <- DEU_analysis(count_f, flat_GTF, target)
  
  message("Save dxr object for visualization ...")
  g <- unique(dxr$dxr_exon[which(dxr$dxr_exon$padj < 0.05), "groupID"])
  dxr_rds <- dxr$dxr_exon[dxr$dxr_exon$groupID %in% g, ]
  saveRDS(dxr_rds, paste0(DEXSeq_o, "dxr_exon.rds"))
  
  message("Post-process DEU results ...")
  final_res <- reformat_DEUres(dxr$dxr_exon, dxr$dxr_gene)
  
  setwd(DEXSeq_o)
  message("Save DEXSeq result ...")
  res <- save_results(final_res$res, final_res$res_flt, final_res$gene_res, fdr_cutoff)
  setwd(wd)
  
  flt <- filter_gene(REF, final_res$gene_res, fdr_cutoff, DEXSeq_o)
  
} else {
  print("Invalid mode specified. Choose 'simulation' or 'case_study'!")
}


##################################################################################
# OK <- requireNamespace("devtools", quietly = TRUE)
# if (!OK) {
#   stop("devtools package required but is not installed (or can't be loaded)")
# }
# # source("/vast/projects/lab_chen/tam/tools/Subread_to_DEXSeq/load_SubreadOutput_modified.R")
# source("../../tools/Subread_to_DEXSeq/load_SubreadOutput_modified.R")
# library(dplyr)
# library(DEXSeq)
# library(BiocParallel)
# 
# # Getting parameters
# args <- R.utils::commandArgs(asValues = TRUE)
# print(args)
# mode  <- as.character(args[['MODE']])
# DIR <- as.character(args[['DIR']])
# flat_GTF <- as.character(args[['flat_GTF']])
# 
# if (mode == "simulation") {
#   target <- as.character(args[['target']])
#   seed <- as.integer(args[['seed']])
#   noOfSim <- as.integer(args[['noOfSim']])
#   libs <- as.integer(args[['libs']])
#   workers <- as.integer(args[['workers']])
#   simID <- paste0("S", 1:noOfSim)
#   
#   ### 0. Check parameter
#   message("Check valid parameters")
#   print(DIR)
#   print(flat_GTF)
#   print(simID)
#   print(target)
#   print(seed)
#   print(workers)
#   print(libs)
#   
#   ### 01. Register BPPARAM
#   message('Register BPPARAM')
#   BPPARAM <- MulticoreParam(workers = workers, progressbar = TRUE, RNGseed = seed)
#   register(BPPARAM = BPPARAM)
#   
#   ### 02. Run DEXSeq
#   message('Start running DEXSeq ...')
#   bplapply(simID, FUN = function(i){
#     
#     message(paste0("Processing: ", i, " ..."))
#     count_f <- paste0(DIR, "/STAR_output/featureCounts_subread/", i, "/featureCounts_output.txt")
#     print(count_f)
#     DEXSeq_output <- paste0(DIR, "/STAR_output/DEXSeq/", i)
#     print(DEXSeq_output)
#     dir.create(DEXSeq_output, showWarnings = FALSE, recursive = TRUE)
# 
#     message('Reformating target.tsv ...')
#     rownames(target) <- target$sampleID
#     target <- target[-1]
#     colnames(target)[1] <- "condition"
#     
#     ### Create DEXSeq object
#     message("Create DEXSeq object ...")
#     dxd <- DEXSeqDataSetFromFeatureCounts(count_f, flattenedfile = flat_GTF, sampleData = target)
#     
#     ### DEU analysis
#     message("Estimate size factors and dispersion ...")
#     dxd <- estimateSizeFactors(dxd)
#     dxd <- estimateDispersions(dxd, BPPARAM=MulticoreParam(workers=4))
#     message("Test for DEU ...")
#     dxd  <- testForDEU(dxd, BPPARAM=MulticoreParam(workers=4))
#     dxd <- estimateExonFoldChanges(dxd, fitExpToVar="condition", BPPARAM=MulticoreParam(workers=4))
#     dxr1 <- DEXSeqResults(dxd)
#     
#     message("Save DEXSeq result ...")
#     dxr1.result  <- as.data.frame(dxr1)
#     dxr1.result$transcripts <- sapply(dxr1.result$transcripts, function(x) paste(x, collapse = ","))
#     write.table(dxr1.result, file.path(DEXSeq_output, paste0("DEUs_detected.allGenes.tsv")), sep="\t", quote=FALSE, row.names=FALSE)
#     
#     if (length(which(dxr1$padj < 0.05)) > 0) {
#       message("Extract significant DEUs ...")
#       sigGenes.result <- as.data.frame(dxr1[which(dxr1$padj < 0.05), ])
#       sigGenes.result$transcripts <- sapply(sigGenes.result$transcripts, function(x) paste(x, collapse = ","))
#       
#       ### Remove ambiguous genes (genes have common exons)
#       message("Remove ambiguous genes (genes have common exons) ...")
#       amb.genes <- grep("\\+", sigGenes.result$groupID)
#       sigGenes.result.flt <- sigGenes.result[-amb.genes, ]
#       
#       message("Save final DEU results ...")
#       write.table(sigGenes.result, file.path(DEXSeq_output, paste0("DEUs_detected.0.05.tsv")), sep="\t", quote=FALSE, row.names=FALSE)
#       write.table(sigGenes.result.flt, file.path(DEXSeq_output, paste0("DEUs_detected.0.05.flt.tsv")), sep="\t", quote=FALSE, row.names=FALSE)
#     }
#     
#   }, BPPARAM = BPPARAM)
#    
# } else if(mode == "case_study") {
#   
#   target <- as.character(args[['target_gp']])
#   count_f <- paste0(DIR, "/STAR_output/featureCounts_subread/featureCounts_output_", gp, ".txt")
#   DEXSeq_output <- paste0(DIR, "/STAR_output/DEXSeq/")
#   dir.create(DEXSeq_output, recursive = T, showWarnings = F)
# 
#   message('Reformating target.tsv ...')
#   rownames(target) <- target$sampleID
#   target <- target[-1]
#   colnames(target)[1] <- "condition"
#   
#   message('Create DEXSeq object ...')
#   dxd <- DEXSeqDataSetFromFeatureCounts(count_f, flattenedfile = flat_GTF, sampleData = target)
#   
#   message('DEU analysis ...')
#   dxd <- estimateSizeFactors(dxd)
#   dxd <- estimateDispersions(dxd, BPPARAM=MulticoreParam(workers=4))
#   dxd  <- testForDEU(dxd, BPPARAM=MulticoreParam(workers=4))
#   dxd <- estimateExonFoldChanges(dxd, fitExpToVar="condition", BPPARAM=MulticoreParam(workers=4))
#   dxr1 <- DEXSeqResults(dxd)
#   
#   message('Save DEXSeq results ...')
#   dxr1.result  <- as.data.frame(dxr1)
#   dxr1.result$transcripts <- sapply(dxr1.result$transcripts, function(x) paste(x, collapse = ","))
#   write.table(dxr1.result, 
#               file.path(paste0(DEXSeq_output, gp), paste0("DEUs_detected_", gp, ".allGenes.tsv")), 
#               sep="\t", quote=FALSE, row.names=FALSE)
#   write.csv(dxr1.result, 
#             file.path(paste0(DEXSeq_output, gp), paste0("DEUs_detected_", gp, ".allGenes.csv")),
#             row.names=FALSE)
#   
#   if (length(which(dxr1$padj < 0.05)) > 0) {
#     
#     message("Extract significant DEUs ...")
#     sigGenes.result <- as.data.frame(dxr1[which(dxr1$padj < 0.05), ])
#     sigGenes.result$transcripts <- sapply(sigGenes.result$transcripts, function(x) paste(x, collapse = ","))
#     
#     message("Remove ambiguous genes (genes have common exons) ...")
#     amb.genes <- grep("\\+", sigGenes.result$groupID)
#     sigGenes.result.flt <- sigGenes.result[-amb.genes, ]
#     
#     message('Remove genes without gene symbol ...')
#     m <- match(sigGenes.result.flt$groupID, gene_info$GeneID)
#     sigGenes.result.flt$symbol <- gene_info$Symbol[m]
#     keep <- !is.na(sigGenes.result.flt$symbol)
#     table(keep)
#     sigGenes.result.flt <- sigGenes.result.flt[keep, ]
#     
#     message("Save final significant DEU results ...")
#     write.table(sigGenes.result, 
#                 file.path(paste0(DEXSeq_output, gp), paste0("DEUs_detected_", gp, ".0.05.tsv")), 
#                 sep="\t", quote=FALSE, row.names=FALSE)
#     write.csv(sigGenes.result, 
#               file.path(paste0(DEXSeq_output, gp), paste0("DEUs_detected_", gp, ".0.05.csv")), 
#               row.names=FALSE)
#     write.table(sigGenes.result.flt, 
#                 file.path(paste0(DEXSeq_output, gp), paste0("DEUs_detected_", gp, ".0.05.flt.tsv")), 
#                 sep="\t", quote=FALSE, row.names=FALSE)
#     write.csv(sigGenes.result.flt, 
#               file.path(paste0(DEXSeq_output, gp), paste0("DEUs_detected_", gp, ".0.05.flt.csv")), 
#               row.names=FALSE)
#     
#     message("Save DEXSeq object for visualization ...")
#     DEUs_genes <- unique(dxr1[which(dxr1$padj < 0.05), "groupID"])
#     print(length(DEUs_genes))
#     dxr1_rds <- dxr1[dxr1$groupID %in% DEUs_genes, ]
#     saveRDS(dxr1_rds, file = paste0(DEXSeq_output, gp, "/dxr1_", gp, ".rds"))
#   }
#   
# } else {
#   
#   print("Invalid mode specified. Choose 'simulation' or 'case_study'!")
# }
