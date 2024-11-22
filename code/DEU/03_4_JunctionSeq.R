OK <- requireNamespace("devtools", quietly = TRUE)
if (!OK) {
  stop("devtools package required but is not installed (or can't be loaded)")
}
library(JunctionSeq)

read_target <- function(targetp, pair) {
  p <- strsplit(pair, "-")[[1]]
  target <- read.table(targetp, header=TRUE, stringsAsFactors=FALSE)
  design <- data.frame(condition = factor(target$group.ID, levels=c(p[2], p[1])))
  
  return(list("target"=target, "design"=design))
}

DEU_analysis <- function(count_files, target, design, GFF, JunctionSeq_o, fdr_cutoff, mode, ncores) {
  
  message("readJunctionSeqCounts ...")
  jscs <- readJunctionSeqCounts(countfiles = count_files,
                                samplenames = target$sample.ID,
                                design = design,
                                flat.gff.file = GFF)
  
  message("Estimate size factors and dispersion ...")
  jscs <- estimateJunctionSeqSizeFactors(jscs)
  jscs <- estimateJunctionSeqDispersions(jscs, nCores = ncores)
  jscs <- fitJunctionSeqDispersionFunction(jscs)
  
  message("Test for DEU ...")
  jscs <- testForDiffUsage(jscs, 
                           nCores = ncores, 
                           optimizeFilteringForAlpha = 0.05)
  jscs <- estimateEffectSizes( jscs, nCores = ncores)
  
  if (mode == "case_study") saveRDS(jscs, paste0(JunctionSeq_o, "jscs.rds"))
  
  message("Write DEU results ...")
  writeCompleteResults(jscs,
                       outfile.prefix = JunctionSeq_o,
                       FDR.threshold = fdr_cutoff,
                       save.allGenes = TRUE,
                       save.sigGenes = TRUE,
                       save.bedTracks = FALSE)
}

filter_res <- function(REF, JunctionSeq_o, fdr_cutoff) {
  
  message("Remove genes without gene symbol ...")
  gene_info_p <- paste0(REF, "/gene_info_GRCm39_M32_ensembl109.tsv")
  gene_info <- read.table(gene_info_p, header=T, na.strings="." )
  DEUg <- read.table(paste0(JunctionSeq_o, "sigGenes.genewiseResults.txt.gz"), header=TRUE)
  m <- match(DEUg$geneID, gene_info$GeneID)
  DEUg$symbol <- gene_info$Symbol[m]
  keep <- !is.na(DEUg$symbol)
  DEUg.flt <- DEUg[keep,]
  write.table(DEUg.flt, 
              paste0(JunctionSeq_o, "sigGenes_", fdr_cutoff, ".genewiseResults.flt.withGeneSymbol.txt"), 
              sep="\t", quote=FALSE, row.names=FALSE)
  
}

# Getting parameters
args <- R.utils::commandArgs(asValues = TRUE)
print(args)

DIR <- as.character(args[['DIR']])
REF <- as.character(args[['REF']])
targetp <- as.character(args[['target']])
ncores <- as.integer(args[['ncores']])
mode <- as.character(args[['MODE']])
fdr_cutoff <- as.numeric(args[['fdr_cutoff']])
pair <- as.character(args[['pair']])
target <- read_target(targetp, pair)

message("Check parameters ...")
print(DIR)
print(REF)
print(targetp)
print(target$target)
print(target$design)
print(pair)
print(ncores)
print(fdr_cutoff)

if (mode == "simulation") {
  noOfSim <- as.integer(args[['noOfSim']])
  simID <- paste0("S", 1:noOfSim)
  print(simID)
  
  for (i in simID) {
    message(paste0('Processing: ', i, " ..."))
    JunctionSeq_o <- paste0(DIR, "JunctionSeq/", i, "/")
    dir.create(JunctionSeq_o, showWarnings = FALSE, recursive = TRUE)
    QoRTs_o <- paste0(DIR, "QoRTs/", i, "/")
    
    count_files <- paste0(QoRTs_o, target$target$sample.ID, 
                          "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz")
    GFF <- paste0(QoRTs_o, "/withNovel.forJunctionSeq.gff.gz")
    
    print(count_files)
    print(GFF)
    
    message("Run DEU analysis ...")
    res <- DEU_analysis(count_files, target$target, target$design, GFF, JunctionSeq_o, fdr_cutoff, mode, ncores)
  }
  
} else if (mode == "case_study") {
  
  JunctionSeq_o <- paste0(DIR, "JunctionSeq/", pair, "/")
  dir.create(JunctionSeq_o, showWarnings = FALSE, recursive = TRUE)
  QoRTs_o <- paste0(DIR, "QoRTs/")
  
  count_files <- paste0(QoRTs_o, target$target$sample.ID, 
                        "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz")
  GFF <- paste0(QoRTs_o, "/withNovel.forJunctionSeq.gff.gz")
  
  print(count_files)
  print(GFF)
  
  message("Run DEU analysis ...")
  res <- DEU_analysis(count_files, target$target, target$design, GFF, JunctionSeq_o, fdr_cutoff, mode, ncores)
  
  message("Filter DEU genes without gene symbol ...")
  res_flt <- filter_res(REF, JunctionSeq_o, fdr_cutoff)
  
} else {
  
  print("Invalid mode specified. Choose 'simulation' or 'case_study'!")
  
}

###############################################################################
# OK <- requireNamespace("devtools", quietly = TRUE)
# if (!OK) {
#   stop("devtools package required but is not installed (or can't be loaded)")
# }
# library(JunctionSeq)
# 
# # Getting parameters
# args <- R.utils::commandArgs(asValues = TRUE)
# print(args)
# 
# # DIR <- "/vast/projects/Spatial/tam/Differential_splicing/simulation_3vs3_mouse_genome/RNA-seq/DEU_mix/unbalanced_3_75_3_sample/"
# # targetp <- "/vast/projects/lab_chen/tam/Differential_splicing/simulation_3vs3_mouse_genome/RNA-seq/config/decoder.bySample.3vs3.txt"
# # ncores <- 8
# 
# DIR <- as.character(args[['DIR']])
# targetp <- as.character(args[['target']])
# ncores <- as.integer(args[['ncores']])
# mode <- as.integer(args[['MODE']])
# 
# target <- read.table(targetp, header=TRUE, stringsAsFactors=FALSE)
# design <- data.frame(condition = factor(target$group.ID))
# 
# if (mode == "simulation") {
#   noOfSim <- as.integer(args[['noOfSim']])
#   simID <- paste0("S", 1:noOfSim)
#   
#   ### 0. Check parameter
#   message("Check valid parameters")
#   print(DIR)
#   print(targetp)
#   print(target)
#   print(simID)
#   print(design)
#   print(ncores)
#   
#   for (i in simID) {
#     message(paste0('Processing: ', i, " ..."))
#     count_files <- paste0(DIR, "STAR_output/QoRTs_fit_RL_75/", i, "/", target$sample.ID, "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz")
#     GFF <- paste0(DIR, "STAR_output/QoRTs_fit_RL_75/", i, "/withNovel.forJunctionSeq.gff.gz")
#     OUTPUT <- paste0(DIR, "STAR_output/JunctionSeq_fit_RL_75/", i, "/")
#     dir.create(OUTPUT, showWarnings = FALSE, recursive = TRUE)
#     
#     print(count_files)
#     print(GFF)
#     print(OUTPUT)
#     
#     ### DEU analysis
#     message("readJunctionSeqCounts ...")
#     jscs <- readJunctionSeqCounts(countfiles = count_files,
#                                   samplenames = target$sample.ID,
#                                   design = design,
#                                   flat.gff.file = GFF)
#     message("Estimate size factors and dispersion ...")
#     jscs <- estimateJunctionSeqSizeFactors(jscs)
#     jscs <- estimateJunctionSeqDispersions(jscs, nCores = ncores)
#     jscs <- fitJunctionSeqDispersionFunction(jscs)
#     message("Test for DEU ...")
#     jscs <- testForDiffUsage(jscs, nCores = ncores, optimizeFilteringForAlpha = 0.05)
#     jscs <- estimateEffectSizes( jscs, nCores = ncores)
#     message("Write DEU results ...")
#     writeCompleteResults(jscs,
#                          outfile.prefix = OUTPUT,
#                          FDR.threshold = 0.05,
#                          save.allGenes = TRUE,
#                          save.sigGenes = TRUE,
#                          save.bedTracks = FALSE)
#   }
#   
# } else if (mode == "case_study") {
#   
#   gp <- as.character(args[['gp']])
#   REF <- as.character(args[['REF']])
#   
#   currentDir <- getwd()
#   
#   gene_info_p <- paste0(REF, "/gene_info_GRCm39_M32_ensembl109.tsv")
#   gene_info <- read.table(gene_info_p, header=T, na.strings="." )
#     
#   GFF <- paste0(DIR, "STAR_output/QoRTs/withNovel.forJunctionSeq.gff.gz")
#   OUTPUT <- paste0(DIR, "STAR_output/JunctionSeq/", gp, "/")
#   dir.create(OUTPUT, recursive=TRUE)
#   
#   count_files <- c()
#   for (sample in target$sample.ID) {
#     count_files <- c(count_files, paste0(DIR, "STAR_output/QoRTs/", 
#                                          sample, "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz"))
#   }
#   
#   # Check parameters
#   print(gp)
#   print(head(gene_info))
#   print(count_files)
#   
#   message('DEU analysis ...')
#   jscs <- readJunctionSeqCounts(countfiles = count_files,
#                                 samplenames = target$sample.ID,
#                                 design = design,
#                                 flat.gff.file = GFF)
#   jscs <- estimateJunctionSeqSizeFactors(jscs)
#   jscs <- estimateJunctionSeqDispersions(jscs, nCores = 4)
#   jscs <- fitJunctionSeqDispersionFunction(jscs)
#   jscs <- testForDiffUsage(jscs, nCores = 4)
#   jscs <- estimateEffectSizes( jscs, nCores = 4)
#   writeCompleteResults(jscs,
#                        outfile.prefix = OUTPUT,
#                        FDR.threshold = 0.05,
#                        save.allGenes = TRUE, 
#                        save.sigGenes = TRUE,
#                        save.bedTracks = FALSE)
#   saveRDS(jscs, file=paste0(OUTPUT, "jscs.rds"))
#   
#   message('Remove genes without gene symbol ...')
#   DEU_detected <- read.table(file.path(OUTPUT, "sigGenes.results.txt.gz"), header = TRUE)
#   m <- match(DEU_detected$geneID, gene_info$GeneID)
#   DEU_detected$symbol <- gene_info$Symbol[m]
#   keep <- !is.na(DEU_detected$symbol)
#   table(keep)
#   DEU_detected <- DEU_detected[keep,]
#   length(unique(DEU_detected$geneID))
#   write.table(DEU_detected, file.path(OUTPUT, "sigGenes.withGeneName.results.txt"), sep="\t", quote=FALSE, row.names=FALSE)
# 
# } else {
#   
#   print("Invalid mode specified. Choose 'simulation' or 'case_study'!")
#   
# }
