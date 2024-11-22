OK <- requireNamespace("devtools", quietly = TRUE)
if (!OK) {
  stop("devtools package required but is not installed (or can't be loaded)")
}
library(Rsubread)
library(edgeR)
library(BiocParallel)

quantify_IE_count <- function(BAM_files, isPairedEnd, flat_exon, featureCounts_o, wd) {
  
  message("Running featureCounts with nonSplitOnly=TRUE and splitOnly=FALSE ...")
  internal_exon_count <- Rsubread::featureCounts(BAM_files,
                                                 annot.ext=flat_exon,
                                                 tmpDir=featureCounts_o,
                                                 nonSplitOnly=TRUE,
                                                 splitOnly=FALSE,
                                                 useMetaFeatures=FALSE,
                                                 allowMultiOverlap=TRUE,
                                                 isPairedEnd=isPairedEnd,
                                                 nthreads=8,
                                                 minMQS=255)
  
  message("Saving internal exon read counts and annotations ...")
  setwd(featureCounts_o)
  write.table(internal_exon_count$counts, 
              "internal_exon_count.tsv", 
              sep="\t", quote=FALSE, row.names=FALSE)
  write.table(internal_exon_count$annotation, 
              "internal_exon_count_annotation.tsv", 
              sep="\t", quote=FALSE, row.names=FALSE)
  write.table(internal_exon_count$stat, 
              "internal_exon_count_stat.tsv", 
              sep="\t", quote=FALSE, row.names=FALSE)
  setwd(wd)
}

quantify_E_J_count <- function(BAM_files, isPairedEnd, flat_exon, FASTA, featureCounts_o, wd) {
  
  message("Running featureCounts with juncCounts=TRUE ...")
  junction_count <- Rsubread::featureCounts(BAM_files,
                                            annot.ext=flat_exon,
                                            tmpDir=featureCounts_o,
                                            genome=FASTA,
                                            juncCounts=TRUE,
                                            useMetaFeatures=FALSE,
                                            allowMultiOverlap=TRUE,
                                            isPairedEnd=isPairedEnd,
                                            nthreads=8,
                                            minMQS=255)
  
  message("Saving exon-level read counts and annotations with duplicated junction read counts ...")
  setwd(featureCounts_o)
  write.table(junction_count$counts, 
              "exon_count.tsv", 
              sep="\t", quote=FALSE, row.names=FALSE)
  write.table(junction_count$annotation, 
              "exon_count_annotation.tsv", 
              sep="\t", quote=FALSE, row.names=FALSE)
  write.table(junction_count$stat, 
              "exon_count_stat.tsv", 
              sep="\t", quote=FALSE, row.names=FALSE)
  
  message("Saving junction read counts and annotations ...")  
  write.table(junction_count$counts_junction[!is.na(junction_count$counts_junction$PrimaryGene),], 
              "junction_count.tsv", 
              sep="\t", quote=FALSE, row.names=FALSE)
  setwd(wd)
  
}

filter_gene <- function(y, REF) {
  
  message("Removing genes without gene symbol for real data ...")
  gene_info_p <- paste0(REF, "gene_info_GRCm39_M32_ensembl109.tsv")
  gene_info <- read.table(gene_info_p, header=T, na.strings="." )
  m <- match(y$genes$GeneID, gene_info$GeneID)
  y$genes$Symbol <- gene_info$Symbol[m]
  keep <- !is.na(y$genes$Symbol)
  table(keep)
  y <- y[keep, , keep.lib.sizes=FALSE]  
  return(y)
}

construct_model_matrix <- function(targetp) {
  
  message("Reading target ...")
  target <- read.table(targetp, sep="\t", header=TRUE, stringsAsFactors=FALSE)
  group <- factor(target$group)
  samples <- target$sampleID
  
  message("Constructing model matrix ...")
  design <- model.matrix(~ 0 + group)
  colnames(design) <- gsub("group", "", colnames(design))
  
  return(list('group'=group, 'samples'=samples, 'design'=design))
}

process_E_count_withDupJ <- function(featureCounts_o) {
  
  E_count <- read.table(paste0(featureCounts_o, "exon_count.tsv"), header=TRUE)
  E_count <- data.frame(E_count[,1:ncol(E_count)], row.names=NULL, check.names = FALSE)
  E_annot <- read.table(paste0(featureCounts_o, "exon_count_annotation.tsv"), header=TRUE)
  
  return(list('E_count'=E_count,
              'E_annot'=E_annot))
}

reassign_SJ <- function(J_count, SJ_database) {
  
  message("Only consider unique junctions ...")
  uniq_SJ <- SJ_database[SJ_database$freq==1,]
  m1 <- match(J_count$juncID, uniq_SJ$juncID)
  J_count$PrimaryGene <- ifelse(!is.na(m1), uniq_SJ$geneID[m1], J_count$PrimaryGene)
  
  message("Whether junctions are annotated or not ...")
  m2 <- match(J_count$juncID, SJ_database$juncID)
  J_count$annotated <- ifelse(!is.na(m2), 1, 0)
  
  return(J_count) 
}

process_IE_J_count <- function(featureCounts_o, SJ_database) {
  
  message("Analyzing internal exon read counts ...")
  IE_count <- read.table(paste0(featureCounts_o, "internal_exon_count.tsv"), header=TRUE)
  IE_count <- data.frame(IE_count[,1:ncol(IE_count)], row.names=NULL, check.names = FALSE)
  IE_annot <- read.table(paste0(featureCounts_o, "internal_exon_count_annotation.tsv"), header=TRUE)
  IE_annot <- cbind(IE_annot, Region="Exon", annotated=1)
  
  message("Analyzing junction read counts ...")
  J_count <- read.table(paste0(featureCounts_o, "junction_count.tsv"), header=TRUE)
  J_count$juncID <- paste(J_count$Site1_chr, 
                          J_count$Site1_location, 
                          J_count$Site2_location, 
                          sep="_")
  
  message("Re-assigning junctions to correct genes using junction database ...")
  J_count <- reassign_SJ(J_count, SJ_database)
  
  message("Processing final junction count table ...")
  J_annot <- data.frame(
    GeneID=J_count$PrimaryGene,
    Chr=J_count$Site1_chr,
    Start=J_count$Site1_location,
    End=J_count$Site2_location
  )
  m <- match(J_annot$GeneID, IE_annot$GeneID)
  Strand <- IE_annot$Strand[m]
  J_annot <- cbind(J_annot, Strand=Strand, Length=1, Region="Junction", annotated=J_count$annotated)
  
  J_count <- data.frame(J_count[, 9:(ncol(J_count)-2)], row.names = NULL, check.names = FALSE)
  
  message("Combining internal exon and junction read count table ...")  
  IE_J_count <- rbind(IE_count, J_count)
  IE_J_annot <- rbind(IE_annot, J_annot)
  
  return(list('IE_J_count'=IE_J_count,
              'IE_J_annot'=IE_J_annot))
}

DEU_analysis <- function(r_count, annot, group, design, mode) {
  
  message("Constructing DGElist object ...")  
  y <- DGEList(counts=r_count, genes=annot, group=group)
  colnames(y) <- gsub("[.].*$", "", colnames(y))
  
  if (mode == "case_study") y <- filter_gene(y, REF)
  
  message("Filtering exons with low mapping reads ...")
  keep <- filterByExpr(y, group=group)
  y <- y[keep, , keep.lib.sizes=FALSE]
  
  message("Normalizing lib sizes ...")
  y <- normLibSizes(y)
  
  message("Estimating dispersion ...")
  y <- estimateDisp(y, design, robust=TRUE)
  
  message("Fitting GLM-QL model for design matrix ...")
  fit <- glmQLFit(y, design, robust=TRUE)
  
  return(fit)
}

DEU_res <- function(fit, design, p, m, fdr_cutoff, mode) {
  
  message("Running diffSpliceDGE ..")
  cmd <- paste("makeContrasts(", p, ",levels=design)", sep="")
  contr <- eval(parse(text=cmd))
  print(contr)
  sp <- diffSpliceDGE(fit, contrast=contr, geneid="GeneID", exonid="Start")
  DEU_simes <- topSpliceDGE(sp, test="Simes", number=Inf)
  DEU_F <- topSpliceDGE(sp, test="gene", number=Inf)
  DEU_exon <- topSpliceDGE(sp, test="exon", number=Inf)
  
  message("Extracting significant DEU genes - FDR cutoff 0.05 ...")  
  sigDEU_simes <- DEU_simes[DEU_simes$FDR<=fdr_cutoff,]
  sigDEU_F <- DEU_F[DEU_F$FDR<=fdr_cutoff,]
  sigDEU_exon <- DEU_exon[DEU_exon$FDR<=fdr_cutoff,]
  
  message("Saving all DEU results ...")
  write.table(DEU_simes, 
              paste0(m, "_simes_test.allGenes.tsv"), 
              sep="\t", quote=FALSE, row.names=FALSE)
  write.table(DEU_F, 
              paste0(m, "_F_test.allGenes.tsv"), 
              sep="\t", quote=FALSE, row.names=FALSE)
  write.table(DEU_exon, 
              paste0(m, "_exon_test.allGenes.tsv"), 
              sep="\t", quote=FALSE, row.names=FALSE)
  
  message("Saving significant DEU results ...")
  write.table(sigDEU_simes, 
              paste0(m, "_simes_test.sigGenes_", fdr_cutoff, ".tsv"), 
              sep="\t", 
              quote=FALSE, row.names=FALSE)  
  write.table(sigDEU_F, 
              paste0(m, "_F_test.sigGenes_", fdr_cutoff, ".tsv"), 
              sep="\t", quote=FALSE, row.names=FALSE)  
  write.table(sigDEU_exon, 
              paste0(m, "_exon_test.sigGenes_", fdr_cutoff, ".tsv"), 
              sep="\t", quote=FALSE, row.names=FALSE) 
  
  if (mode == "case_study") {
    message("Saving DEU results in csv ...")
    write.csv(DEU_simes, 
              paste0(m, "_simes_test.allGenes.csv"), 
              row.names=FALSE)
    write.csv(DEU_F, 
              paste0(m, "_F_test.allGenes.csv"), 
              row.names=FALSE)
    write.csv(DEU_exon, 
              paste0(m, "_exon_test.allGenes.csv"), 
              row.names=FALSE)
    
    write.csv(sigDEU_simes, 
              paste0(m, "_simes_test.sigGenes_", fdr_cutoff, ".csv"), 
              row.names=FALSE)
    write.csv(sigDEU_F, 
              paste0(m, "_F_test.sigGenes_", fdr_cutoff, ".csv"), 
              row.names=FALSE)
    write.csv(sigDEU_exon, 
              paste0(m, "_exon_test.sigGenes_", fdr_cutoff, ".csv"), 
              row.names=FALSE)
  }
  
}

# Getting parameters
args <- R.utils::commandArgs(asValues = TRUE)
print(args)

wd <- getwd()
DIR <- as.character(args[['DIR']])
REF <- as.character(args[['REF']])
targetp <- as.character(args[['target']])
pair <- as.character(args[['pair']])
mode <- as.character(args[['mode']])
isPairedEnd <- as.logical(args[['isPairedEnd']])
fdr_cutoff <- as.numeric(args[['fdr_cutoff']])

# DIR <- "/vast/projects/Spatial/tam/Differential_splicing/github/data/simulation/unbalanced_3_75_3/"
# REF <- "/vast/projects/Spatial/tam/Differential_splicing/github/annotation/"
# targetp <- "/vast/projects/Spatial/tam/Differential_splicing/github/data/simulation/target/target.3vs3.tsv"
# p <- "Group_2-Group_1"
# mode <- "simulation"
# isPairedEnd <- TRUE
# fdr_cutoff <- 0.05
# noOfSim <- 1
# seed <- 2024
# workers <- 4

FASTA <- paste0(REF, "GRCm39.primary_assembly.genome.fa.gz")
flat_exon <- paste0(REF, "gencode.vM32.annotation.flattened.exon.saf")
SJ <- paste0(REF, "gencode.vM32.primary_assembly.annotation.SJdatabase.tsv")
SJ_database <- read.table(SJ, header=TRUE)

message("Check parameters ...")
print(wd)
print(DIR)
print(REF)
print(targetp)
print(FASTA)
print(flat_exon)
print(SJ)
print(pair)
print(mode)
print(isPairedEnd)
print(fdr_cutoff)

if (mode == "simulation") {
  
  seed <- as.integer(args[['seed']])
  workers <- as.integer(args[['workers']])
  noOfSim <- as.integer(args[['noOfSim']])
  simID <- paste0("S", 1:noOfSim)
  
  print(simID)
  print(seed)
  print(workers)
  
  ### 01. Register BPPARAM
  message('Register BPPARAM')
  BPPARAM <- MulticoreParam(workers = workers, progressbar = TRUE, RNGseed = seed)
  register(BPPARAM = BPPARAM)
  
  ### 02. Run edgeR
  message('Start running edgeR ...')
  bplapply(simID, FUN = function(i){
    
    message(paste0('Processing: ', i, " ..."))
    
    message("Creating output folders to store read counts and DEU results ...")
    featureCounts_o <- paste0(DIR, "featureCounts/", i, "/")
    dir.create(featureCounts_o, showWarnings = FALSE, recursive = TRUE)
    DEU_analysis_o <- paste0(DIR, "edgeR_diffSpliceDGE/", i, "/")
    dir.create(DEU_analysis_o, showWarnings = FALSE, recursive = TRUE)
    
    message("Constructing design matrix for contrast ...")
    mat <- construct_model_matrix(targetp)
    
    BAM_files <- paste0(DIR, 
                        "aligned_pass2/", 
                        i, "/", 
                        mat$samples, "/", 
                        mat$samples, 
                        ".Aligned.sortedByCoord.out.bam")
    
    print(BAM_files)
    
    message("Checking if running featureCounts yet ...")
    f <- paste0(featureCounts_o, c("internal_exon_count.tsv", 
                                   "internal_exon_count_annotation.tsv", 
                                   "internal_exon_count_stat.tsv",
                                   "exon_count.tsv", 
                                   "exon_count_annotation.tsv", 
                                   "exon_count_stat.tsv",
                                   "junction_count.tsv"))
    print(any(!file.exists(f)))
    if (any(!file.exists(f))) {
      
      message("Quantifing internal exon read count ...")
      IE_count_table <- quantify_IE_count(BAM_files, isPairedEnd, flat_exon, featureCounts_o, wd)
      
      message("Quantifying exon-level read with duplicated junction count and junction read count ...")
      E_J_count_table <- quantify_E_J_count(BAM_files, isPairedEnd, flat_exon, FASTA, featureCounts_o, wd)
      
    }
    
    message("Constructing a matrix for exon read counts with duplicated junction count ...")
    E <- process_E_count_withDupJ(featureCounts_o)
    
    message("Starting DEU analysis for DEU-edgeR ...")
    DEU_fit <- DEU_analysis(E$E_count,
                            E$E_annot,
                            mat$group, mat$design, mode)

    message("Combining internal exon and junction read count matrix ...")  
    IE_J <- process_IE_J_count(featureCounts_o, SJ_database)
    
    message("Starting DEJU analysis for DEJU-edgeR ...")
    DEJU_fit <- DEU_analysis(IE_J$IE_J_count,
                             IE_J$IE_J_annot,
                             mat$group, mat$design, mode)

    setwd(DEU_analysis_o)
    message("Running diffSpliceDGE for DEU-edgeR ...")
    res <- DEU_res(DEU_fit, mat$design, pair, "DEU", fdr_cutoff, mode)  

    message("Running diffSpliceDGE for DEJU-edgeR ...")
    res <- DEU_res(DEJU_fit, mat$design, pair, "DEJU", fdr_cutoff, mode)  
    setwd(wd)
    
  }, BPPARAM = BPPARAM)
  
} else if (mode == "case_study") {
  
  message("Creating output folders to store read counts and DEU results ...")
  featureCounts_o <- paste0(DIR, "featureCounts/")
  dir.create(featureCounts_o, recursive = TRUE, showWarnings = FALSE)
  
  message("Constructing design matrix for contrast ...")
  mat <- construct_model_matrix(targetp)
  
  BAM_files <- paste0(DIR, 
                      "aligned_pass2/", 
                      mat$samples, "/", 
                      mat$samples, 
                      ".Aligned.sortedByCoord.out.bam")
  
  print(BAM_files)
  
  message("Checking if running featureCounts yet ...")
  f <- paste0(featureCounts_o, c("internal_exon_count.tsv", 
                                 "internal_exon_count_annotation.tsv", 
                                 "internal_exon_count_stat.tsv",
                                 "exon_count.tsv", 
                                 "exon_count_annotation.tsv", 
                                 "exon_count_stat.tsv",
                                 "junction_count.tsv"))
  print(any(!file.exists(f)))
  if (any(!file.exists(f))) {
    
    message("Quantifing internal exon read count ...")
    IE_count_table <- quantify_IE_count(BAM_files, isPairedEnd, flat_exon, featureCounts_o, wd)
    
    message("Quantifying exon-level read with duplicated junction count and junction read count ...")
    E_J_count_table <- quantify_E_J_count(BAM_files, isPairedEnd, flat_exon, FASTA, featureCounts_o, wd)
    
  }
  
  message("Constructing a matrix for exon read counts with duplicated junction count ...")
  E <- process_E_count_withDupJ(featureCounts_o)
  
  message("Start DEU analysis for DEU-edgeR for ...")
  DEU_fit <- DEU_analysis(E$E_count,
                          E$E_annot,
                          mat$group, mat$design, mode)
 
  message("Combining internal exon and junction read count matrix ...")  
  IE_J <- process_IE_J_count(featureCounts_o, SJ_database)

  message("Start DEJU analysis for DEU-edgeR ...")
  DEJU_fit <- DEU_analysis(IE_J$IE_J_count,
                           IE_J$IE_J_annot,
                           mat$group, mat$design, mode)
  
  pair <- strsplit(pair, ",")[[1]]
  for (p in pair) {
    
    DEU_analysis_o <- paste0(DIR, "edgeR_diffSpliceDGE/", p, "/")
    dir.create(DEU_analysis_o, recursive = TRUE, showWarnings = FALSE) 
    
    setwd(DEU_analysis_o)
    message("Running diffSpliceDGE for DEU-edgeR ...")
    res <- DEU_res(DEU_fit, mat$design, p, "DEU", fdr_cutoff, mode)  
    
    message("Running diffSpliceDGE for DEJU-edgeR ...")
    res <- DEU_res(DEJU_fit, mat$design, p, "DEJU", fdr_cutoff, mode)  
    setwd(wd)
  }

} else {
  
  print("Invalid mode specified. Choose 'simulation' or 'case_study'!")
  
}

###############################################################################
# OK <- requireNamespace("devtools", quietly = TRUE)
# if (!OK) {
#   stop("devtools package required but is not installed (or can't be loaded)")
# }
# library(Rsubread)
# library(edgeR)
# library(BiocParallel)
# 
# quantify_IE_count <- function(BAM_files, isPairedEnd, flat_exon, featureCounts_o, wd) {
#   
#   message("Running featureCounts with nonSplitOnly=TRUE and splitOnly=FALSE ...")
#   internal_exon_count <- Rsubread::featureCounts(BAM_files,
#                                                  annot.ext=flat_exon,
#                                                  tmpDir=featureCounts_o,
#                                                  nonSplitOnly=TRUE,
#                                                  splitOnly=FALSE,
#                                                  useMetaFeatures=FALSE,
#                                                  allowMultiOverlap=TRUE,
#                                                  isPairedEnd=isPairedEnd,
#                                                  nthreads=8,
#                                                  minMQS=255)
#   
#   message("Saving internal exon read counts and annotations ...")
#   setwd(featureCounts_o)
#   write.table(internal_exon_count$counts, 
#               "internal_exon_count.tsv", 
#               sep="\t", quote=FALSE, row.names=FALSE)
#   write.table(internal_exon_count$annotation, 
#               "internal_exon_count_annotation.tsv", 
#               sep="\t", quote=FALSE, row.names=FALSE)
#   write.table(internal_exon_count$stat, 
#               "internal_exon_count_stat.tsv", 
#               sep="\t", quote=FALSE, row.names=FALSE)
#   setwd(wd)
# }
# 
# quantify_E_J_count <- function(BAM_files, isPairedEnd, flat_exon, FASTA, featureCounts_o, wd) {
#   
#   message("Running featureCounts with juncCounts=TRUE ...")
#   junction_count <- Rsubread::featureCounts(BAM_files,
#                                             annot.ext=flat_exon,
#                                             tmpDir=featureCounts_o,
#                                             genome=FASTA,
#                                             juncCounts=TRUE,
#                                             useMetaFeatures=FALSE,
#                                             allowMultiOverlap=TRUE,
#                                             isPairedEnd=isPairedEnd,
#                                             nthreads=8,
#                                             minMQS=255)
#   
#   message("Saving exon-level read counts and annotations with duplicated junction read counts ...")
#   setwd(featureCounts_o)
#   write.table(junction_count$counts, 
#               "exon_count.tsv", 
#               sep="\t", quote=FALSE, row.names=FALSE)
#   write.table(junction_count$annotation, 
#               "exon_count_annotation.tsv", 
#               sep="\t", quote=FALSE, row.names=FALSE)
#   write.table(junction_count$stat, 
#               "exon_count_stat.tsv", 
#               sep="\t", quote=FALSE, row.names=FALSE)
#   
#   message("Saving junction read counts and annotations ...")  
#   write.table(junction_count$counts_junction[!is.na(junction_count$counts_junction$PrimaryGene),], 
#               "junction_count.tsv", 
#               sep="\t", quote=FALSE, row.names=FALSE)
#   setwd(wd)
#   
# }
# 
# filter_gene <- function(y, REF) {
#   
#   message("Removing genes without gene symbol ...")
#   gene_info_p <- paste0(REF, "gene_info_GRCm39_M32_ensembl109.tsv")
#   gene_info <- read.table(gene_info_p, header=T, na.strings="." )
#   m <- match(y$genes$GeneID, gene_info$GeneID)
#   y$genes$Symbol <- gene_info$Symbol[m]
#   keep <- !is.na(y$genes$Symbol)
#   table(keep)
#   y <- y[keep, , keep.lib.sizes=FALSE]  
#   return(y)
# }
# 
# construct_model_matrix <- function(targetp, p) {
#   
#   message("Reading target ...")
#   target <- read.table(targetp, sep="\t", header=TRUE, stringsAsFactors=FALSE)
#   group <- factor(target$group)
#   samples <- target$sampleID
#   
#   message("Constructing model matrix ...")
#   design <- model.matrix(~ 0 + group)
#   colnames(design) <- gsub("group", "", colnames(design))
#   
#   message("Making contrast ...")
#   contr <- limma::makeContrasts(p, levels = design)
#   
#   return(list('group'=group, 'samples'=samples, 'design'=design, 'contr'=contr))
# }
# 
# process_E_count_withDupJ <- function(featureCounts_o) {
#   
#   E_count <- read.table(paste0(featureCounts_o, "exon_count.tsv"), header=TRUE)
#   E_count <- data.frame(E_count[,1:ncol(E_count)], row.names=NULL, check.names = FALSE)
#   E_annot <- read.table(paste0(featureCounts_o, "exon_count_annotation.tsv"), header=TRUE)
#   
#   return(list('E_count'=E_count,
#               'E_annot'=E_annot))
# }
# 
# reassign_SJ <- function(J_count, SJ_database) {
#   
#   message("Only consider unique junctions ...")
#   uniq_SJ <- SJ_database[SJ_database$freq==1,]
#   m1 <- match(J_count$juncID, uniq_SJ$juncID)
#   J_count$PrimaryGene <- ifelse(!is.na(m1), uniq_SJ$geneID[m1], J_count$PrimaryGene)
#   
#   message("Whether junctions are annotated or not ...")
#   m2 <- match(J_count$juncID, SJ_database$juncID)
#   J_count$annotated <- ifelse(!is.na(m2), 1, 0)
#   
#   return(J_count) 
# }
# 
# process_IE_J_count <- function(featureCounts_o, SJ_database) {
#   
#   message("Analyzing internal exon read counts ...")
#   IE_count <- read.table(paste0(featureCounts_o, "internal_exon_count.tsv"), header=TRUE)
#   IE_count <- data.frame(IE_count[,1:ncol(IE_count)], row.names=NULL, check.names = FALSE)
#   IE_annot <- read.table(paste0(featureCounts_o, "internal_exon_count_annotation.tsv"), header=TRUE)
#   IE_annot <- cbind(IE_annot, Region="Exon", annotated=1)
#   
#   message("Analyzing junction read counts ...")
#   J_count <- read.table(paste0(featureCounts_o, "junction_count.tsv"), header=TRUE)
#   J_count$juncID <- paste(J_count$Site1_chr, 
#                           J_count$Site1_location, 
#                           J_count$Site2_location, 
#                           sep="_")
#   
#   message("Re-assigning junctions to correct genes using junction database ...")
#   J_count <- reassign_SJ(J_count, SJ_database)
#   
#   message("Processing final junction count table ...")
#   J_annot <- data.frame(
#     GeneID=J_count$PrimaryGene,
#     Chr=J_count$Site1_chr,
#     Start=J_count$Site1_location,
#     End=J_count$Site2_location
#   )
#   m <- match(J_annot$GeneID, IE_annot$GeneID)
#   Strand <- IE_annot$Strand[m]
#   J_annot <- cbind(J_annot, Strand=Strand, Length=1, Region="Junction", annotated=J_count$annotated)
#   
#   J_count <- data.frame(J_count[, 9:(ncol(J_count)-2)], row.names = NULL, check.names = FALSE)
#   
#   message("Combining internal exon and junction read count table ...")  
#   IE_J_count <- rbind(IE_count, J_count)
#   IE_J_annot <- rbind(IE_annot, J_annot)
#   
#   return(list('IE_J_count'=IE_J_count,
#               'IE_J_annot'=IE_J_annot))
# }
# 
# DEU_analysis <- function(r_count, annot, group, design, contr, mode, fdr_cutoff) {
#   
#   message("Constructing DGElist object ...")  
#   y <- DGEList(counts=r_count, genes=annot, group=group)
#   colnames(y) <- gsub("[.].*$", "", colnames(y))
#   
#   message("Filtering genes without gene symbol for real data ...")
#   if (mode == "case_study") y <- filter_gene(y, REF)
#   
#   message("Filtering exons with low mapping reads ...")
#   keep <- filterByExpr(y, group=group)
#   y <- y[keep, , keep.lib.sizes=FALSE]
#   
#   message("Normalizing lib sizes ...")
#   y <- normLibSizes(y)
#   
#   message("Estimating dispersion ...")
#   y <- estimateDisp(y, design, robust=TRUE)
#   
#   message("Fitting GLM model for design matrix ...")
#   fit <- glmQLFit(y, design, robust=TRUE)
#   
#   message("Running diffSpliceDGE ..")
#   sp <- diffSpliceDGE(fit, contrast=contr, geneid="GeneID", exonid="Start")
#   DEU_simes <- topSpliceDGE(sp, test="Simes", number=Inf)
#   DEU_F <- topSpliceDGE(sp, test="gene", number=Inf)
#   DEU_exon <- topSpliceDGE(sp, test="exon", number=Inf)
#   
#   message("Extracting significant DEU genes - FDR cutoff 0.05 ...")  
#   sigDEU_simes <- DEU_simes[DEU_simes$FDR<=fdr_cutoff,]
#   sigDEU_F <- DEU_F[DEU_F$FDR<=fdr_cutoff,]
#   sigDEU_exon <- DEU_exon[DEU_exon$FDR<=fdr_cutoff,]
#   
#   return(list('simes'=DEU_simes,
#               'F'=DEU_F,
#               'exon'=DEU_exon,
#               'sig_simes'=sigDEU_simes,
#               'sig_F'=sigDEU_F,
#               'sig_exon'=sigDEU_exon))
# }
# 
# save_results <- function(DEU, DEJU, fdr_cutoff) {
#   
#   message("Saving all DEU results for DEU-edgeR")
#   write.table(DEU$simes, 
#               "DEU_simes_test.allGenes.tsv", 
#               sep="\t", quote=FALSE, row.names=FALSE)
#   write.table(DEU$F, 
#               "DEU_F_test.allGenes.tsv", 
#               sep="\t", quote=FALSE, row.names=FALSE)
#   write.table(DEU$exon, 
#               "DEU_exon_test.allGenes.tsv", 
#               sep="\t", quote=FALSE, row.names=FALSE)
#   
#   message("Saving significant DEU results for DEU-edgeR ...")
#   write.table(DEU$sig_simes, 
#               paste0("DEU_simes_test.sigGenes_", fdr_cutoff, ".tsv"), 
#               sep="\t", 
#               quote=FALSE, row.names=FALSE)  
#   write.table(DEU$sig_F, 
#               paste0("DEU_F_test.sigGenes_", fdr_cutoff, ".tsv"), 
#               sep="\t", quote=FALSE, row.names=FALSE)  
#   write.table(DEU$sig_exon, 
#               paste0("DEU_exon_test.sigGenes_", fdr_cutoff, ".tsv"), 
#               sep="\t", quote=FALSE, row.names=FALSE)  
#   
#   message("Saving all DEU results for DEJU-edgeR ...")
#   write.table(DEJU$simes, 
#               "DEJU_simes_test.allGenes.tsv", 
#               sep="\t", quote=FALSE, row.names=FALSE) 
#   write.table(DEJU$F, 
#               "DEJU_F_test.allGenes.tsv", 
#               sep="\t", quote=FALSE, row.names=FALSE) 
#   write.table(DEJU$exon, 
#               "DEJU_exon_test.allGenes.tsv", 
#               sep="\t", quote=FALSE, row.names=FALSE) 
#   
#   message("Saving significant DEU results for DEJU-edgeR...")  
#   write.table(DEJU$sig_simes, 
#               paste0("DEJU_simes_test.sigGenes_", fdr_cutoff, ".tsv"), 
#               sep="\t", quote=FALSE, row.names=FALSE)  
#   write.table(DEJU$sig_F, 
#               paste0("DEJU_F_test.sigGenes_", fdr_cutoff, ".tsv"), 
#               sep="\t", quote=FALSE, row.names=FALSE)  
#   write.table(DEJU$sig_exon, 
#               paste0("DEJU_exon_test.sigGenes_", fdr_cutoff, ".tsv"), 
#               sep="\t", quote=FALSE, row.names=FALSE)
#   
# }
# 
# save_csv <- function(DEU, DEJU, fdr_cutoff) {
#   message("Saving all DEU results for DEU-edgeR ...")
#   write.csv(DEU$simes, 
#             "DEU_simes_test.allGenes.csv", 
#             row.names=FALSE)
#   write.csv(DEU$F, 
#             "DEU_F_test.allGenes.csv", 
#             row.names=FALSE)
#   write.csv(DEU$exon, 
#             "DEU_exon_test.allGenes.csv", 
#             row.names=FALSE)
# 
#   message("Saving significant DEU results for DEU-edgeR ...")
#   write.csv(DEU$sig_simes, 
#             paste0("DEU_simes_test.sigGenes_", fdr_cutoff, ".csv"), 
#             row.names=FALSE)
#   write.csv(DEU$sig_F, 
#             paste0("DEU_F_test.sigGenes_", fdr_cutoff, ".csv"), 
#             row.names=FALSE)
#   write.csv(DEU$sig_exon, 
#             paste0("DEU_exon_test.sigGenes_", fdr_cutoff, ".csv"), 
#             row.names=FALSE)
#   
#   message("Saving all DEU results for DEJU-edgeR ...")  
#   write.csv(DEJU$simes, 
#             "DEJU_simes_test.allGenes.csv", 
#             row.names=FALSE)
#   write.csv(DEJU$F, 
#             "DEJU_F_test.allGenes.csv", 
#             row.names=FALSE)
#   write.csv(DEJU$exon, 
#             "DEJU_exon_test.allGenes.csv", 
#             row.names=FALSE)
#   
#   message("Saving significant DEU results for DEJU-edgeR  ...")  
#   write.csv(DEJU$sig_simes, 
#             paste0("DEJU_simes_test.sigGenes_", fdr_cutoff, ".csv"), 
#             row.names=FALSE)
#   write.csv(DEJU$sig_F, 
#             paste0("DEJU_F_test.sigGenes_", fdr_cutoff, ".csv"), 
#             row.names=FALSE)
#   write.csv(DEJU$sig_exon, 
#             paste0("DEJU_exon_test.sigGenes_", fdr_cutoff, ".csv"), 
#             row.names=FALSE)
# }
# 
# # Getting parameters
# args <- R.utils::commandArgs(asValues = TRUE)
# print(args)
# 
# wd <- getwd()
# DIR <- as.character(args[['DIR']])
# REF <- as.character(args[['REF']])
# targetp <- as.character(args[['target']])
# p <- as.character(args[['pair']])
# mode <- as.character(args[['mode']])
# isPairedEnd <- as.logical(args[['isPairedEnd']])
# fdr_cutoff <- as.numeric(args[['fdr_cutoff']])
# 
# # DIR <- "/vast/projects/Spatial/tam/Differential_splicing/github/data/simulation/unbalanced_3_75_3/"
# # REF <- "/vast/projects/Spatial/tam/Differential_splicing/github/annotation/"
# # targetp <- "/vast/projects/Spatial/tam/Differential_splicing/github/data/simulation/target/target.3vs3.tsv"
# # p <- "Group_2-Group_1"
# # mode <- "simulation"
# # isPairedEnd <- TRUE
# # fdr_cutoff <- 0.05
# # noOfSim <- 1
# # seed <- 2024
# # workers <- 4
# 
# FASTA <- paste0(REF, "GRCm39.primary_assembly.genome.fa.gz")
# flat_exon <- paste0(REF, "gencode.vM32.annotation.flattened.exon.saf")
# SJ <- paste0(REF, "gencode.vM32.primary_assembly.annotation.SJdatabase.tsv")
# SJ_database <- read.table(SJ, header=TRUE)
# 
# message("Check parameters ...")
# print(wd)
# print(DIR)
# print(REF)
# print(targetp)
# print(FASTA)
# print(flat_exon)
# print(SJ)
# print(p)
# print(mode)
# print(isPairedEnd)
# print(fdr_cutoff)
# 
# if (mode == "simulation") {
#   
#   seed <- as.integer(args[['seed']])
#   workers <- as.integer(args[['workers']])
#   noOfSim <- as.integer(args[['noOfSim']])
#   simID <- paste0("S", 1:noOfSim)
#   
#   print(simID)
#   print(seed)
#   print(workers)
#   
#   ### 01. Register BPPARAM
#   message('Register BPPARAM')
#   BPPARAM <- MulticoreParam(workers = workers, progressbar = TRUE, RNGseed = seed)
#   register(BPPARAM = BPPARAM)
#   
#   ### 02. Run edgeR
#   message('Start running edgeR ...')
#   bplapply(simID, FUN = function(i){
#     
#     message(paste0('Processing: ', i, " ..."))
#     
#     message("Creating output folders to store read counts and DEU results ...")
#     featureCounts_o <- paste0(DIR, "featureCounts/", i, "/")
#     dir.create(featureCounts_o, showWarnings = FALSE, recursive = TRUE)
#     DEU_analysis_o <- paste0(DIR, "edgeR_diffSpliceDGE/", i, "/")
#     dir.create(DEU_analysis_o, showWarnings = FALSE, recursive = TRUE)
#     
#     message("Constructing design matrix for contrast ...")
#     mat <- construct_model_matrix(targetp, p)
#     
#     BAM_files <- paste0(DIR, 
#                        "aligned_pass2/", 
#                        i, "/", 
#                        mat$samples, "/", 
#                        mat$samples, 
#                        ".Aligned.sortedByCoord.out.bam")
#     
#     print(BAM_files)
#     
#     message("Checking if running featureCounts yet ...")
#     f <- paste0(featureCounts_o, c("internal_exon_count.tsv", 
#                                    "internal_exon_count_annotation.tsv", 
#                                    "internal_exon_count_stat.tsv",
#                                    "exon_count.tsv", 
#                                    "exon_count_annotation.tsv", 
#                                    "exon_count_stat.tsv",
#                                    "junction_count.tsv"))
#     print(any(!file.exists(f)))
#     if (any(!file.exists(f))) {
#       
#       message("Quantifing internal exon read count ...")
#       IE_count_table <- quantify_IE_count(BAM_files, isPairedEnd, flat_exon, featureCounts_o, wd)
#       
#       message("Quantifying exon-level read with duplicated junction count and junction read count ...")
#       E_J_count_table <- quantify_E_J_count(BAM_files, isPairedEnd, flat_exon, FASTA, featureCounts_o, wd)
#       
#     }
#     
#     message("Constructing a matrix for exon read counts with duplicated junction count ...")
#     E <- process_E_count_withDupJ(featureCounts_o)
#     
#     message("Starting DEU analysis for DEU-edgeR method ...")
#     DEU <- DEU_analysis(E$E_count, 
#                         E$E_annot, 
#                         mat$group, mat$design, mat$contr, mode, fdr_cutoff)
#     
#     message("Combining internal exon and junction read count matrix ...")  
#     IE_J <- process_IE_J_count(featureCounts_o, SJ_database)
#     
#     message("Starting DEJU analysis for DEJU-edgeR ...")
#     DEJU <- DEU_analysis(IE_J$IE_J_count, 
#                          IE_J$IE_J_annot, 
#                          mat$group, mat$design, mat$contr, mode, fdr_cutoff)
#     
#     setwd(DEU_analysis_o)
#     message("Saving final DEU results ...")
#     res <- save_results(DEU, DEJU, fdr_cutoff)
#     setwd(wd)
#     
#   }, BPPARAM = BPPARAM)
#   
# } else if (mode == "case_study") {
#   
#   message("Creating output folders to store read counts and DEU results ...")
#   featureCounts_o <- paste0(DIR, "featureCounts/")
#   dir.create(featureCounts_o, recursive = TRUE, showWarnings = FALSE)
#   DEU_analysis_o <- paste0(DIR, "edgeR_diffSpliceDGE/", p, "/")
#   dir.create(DEU_analysis_o, recursive = TRUE, showWarnings = FALSE)
#   
#   message("Constructing design matrix for contrast ...")
#   mat <- construct_model_matrix(targetp, p)
#   
#   BAM_files <- paste0(DIR, 
#                       "aligned_pass2/", 
#                       mat$samples, "/", 
#                       mat$samples, 
#                       ".Aligned.sortedByCoord.out.bam")
#   
#   print(BAM_files)
#   
#   message("Checking if running featureCounts yet ...")
#   f <- paste0(featureCounts_o, c("internal_exon_count.tsv", 
#                                  "internal_exon_count_annotation.tsv", 
#                                  "internal_exon_count_stat.tsv",
#                                  "exon_count.tsv", 
#                                  "exon_count_annotation.tsv", 
#                                  "exon_count_stat.tsv",
#                                  "junction_count.tsv"))
#   print(any(!file.exists(f)))
#   if (any(!file.exists(f))) {
#     
#     message("Quantifing internal exon read count ...")
#     IE_count_table <- quantify_IE_count(BAM_files, isPairedEnd, flat_exon, featureCounts_o, wd)
#     
#     message("Quantifying exon-level read with duplicated junction count and junction read count ...")
#     E_J_count_table <- quantify_E_J_count(BAM_files, isPairedEnd, flat_exon, FASTA, featureCounts_o, wd)
#     
#   }
#   
#   message("Constructing a matrix for exon read counts with duplicated junction count ...")
#   E <- process_E_count_withDupJ(featureCounts_o)
#   
#   message("Starting DEU analysis for DEU-edgeR method ...")
#   DEU <- DEU_analysis(E$E_count, 
#                       E$E_annot, 
#                       mat$group, mat$design, mat$contr, mode, fdr_cutoff)
#   
#   message("Combining internal exon and junction read count matrix ...")  
#   IE_J <- process_IE_J_count(featureCounts_o, SJ_database)
#   
#   message("Starting DEJU analysis for DEJU-edgeR ...")
#   DEJU <- DEU_analysis(IE_J$IE_J_count, 
#                        IE_J$IE_J_annot, 
#                        mat$group, mat$design, mat$contr, mode, fdr_cutoff)
#   
#   setwd(DEU_analysis_o)
#   message("Saving final DEU results ...")
#   res <- save_results(DEU, DEJU, fdr_cutoff)
#   setwd(wd)
#   
# } else {
#   
#   print("Invalid mode specified. Choose 'simulation' or 'case_study'!")
#   
# }

###############################################################################
# OK <- requireNamespace("devtools", quietly = TRUE)
# if (!OK) {
#   stop("devtools package required but is not installed (or can't be loaded)")
# }
# library(Rsubread)
# library(edgeR)
# library(BiocParallel)
# 
# # Getting parameters
# args <- R.utils::commandArgs(asValues = TRUE)
# print(args)
# 
# DIR <- as.character(args[['DIR']])
# REF <- as.character(args[['REF']])
# targetp <- as.character(args[['target']])
# 
# # DIR <- "/vast/projects/Spatial/tam/Differential_splicing/github/data/simulation/unbalanced_3_75_3/"
# # REF <- "/vast/projects/Spatial/tam/Differential_splicing/github/annotation/"
# # targetp <- "/vast/projects/lab_chen/tam/Differential_splicing/simulation_3vs3_mouse_genome/RNA-seq/config/target.3vs3.tsv"
# # path_to_sj <- paste0(REF, "gencode.vM32.primary_assembly.annotation.uniqSJInfo.tsv")
# # num_of_occur_sj <- paste0(REF, "gencode.vM32.primary_assembly.annotation.uniqSJInfo.numOfOccur.tsv")
# # annotated_sj <- read.table(path_to_sj, header=TRUE)
# # numOfOccur_sj <- read.table(num_of_occur_sj, header=TRUE)
# # noOfSim <- 2
# # seed <- 2024
# # workers <- 90
# 
# FASTA <- paste0(REF, "GRCm39.primary_assembly.genome.fa.gz")
# flatExonAnno <- paste0(REF, "gencode.vM32.annotation.flattened.exon.saf")
# SJ <- paste0(REF, "gencode.vM32.primary_assembly.annotation.SJdatabase.tsv")
# 
# SJ_database <- read.table(SJ, header=TRUE)
# 
# if (mode == "simulation") {
# 
#   seed <- as.integer(args[['seed']])
#   workers <- as.integer(args[['workers']])
#   noOfSim <- as.integer(args[['noOfSim']])
#   simID <- paste0("S", 1:noOfSim)
#   
#   ### 0. Check parameter
#   message("Check valid parameters")
#   print(DIR)
#   print(REF)
#   print(FASTA)
#   print(flatExonAnno)
#   print(path_to_sj)
#   print(num_of_occur_sj)
#   print(targetp)
#   print(simID)
#   print(seed)
#   print(workers)
#   
#   ### 01. Register BPPARAM
#   message('Register BPPARAM')
#   BPPARAM <- MulticoreParam(workers = workers, progressbar = TRUE, RNGseed = seed)
#   register(BPPARAM = BPPARAM)
#   
#   ### 02. Run edgeR
#   message('Start running edgeR ...')
#   bplapply(simID, FUN = function(i){
#     
#     message(paste0('Processing: ', i, " ..."))
#     OUTPUT <- paste0(DIR, "edgeR/", i, "/")
#     readCount_output <- paste0(DIR, "featureCounts/", i, "/")
#     dir.create(OUTPUT, showWarnings = FALSE, recursive = TRUE)
#     dir.create(readCount_output, showWarnings = FALSE, recursive = TRUE)
#     
#     ### 021. Read count quantification
#     message('Quantifying read counts at exon, junction ...')
#     target <- read.table(targetp, sep="\t", header=TRUE, stringsAsFactors=FALSE)
#     group <- factor(target$Group, levels=c("Group_1","Group_2"))
#     design <- model.matrix(~ group)
#     PE_bam_files <- paste0(DIR, 
#                            "STAR_aligned_pass2_minUniqSJReads_3/", 
#                            i, "/", 
#                            target$Sample, "/", 
#                            target$Sample, 
#                            ".Aligned.sortedByCoord.out.bam")
#     
#     exon_count <- Rsubread::featureCounts(PE_bam_files, annot.ext=flatExonAnno,
#                                           tmpDir=readCount_output,
#                                           nonSplitOnly=TRUE, splitOnly=FALSE,
#                                           useMetaFeatures=FALSE, allowMultiOverlap=TRUE,
#                                           isPairedEnd=TRUE, nthreads=8, minMQS=255)
#     write.table(exon_count$counts, "internal_exon_count.tsv", 
#                 sep="\t", quote=FALSE, row.names=FALSE)
#     write.table(exon_count$annotation, "internal_exon_count_annotation.tsv", 
#                 sep="\t", quote=FALSE, row.names=FALSE)
#     write.table(exon_count$stat, "internal_exon_count_stat.tsv", 
#                 sep="\t", quote=FALSE, row.names=FALSE)
#     
#     
#     junction_count <- Rsubread::featureCounts(PE_bam_files, annot.ext=flatExonAnno,
#                                               tmpDir=readCount_output,
#                                               genome=FASTA, juncCounts=TRUE,
#                                               useMetaFeatures=FALSE, allowMultiOverlap=TRUE,
#                                               isPairedEnd=TRUE, nthreads=8, minMQS=255)
#     write.table(junction_count$counts, "exon_count.tsv", 
#                 sep="\t", quote=FALSE, row.names=FALSE)
#     write.table(junction_count$counts_junction[!is.na(junction_count$counts_junction$PrimaryGene),], 
#                 "junction_count.tsv", 
#                 sep="\t", quote=FALSE, row.names=FALSE)
#     write.table(junction_count$annotation, "exon_count_annotation.tsv", 
#                 sep="\t", quote=FALSE, row.names=FALSE)
#     write.table(junction_count$stat, "exon_count_stat.tsv", 
#                 sep="\t", quote=FALSE, row.names=FALSE)
#     
#     ### 022. Loading read counts
#     message('Loading read counts ...')
#     
#     # Exon-level read counts with double-counted junctions
#     exon_count_f <- read.table("exon_count.tsv", header=TRUE)
#     exon_count <- data.frame(exon_count_f[,1:ncol(exon_count_f)], row.names=NULL, check.names = FALSE)
#     exon_annot <- read.table("exon_count_annotation.tsv", header=TRUE)
#     
#     # Internal exon read counts
#     IE_count_f <- read.table("internal_exon_count.tsv", header=TRUE)
#     IE_count <- data.frame(IE_count_f[,1:ncol(IE_count_f)], row.names=NULL, check.names = FALSE)
#     IE_annot <- read.table("internal_exon_count_annotation.tsv", header=TRUE)
#     IE_annot <- cbind(IE_annot, Region="Exon", annotated=1)
#     
#     # Junction read counts
#     J_count <- read.table("junction_count.tsv", header=TRUE)
#     J_count$juncID <- paste(J_count$Site1_chr, 
#                               J_count$Site1_location, 
#                               J_count$Site2_location, 
#                               sep="_")
#     m_numOfOccur <- match(annotated_sj$juncID, numOfOccur_sj$juncID)
#     annotated_sj$numOfOccur <- numOfOccur_sj$numOfOccur[m_numOfOccur]
#     annotated_sj_numOfOccur_1 <- annotated_sj[annotated_sj$numOfOccur==1,]
#     m1 <- match(J_count$juncID, annotated_sj_numOfOccur_1$juncID)
#     J_count$PrimaryGene <- ifelse(!is.na(m1), annotated_sj_numOfOccur_1$geneID[m1], J_count$PrimaryGene)
#     m2 <- match(J_count$juncID, annotated_sj$juncID)
#     J_count$annotated <- ifelse(!is.na(m2), 1, 0)
#     J_count <- data.frame(J_count[, 9:(ncol(J_count_f)-2)], row.names = NULL, check.names = FALSE)
#     J_annot <- data.frame(
#       GeneID=J_count_f$PrimaryGene,
#       Chr=J_count_f$Site1_chr,
#       Start=J_count_f$Site1_location,
#       End=J_count_f$Site2_location
#     )
#     m <- match(J_annot$GeneID, IE_annot$GeneID)
#     Strand <- IE_annot$Strand[m]
#     J_annot <- cbind(J_annot, Strand=Strand, Length=1, Region="Junction", annotated=J_count_f$annotated)
#     
#     IE_J_count <- rbind(IE_count, J_count)
#     IE_J_annot <- rbind(IE_annot, J_annot)
#     
#     ### 023. Preliminary analysis
#     message('Preliminary analysis ...')
#     
#     ## Existing approach
#     y_DEU <- DGEList(counts=exon_count, genes=exon_annot, group=group)
#     colnames(y_DEU) <- gsub("[.].*$", "", colnames(y_DEU))
#     keep <- filterByExpr(y_DEU, group=group)
#     y_DEU <- y_DEU[keep, , keep.lib.sizes=FALSE]
#     y_DEU <- normLibSizes(y_DEU)
#     y_DEU <- estimateDisp(y_DEU, design, robust=TRUE)
#     print(y_DEU)
#     fit_DEU <- glmQLFit(y_DEU, design, robust=TRUE)
#     
#     ## Junction approach (internal exon read + junction read)
#     y_DEJU <- DGEList(counts=internal_exon_junc_count, genes=internal_exon_junc_annot, group=group)
#     colnames(y_DEJU) <- gsub("[.].*$", "", colnames(y_DEJU))
#     keep <- filterByExpr(y_DEJU, group=group)
#     y_DEJU <- y_DEJU[keep, , keep.lib.sizes=FALSE]
#     y_DEJU <- normLibSizes(y_DEJU)
#     y_DEJU <- estimateDisp(y_DEJU, design, robust=TRUE)
#     print(y_DEJU)
#     fit_DEJU <- glmQLFit(y_DEJU, design, robust=TRUE)
#     
#     ### 024. Differential alternative splicing analysis
#     message('DEU analysis ...')
#     
#     message('Save diffSpliceDGE results ...')
#     
#     ## Existing approach
#     sp_DEU <- diffSpliceDGE(fit_DEU, coef=2, geneid="GeneID", exonid="Start")
#     DEU_simes <- topSpliceDGE(sp_DEU, test="Simes", n=Inf)
#     DEU_gene <- topSpliceDGE(sp_DEU, test="gene", n=Inf)
#     DEU_exon <- topSpliceDGE(sp_DEU, test="exon", n=Inf)
#     write.table(DEU_simes, "DEU_simes_test.allGenes.tsv",
#                 row.names=FALSE, sep="\t", quote=FALSE)
#     write.table(DEU_gene, "DEU_F_test.allGenes.tsv", 
#                 row.names=FALSE, sep="\t", quote=FALSE)
#     write.table(DEU_exon, "DEU_exon_test.allGenes.tsv", 
#                 row.names=FALSE, sep="\t", quote=FALSE)
#     
#     ## Junction approach
#     sp_DEJU <- diffSpliceDGE(fit_DEJU, coef=2, geneid="GeneID", exonid="Start")
#     DEJU_simes <- topSpliceDGE(sp_DEJU, test="Simes", n=Inf)
#     DEJU_gene <- topSpliceDGE(sp_DEJU, test="gene", n=Inf)
#     DEJU_exon <- topSpliceDGE(sp_DEJU, test="exon", n=Inf)
#     write.table(DEJU_simes, "DEJU_simes_test.allGenes.tsv", 
#                 row.names=FALSE, sep="\t", quote=FALSE)
#     write.table(DEJU_gene, "DEJU_simes_test.allGenes.tsv", 
#                 row.names=FALSE, sep="\t", quote=FALSE)
#     write.table(DEJU_exon, "DEJU_simes_test.allGenes.tsv", 
#                 row.names=FALSE, sep="\t", quote=FALSE)
#     
#     message('Save significant DEUs ...')  
#     ## Existing approach
#     DEU_simes <- topSpliceDGE(sp_DEU, test="Simes", n=NULL, FDR=0.05)
#     DEU_gene <- topSpliceDGE(sp_DEU, test="gene", n=NULL, FDR=0.05)
#     write.table(DEU_simes, "DEU_simes_test.sigGenes_0.05.tsv", 
#                 row.names=FALSE, sep="\t", quote=FALSE)
#     write.table(DEU_gene, "DEU_F_test.sigGenes_0.05.tsv", 
#                 row.names=FALSE, sep="\t", quote=FALSE)
#     
#     ## Junction approach
#     DEJU_simes <- topSpliceDGE(sp_DEJU, test="Simes", n=NULL, FDR=0.05)
#     DEJU_gene <- topSpliceDGE(sp_DEJU, test="gene", n=NULL, FDR=0.05)
#     write.table(DEJU_simes, "DEJU_simes_test.sigGenes_0.05.tsv", 
#                 row.names=FALSE, sep="\t", quote=FALSE)
#     write.table(DEJU_gene, "DEJU_F_test.sigGenes_0.05.tsv", 
#                 row.names=FALSE, sep="\t", quote=FALSE)
#     
#   }, BPPARAM = BPPARAM)
#   
# } else if (mode == "case_study") {
# 
#   message("Loading groups to be compared ...")
#   p <- as.integer(args[['p']])
#   # p <- "Basal_LP"
#   
#   message("Creating output folders to store read counts and DEU results ...")
#   readCount_output <- paste0(DIR, "STAR_output/featureCounts/")
#   OUTPUT <- paste0(DIR, "STAR_output/edgeR/", p, "/")
#   
#   dir.create(OUTPUT, recursive = TRUE, showWarnings = FALSE)
#   dir.create(readCount_output, recursive = TRUE, showWarnings = FALSE)
# 
#   PE_bam_files <- list.files(path=paste0(DIR, "/aligned_pass2"), pattern="*.bam", full.names=TRUE, recursive=TRUE)
#   
#   message("Quantifing internal exon read count ...")
#   internal_exon_count <- featureCounts(PE_bam_files,
#                                        annot.ext=flatExonAnno,
#                                        tmpDir=readCount_output,,
#                                        nonSplitOnly=TRUE, 
#                                        splitOnly=FALSE,
#                                        useMetaFeatures=FALSE, 
#                                        allowMultiOverlap=TRUE,
#                                        isPairedEnd=TRUE, 
#                                        nthreads=8, 
#                                        minMQS=255)
#   
#   write.table(internal_exon_count$counts, 
#               "internal_exon_count.tsv", sep="\t", quote=FALSE, row.names=FALSE)
#   write.table(internal_exon_count$annotation, 
#               "internal_exon_count_annotation.tsv", sep="\t", quote=FALSE, row.names=FALSE)
#   write.table(internal_exon_count$stat, 
#               "internal_exon_count_stat.tsv", sep="\t", quote=FALSE, row.names=FALSE)
#   
#   message("Quantifying exon-level read with duplicated junction count and junction read count ...")
#   junction_count <- featureCounts(PE_bam_files, 
#                                   annot.ext=flatExonAnno,
#                                   tmpDir=readCount_output,
#                                   genome=path_to_ref, 
#                                   juncCounts=TRUE,
#                                   useMetaFeatures=FALSE, 
#                                   allowMultiOverlap=TRUE,
#                                   isPairedEnd=TRUE, 
#                                   nthreads=8, 
#                                   minMQS=255)
#   
#   write.table(junction_count$counts, 
#               "exon_count.tsv", sep="\t", quote=FALSE, row.names=FALSE)
#   write.table(junction_count$counts_junction[!is.na(junction_count$counts_junction$PrimaryGene),], 
#               "junction_count.tsv", sep="\t", quote=FALSE, row.names=FALSE)
#   write.table(junction_count$annotation, 
#               "exon_count_annotation.tsv", sep="\t", quote=FALSE, row.names=FALSE)
#   write.table(junction_count$stat, 
#               "exon_count_stat.tsv", sep="\t", quote=FALSE, row.names=FALSE)
# 
#   message("Loading gene annotation ...")
#   gene_info_p <- paste0(REF, "gene_info_GRCm39_M32_ensembl109.tsv")
#   gene_info <- read.table(gene_info_p, header=T, na.strings="." )
#   
#   message("Loading target info ...")
#   target <- read.table(targetp, sep="\t", header=TRUE, stringsAsFactors=FALSE)
#   group <- factor(target$Group)
#   
#   message("Constructing design matrix for contrast ...")
#   design <- model.matrix(~ 0 + group)
#   colnames(design) <- gsub("group", "", colnames(design))
#   g <- unlist(strsplit(p, "-"))
#   contr <- makeContrasts(g[1] - g[2], levels = design)  
#   
#   message("1. Starting DEU analysis for DEU-edgeR method ...")
#   message("Constructing a matrix for exon read counts with duplicated junction count ...")
#   E_count_f <- read.table("exon_count.tsv", header=TRUE)
#   E_count <- data.frame(E_count_f[,1:ncol(E_count_f)], row.names=NULL, check.names = FALSE)
#   E_annot <- read.table("exon_count_annotation.tsv", header=TRUE)
#   
#   message("Constructing DGElist object ...")  
#   y_DEU <- DGEList(counts=E_count, genes=E_annot, group=group)
#   colnames(y_DEU) <- gsub("[.].*$", "", colnames(y_DEU))
#   
#   message("Removing genes without gene symbol ...")
#   m <- match(y_DEU$genes$GeneID, gene_info$GeneID)
#   y_DEU$genes$Symbol <- gene_info$Symbol[m]
#   keep <- !is.na(y_DEU$genes$Symbol)
#   table(keep)
#   y_DEU <- y_DEU[keep, , keep.lib.sizes=FALSE]
#   
#   message("Filtering exons with low mapping reads ...")
#   keep <- filterByExpr(y_DEU, group=group)
#   table(keep)
#   y_DEU <- y_DEU[keep, , keep.lib.sizes=FALSE]
#   
#   message("Normalizing lib sizes ...")
#   y_DEU <- normLibSizes(y_DEU)
#   y_DEU$samples
#   
#   message("Estimating dispersion ...")
#   y_DEU <- estimateDisp(y_DEU, design, robust=TRUE)
#   y_DEU$common.dispersion
#   
#   message("Fitting GLM model for design matrix ...")
#   fit_DEU <- glmQLFit(y_DEU, design, robust=TRUE)
# 
#   message("Running diffSpliceDGE ..")
#   sp_DEU <- diffSpliceDGE(fit_DEU, contrast=contr, geneid="GeneID", exonid="Start")
#   DEU_simes <- topSpliceDGE(sp_DEU, test="Simes", number=Inf)
#   DEU_F <- topSpliceDGE(sp_DEU, test="gene", number=Inf)
#   DEU_exon <- topSpliceDGE(sp_DEU, test="exon", number=Inf)
# 
#   message("Saving all DEU results ...")  
#   write.csv(DEU_simes, paste0("DEU_simes_test.", p, ".allGenes.csv"), row.names=FALSE)
#   write.table(DEU_simes, paste0("DEU_simes_test.", p, ".allGenes.tsv"), 
#               sep="\t", quote=FALSE, row.names=FALSE)
#   write.csv(DEU_F, paste0("DEU_F_test.", p, ".allGenes.csv"), row.names=FALSE)
#   write.table(DEU_F, paste0("DEU_F_test.", p, ".allGenes.tsv"), 
#               sep="\t", quote=FALSE, row.names=FALSE)
#   write.csv(DEU_exon, paste0("DEU_exon_test.", p, ".allGenes.csv"), row.names=FALSE)
#   write.table(DEU_exon, paste0("DEU_exon_test.", p, ".allGenes.tsv"), 
#               sep="\t", quote=FALSE, row.names=FALSE)
# 
#   message("Extracting significant DEU genes - FDR cutoff 0.05 ...")  
#   DEU_simes_0.05 <- DEU_simes[DEU_simes$FDR<=0.05,]
#   DEU_F_0.05 <- DEU_F[DEU_F$FDR<=0.05,]
#   DEU_exon_0.05 <- DEU_exon[DEU_exon$FDR<=0.05,]
# 
#   message("Saving significant DEU results ...")
#   write.csv(DEU_simes_0.05, paste0("DEU_simes_test.", p, ".sigGenes_0.05.csv"), row.names=FALSE)
#   write.table(DEU_simes_0.05, paste0("DEU_simes_test.", p, ".sigGenes_0.05.tsv"), 
#               sep="\t", quote=FALSE, row.names=FALSE)  
#   
#   write.csv(DEU_F_0.05, paste0("DEU_F_test.", p, ".sigGenes_0.05.csv"), row.names=FALSE)
#   write.table(DEU_F_0.05, paste0("DEU_F_test.", p, ".sigGenes_0.05.tsv"), 
#               sep="\t", quote=FALSE, row.names=FALSE)  
#   
#   write.csv(DEU_exon_0.05, paste0("DEU_exon_test.", p, ".sigGenes_0.05.csv"), row.names=FALSE)
#   write.table(DEU_exon_0.05, paste0("DEU_exon_test.", p, ".sigGenes_0.05.tsv"), 
#               sep="\t", quote=FALSE, row.names=FALSE)  
#   
#   message("2. Starting DEJU analysis for DEJU-edgeR ...")  
#   message("Analyzing internal exon read counts ...")
#   IE_count_f <- read.table("internal_exon_count.tsv", header=TRUE)
#   IE_count <- data.frame(IE_count_f[,1:ncol(IE_count_f)], row.names=NULL, check.names = FALSE)
#   IE_annot <- read.table("internal_exon_count_annotation.tsv", header=TRUE)
#   IE_annot <- cbind(IE_annot, Region="Exon", annotated=1)
#   
#   message("Analyzing junction read counts ...")
#   J_count_f <- read.table("junction_count.tsv", header=TRUE)
#   J_count_f$juncID <- paste(J_count_f$Site1_chr, 
#                             J_count_f$Site1_location, 
#                             J_count_f$Site2_location, 
#                             sep="_")
#   
#   message("Re-assigning junctions to correct genes using junction database ...")
#   m_numOfOccur <- match(annotated_sj$juncID, numOfOccur_sj$juncID)
#   annotated_sj$numOfOccur <- numOfOccur_sj$numOfOccur[m_numOfOccur]
#   
#   message("Only consider unique junctions ...")
#   annotated_sj_numOfOccur_1 <- annotated_sj[annotated_sj$numOfOccur==1,]
#   m1 <- match(junc_count_f$juncID, annotated_sj_numOfOccur_1$juncID)
#   J_count_f$PrimaryGene <- ifelse(!is.na(m1), annotated_sj_numOfOccur_1$geneID[m1], J_count_f$PrimaryGene)
#   m2 <- match(junc_count_f$juncID, annotated_sj$juncID)
#   J_count_f$annotated <- ifelse(!is.na(m2), 1, 0)
#   
#   J_count <- data.frame(J_count_f[, 9:(ncol(J_count_f)-2)], row.names = NULL, check.names = FALSE)
#   J_annot <- data.frame(
#     GeneID=J_count_f$PrimaryGene,
#     Chr=J_count_f$Site1_chr,
#     Start=J_count_f$Site1_location,
#     End=J_count_f$Site2_location
#   )
#   m <- match(J_annot$GeneID, IE_annot$GeneID)
#   Strand <- IE_annot$Strand[m]
#   junc_annot <- cbind(J_annot, Strand=Strand, Length=1, Region="Junction", annotated=J_count_f$annotated)
# 
#   message("Combining internal exon and junction read count matrix ...")  
#   IE_J_count <- rbind(IE_count, J_count)
#   IE_J_annot <- rbind(IE_annot, J_annot)
#   
#   message("Constructing DGElist object ...")
#   y_DEJU <- DGEList(counts=IE_J_count, genes=IE_J_annot, group=group)
#   colnames(y_DEJU) <- gsub("[.].*$", "", colnames(y_DEJU))
#   
#   message("Removing genes without gene symbol ...")
#   m <- match(y_DEJU$genes$GeneID, gene_info$GeneID)
#   y_DEJU$genes$Symbol <- gene_info$Symbol[m]
#   keep <- !is.na(y_DEJU$genes$Symbol)
#   table(keep)
#   y_DEJU <- y_DEJU[keep, , keep.lib.sizes=FALSE]
#   
#   message("Filtering exons/junctions with low mapping reads ...")
#   keep <- filterByExpr(y_DEJU, group=group)
#   y_DEJU <- y_DEJU[keep, , keep.lib.sizes=FALSE]
#   
#   message("Normalizing lib sizes ...")
#   y_DEJU <- normLibSizes(y_DEJU)
#   y_DEJU$samples
#   
#   message("Estimating dispersion ...")
#   y_DEJU <- estimateDisp(y_DEJU, design, robust=TRUE) # need "statmod" package
#   
#   message("Fitting GLM model for design matrix ...")
#   fit_DEJU <- glmQLFit(y_DEJU, design, robust=TRUE)
#   
#   message("Running diffSpliceDGE ..")
#   sp_DEJU <- diffSpliceDGE(fit_DEJU, contrast=contr, geneid="GeneID", exonid="Start")
#   DEJU_simes <- topSpliceDGE(sp_DEJU, test="Simes", number=Inf)
#   DEJU_F <- topSpliceDGE(sp_DEJU, test="gene", number=Inf)
#   DEJU_exon <- topSpliceDGE(sp_DEJU, test="exon", number=Inf)
# 
#   message("Saving all DEU results ...")
#   write.csv(DEJU_simes, paste0("DEJU_simes_test.", p, ".allGenes.csv"), row.names=FALSE)
#   write.table(DEJU_simes, paste0("DEJU_simes_test.", p, ".allGenes.tsv"), 
#               sep="\t", quote=FALSE, row.names=FALSE) 
#   write.csv(DEJU_F, paste0("DEJU_F_test.", p, ".allGenes.csv"), row.names=FALSE)
#   write.table(DEJU_F, paste0("DEJU_F_test.", p, ".allGenes.tsv"), 
#               sep="\t", quote=FALSE, row.names=FALSE) 
#   write.csv(DEJU_exon, paste0("DEJU_exon_test.", p, ".allGenes.csv"), row.names=FALSE)
#   write.table(DEJU_exon, paste0("DEJU_exon_test.", p, ".allGenes.tsv"), 
#               sep="\t", quote=FALSE, row.names=FALSE) 
# 
#   message("Extracting significant DEU genes - FDR cutoff 0.05 ...")
#   DEJU_simes_0.05 <- DEJU_simes[DEJU_simes$FDR<=0.05,]
#   DEJU_F_0.05 <- DEJU_F[DEJU_F$FDR<=0.05,]
#   DEJU_exon_0.05 <- DEJU_exon[DEJU_exon$FDR<=0.05,]
# 
#   message("Saving significant DEU results ...")  
#   write.csv(DEJU_simes_0.05, paste0("DEJU_simes_test.", p, ".sigGenes_0.05.csv"), row.names=FALSE)
#   write.table(DEJU_simes_0.05, paste0("DEJU_simes_test.", p, ".sigGenes_0.05.tsv"), 
#               sep="\t", quote=FALSE, row.names=FALSE)  
#   
#   write.csv(DEJU_F_0.05, paste0("DEJU_F_test.", p, ".sigGenes_0.05.csv"), row.names=FALSE)
#   write.table(DEJU_F_0.05, paste0("DEJU_F_test.", p, ".sigGenes_0.05.tsv"), 
#               sep="\t", quote=FALSE, row.names=FALSE)  
#   
#   write.csv(DEJU_exon_0.05, paste0("DEJU_exon_test.", p, ".sigGenes_0.05.csv"), row.names=FALSE)
#   write.table(DEJU_exon_0.05, paste0("DEJU_exon_test.", p, ".sigGenes_0.05.tsv"), 
#               sep="\t", quote=FALSE, row.names=FALSE)
#   
#   
# } else {
#   
#   print("Invalid mode specified. Choose 'simulation' or 'case_study'!")
# 
# }

#################################################################################
