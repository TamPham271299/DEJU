OK <- requireNamespace("devtools", quietly = TRUE)
if (!OK) {
  stop("devtools package required but is not installed (or can't be loaded)")
}
library(Rsubread)
library(edgeR)
library(limma)
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
  
  message("Transforming data for LM ...")
  v <- voom(y, design, plot=FALSE)
  
  message("Fitting LM for design matrix ...")
  fit <- lmFit(v, design)
  
  return(fit)

}

DEU_res <- function(fit, design, p, m, fdr_cutoff, mode) {
  cmd <- paste("makeContrasts(", p, ",levels=design)", sep="")
  contr <- eval(parse(text=cmd))
  print(contr)
  cfit <- contrasts.fit(fit, contr)
  
  message("Running diffSplice ..")
  sp <- diffSplice(cfit, geneid="GeneID", robust=TRUE, exonid="Start")
  DEU_simes <- topSplice(sp, test="simes", number=Inf)
  DEU_F <- topSplice(sp, test="F", number=Inf)
  DEU_exon <- topSplice(sp, test="t", number=Inf)
  
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
  
  ### 02. Run limma
  message('Start running limma ...')
  bplapply(simID, FUN = function(i){
    
    message(paste0('Processing: ', i, " ..."))
    
    message("Creating output folders to store read counts and DEU results ...")
    featureCounts_o <- paste0(DIR, "featureCounts/", i, "/")
    dir.create(featureCounts_o, showWarnings = FALSE, recursive = TRUE)
    DEU_analysis_o <- paste0(DIR, "limma_diffSplice/", i, "/")
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
    
    message("Starting DEU analysis for DEU-limma method ...")
    DEU_fit <- DEU_analysis(E$E_count,
                            E$E_annot,
                            mat$group, mat$design, mode)
    
    message("Combining internal exon and junction read count matrix ...")  
    IE_J <- process_IE_J_count(featureCounts_o, SJ_database)
    
    message("Starting DEJU analysis for DEJU-limma ...")
    DEJU_fit <- DEU_analysis(IE_J$IE_J_count,
                             IE_J$IE_J_annot,
                             mat$group, mat$design, mode)
    
    setwd(DEU_analysis_o)
    message("Running diffSplice for DEU-limma ...")
    res <- DEU_res(DEU_fit, mat$design, pair, "DEU", fdr_cutoff, mode)  
    
    message("Running diffSplice for DEJU-limma ...")
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
  
  message("Starting DEU analysis for DEU-limma method ...")
  DEU_fit <- DEU_analysis(E$E_count,
                          E$E_annot,
                          mat$group, mat$design, mode)
  
  message("Combining internal exon and junction read count matrix ...")  
  IE_J <- process_IE_J_count(featureCounts_o, SJ_database)
  
  message("Starting DEJU analysis for DEJU-limma ...")
  DEJU_fit <- DEU_analysis(IE_J$IE_J_count,
                           IE_J$IE_J_annot,
                           mat$group, mat$design, mode)
  
  pair <- strsplit(pair, ",")[[1]]
  for (p in pair) {
    
    DEU_analysis_o <- paste0(DIR, "limma_diffSplice/", p, "/")
    dir.create(DEU_analysis_o, recursive = TRUE, showWarnings = FALSE) 
    
    setwd(DEU_analysis_o)
    message("Running diffSplice for DEU-limma ...")
    res <- DEU_res(DEU_fit, mat$design, p, "DEU", fdr_cutoff, mode)  
    
    message("Running diffSplice for DEJU-limma ...")
    res <- DEU_res(DEJU_fit, mat$design, p, "DEJU", fdr_cutoff, mode)  
    setwd(wd)
  }
  
} else {
  
  print("Invalid mode specified. Choose 'simulation' or 'case_study'!")
  
}

#################################################################################
# OK <- requireNamespace("devtools", quietly = TRUE)
# if (!OK) {
#   stop("devtools package required but is not installed (or can't be loaded)")
# }
# library(Rsubread)
# library(limma)
# library(BiocParallel)
# 
# quantify_IE_count <- function(BAM_files, isPairedEnd, flat_exon, featureCounts_o) {
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
#   write.table(internal_exon_count$counts, 
#               "internal_exon_count.tsv", sep="\t", quote=FALSE, row.names=FALSE)
#   write.table(internal_exon_count$annotation, 
#               "internal_exon_count_annotation.tsv", sep="\t", quote=FALSE, row.names=FALSE)
#   write.table(internal_exon_count$stat, 
#               "internal_exon_count_stat.tsv", sep="\t", quote=FALSE, row.names=FALSE)  
# }
# 
# quantify_E_J_count <- function(BAM_files, isPairedEnd, flat_exon, FASTA, featureCounts_o) {
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
#   write.table(junction_count$counts, 
#               "exon_count.tsv", sep="\t", quote=FALSE, row.names=FALSE)
#   write.table(junction_count$annotation, 
#               "exon_count_annotation.tsv", sep="\t", quote=FALSE, row.names=FALSE)
#   write.table(junction_count$stat, 
#               "exon_count_stat.tsv", sep="\t", quote=FALSE, row.names=FALSE)
#   
#   message("Saving junction read counts and annotations ...")  
#   write.table(junction_count$counts_junction[!is.na(junction_count$counts_junction$PrimaryGene),], 
#               "junction_count.tsv", sep="\t", quote=FALSE, row.names=FALSE)
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
#   g <- unlist(strsplit(p, "-"))
#   contr <- limma::makeContrasts(contrast = paste(g[1], "-", g[2]), levels = design)
#   
#   return(list('group'=group, 'samples'=samples, 'design'=design, 'contr'=contr))
# }
# 
# process_E_count_withDupJ <- function() {
#   
#   message("Constructing a matrix for exon read counts with duplicated junction count ...")
#   E_count <- read.table("exon_count.tsv", header=TRUE)
#   E_count <- data.frame(E_count[,1:ncol(E_count)], row.names=NULL, check.names = FALSE)
#   E_annot <- read.table("exon_count_annotation.tsv", header=TRUE)
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
# process_IE_J_count <- function(SJ_database) {
#   
#   message("Analyzing internal exon read counts ...")
#   IE_count <- read.table("internal_exon_count.tsv", header=TRUE)
#   IE_count <- data.frame(IE_count[,1:ncol(IE_count)], row.names=NULL, check.names = FALSE)
#   IE_annot <- read.table("internal_exon_count_annotation.tsv", header=TRUE)
#   IE_annot <- cbind(IE_annot, Region="Exon", annotated=1)
#   
#   message("Analyzing junction read counts ...")
#   J_count <- read.table("junction_count.tsv", header=TRUE)
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
#   v <- voom(y, design, plot=FALSE)
#   
#   message("Fitting GLM model for design matrix ...")
#   fit <- lmFit(v, design)
#   
#   message("Running diffSplice ..")
#   sp <- diffSplice(fit, geneid="GeneID", robust=TRUE, exonid="Start")
#   DEU_simes <- topSplice(sp, test="simes", number=Inf)
#   DEU_F <- topSplice(sp, test="F", number=Inf)
#   DEU_exon <- topSplice(sp, test="t", number=Inf)
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
#   message("Saving all DEU results for DEU-limma")
#   write.table(DEU$simes, paste0("DEU_simes_test.allGenes.tsv"), 
#               sep="\t", quote=FALSE, row.names=FALSE)
#   write.table(DEU$F, paste0("DEU_F_test.allGenes.tsv"), 
#               sep="\t", quote=FALSE, row.names=FALSE)
#   write.table(DEU$exon, paste0("DEU_exon_test.allGenes.tsv"), 
#               sep="\t", quote=FALSE, row.names=FALSE)
#   
#   message("Saving significant DEU results for DEU-limma ...")
#   write.table(DEU$sig_simes, paste0("DEU_simes_test.sigGenes_", fdr_cutoff, ".tsv"), 
#               sep="\t", quote=FALSE, row.names=FALSE)  
#   write.table(DEU$sig_F, paste0("DEU_F_test.sigGenes_", fdr_cutoff, ".tsv"), 
#               sep="\t", quote=FALSE, row.names=FALSE)  
#   write.table(DEU$sig_exon, paste0("DEU_exon_test.sigGenes_", fdr_cutoff, ".tsv"), 
#               sep="\t", quote=FALSE, row.names=FALSE)  
#   
#   message("Saving all DEU results for DEJU-limma ...")
#   write.table(DEJU$simes, paste0("DEJU_simes_test.allGenes.tsv"), 
#               sep="\t", quote=FALSE, row.names=FALSE) 
#   write.table(DEJU$F, paste0("DEJU_F_test.allGenes.tsv"), 
#               sep="\t", quote=FALSE, row.names=FALSE) 
#   write.table(DEJU$exon, paste0("DEJU_exon_test.allGenes.tsv"), 
#               sep="\t", quote=FALSE, row.names=FALSE) 
#   
#   message("Saving significant DEU results for DEJU-limma...")  
#   write.table(DEJU$sig_simes, paste0("DEJU_simes_test.sigGenes_", fdr_cutoff, ".tsv"), 
#               sep="\t", quote=FALSE, row.names=FALSE)  
#   write.table(DEJU$sig_F, paste0("DEJU_F_test.sigGenes_", fdr_cutoff, ".tsv"), 
#               sep="\t", quote=FALSE, row.names=FALSE)  
#   write.table(DEJU$sig_exon, paste0("DEJU_exon_test.sigGenes_", fdr_cutoff, ".tsv"), 
#               sep="\t", quote=FALSE, row.names=FALSE)
#   
# }
# 
# save_csv <- function(DEU, DEJU, fdr_cutoff) {
#   message("Saving all DEU results for DEU-limma ...")
#   write.csv(DEU$simes, "DEU_simes_test.allGenes.csv", row.names=FALSE)
#   write.csv(DEU$F, "DEU_F_test.allGenes.csv", row.names=FALSE)
#   write.csv(DEU$exon, "DEU_exon_test.allGenes.csv", row.names=FALSE)
#   
#   message("Saving significant DEU results for DEU-limma ...")
#   write.csv(DEU$sig_simes, paste0("DEU_simes_test.sigGenes_", fdr_cutoff, ".csv"), 
#             row.names=FALSE)
#   write.csv(DEU$sig_F, paste0("DEU_F_test.sigGenes_", fdr_cutoff, ".csv"), 
#             row.names=FALSE)
#   write.csv(DEU$sig_exon, paste0("DEU_exon_test.sigGenes_", fdr_cutoff, ".csv"), 
#             row.names=FALSE)
#   
#   message("Saving all DEU results for DEJU-limma ...")  
#   write.csv(DEJU$simes, "DEJU_simes_test.allGenes.csv", 
#             row.names=FALSE)
#   write.csv(DEJU$F, "DEJU_F_test.allGenes.csv", 
#             row.names=FALSE)
#   write.csv(DEJU$exon, "DEJU_exon_test.allGenes.csv", 
#             row.names=FALSE)
#   
#   message("Saving significant DEU results for DEJU-limma  ...")  
#   write.csv(DEJU$sig_simes, paste0("DEJU_simes_test.sigGenes_", fdr_cutoff, ".csv"), 
#             row.names=FALSE)
#   write.csv(DEJU$sig_F, paste0("DEJU_F_test.sigGenes_", fdr_cutoff, ".csv"), 
#             row.names=FALSE)
#   write.csv(DEJU$sig_exon, paste0("DEJU_exon_test.sigGenes_", fdr_cutoff, ".csv"), 
#             row.names=FALSE)
# }
# 
# # Getting parameters
# args <- R.utils::commandArgs(asValues = TRUE)
# print(args)
# 
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
# f <- c("internal_exon_count.tsv", "internal_exon_count_annotation.tsv", "internal_exon_count_stat.tsv",
#        "exon_count.tsv", "exon_count_annotation.tsv", "exon_count_stat.tsv",
#        "junction_count.tsv")
# 
# message("Check parameters ...")
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
#   ### 02. Run limma
#   message('Start running limma ...')
#   bplapply(simID, FUN = function(i){
#     
#     message(paste0('Processing: ', i, " ..."))
#     
#     message("Creating output folders to store read counts and DEU results ...")
#     featureCounts_o <- paste0(DIR, "featureCounts/", i, "/")
#     dir.create(featureCounts_o, showWarnings = FALSE, recursive = TRUE)
#     DEU_analysis_o <- paste0(DIR, "limma_diffSplice/", i, "/")
#     dir.create(DEU_analysis_o, showWarnings = FALSE, recursive = TRUE)
#     
#     message("Constructing design matrix for contrast ...")
#     mat <- construct_model_matrix(targetp, p)
#     
#     BAM_files <- paste0(DIR, 
#                         "aligned_pass2/", 
#                         i, "/", 
#                         mat$samples, "/", 
#                         mat$samples, 
#                         ".Aligned.sortedByCoord.out.bam")
#     
#     setwd(featureCounts_o)
#     
#     message("Checking if running featureCounts yet ...")
#     if (any(!file.exists(f))) {
#       
#       message("Quantifing internal exon read count ...")
#       IE_count_table <- quantify_IE_count(BAM_files, isPairedEnd, flat_exon, featureCounts_o)
#       
#       message("Quantifying exon-level read with duplicated junction count and junction read count ...")
#       E_J_count_table <- quantify_E_J_count(BAM_files, isPairedEnd, flat_exon, FASTA, featureCounts_o)
#       
#     }
#     
#     message("Constructing a matrix for exon read counts with duplicated junction count ...")
#     E <- process_E_count_withDupJ()
#     
#     message("Starting DEU analysis for DEU-limma method ...")
#     DEU <- DEU_analysis(E$E_count, 
#                         E$E_annot, 
#                         mat$group, mat$design, mat$contr, mode, fdr_cutoff)
#     
#     message("Combining internal exon and junction read count matrix ...")  
#     IE_J <- process_IE_J_count(SJ_database)
#     
#     message("Starting DEJU analysis for DEJU-limma ...")
#     DEJU <- DEU_analysis(IE_J$IE_J_count, 
#                          IE_J$IE_J_annot, 
#                          mat$group, mat$design, mat$contr, mode, fdr_cutoff)
#     
#     setwd(DEU_analysis_o)
#     message("Saving final DEU results ...")
#     res <- save_results(DEU, DEJU, fdr_cutoff)
#     
#   }, BPPARAM = BPPARAM)
#   
# } else if (mode == "case_study") {
#   
#   message("Creating output folders to store read counts and DEU results ...")
#   featureCounts_o <- paste0(DIR, "featureCounts/")
#   dir.create(featureCounts_o, recursive = TRUE, showWarnings = FALSE)
#   DEU_analysis_o <- paste0(DIR, "limma_diffSplice/", p, "/")
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
#   setwd(featureCounts_o)
#   
#   message("Checking if running featureCounts yet ...")
#   if (any(!file.exists(f))) {
#     
#     message("Quantifing internal exon read count ...")
#     IE_count_table <- quantify_IE_count(BAM_files, isPairedEnd, flat_exon, featureCounts_o)
#     
#     message("Quantifying exon-level read with duplicated junction count and junction read count ...")
#     E_J_count_table <- quantify_E_J_count(BAM_files, isPairedEnd, flat_exon, FASTA, featureCounts_o)
#     
#   }
#   
#   message("Constructing a matrix for exon read counts with duplicated junction count ...")
#   E <- process_E_count_withDupJ()
#   
#   message("Starting DEU analysis for DEU-limma method ...")
#   DEU <- DEU_analysis(E$E_count, 
#                       E$E_annot, 
#                       mat$group, mat$design, mat$contr, mode, fdr_cutoff)
#   
#   message("Combining internal exon and junction read count matrix ...")  
#   IE_J <- process_IE_J_count()
#   
#   message("Starting DEJU analysis for DEJU-limma ...")
#   DEJU <- DEU_analysis(IE_J$IE_J_count, 
#                        IE_J$IE_J_annot, 
#                        mat$group, mat$design, mat$contr, mode, fdr_cutoff)
#   
#   setwd(DEU_analysis_o)
#   message("Saving final DEU results ...")
#   res <- save_results(DEU, DEJU, fdr_cutoff)
#   
#   message("Saving final DEU results in .csv ...")
#   res_csv <- save_csv(DEU, DEJU, fdr_cutoff)
#   
# } else {
#   
#   print("Invalid mode specified. Choose 'simulation' or 'case_study'!")
#   
# }

#################################################################################
# OK <- requireNamespace("devtools", quietly = TRUE)
# if (!OK) {
#   stop("devtools package required but is not installed (or can't be loaded)")
# }
# library(limma)
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
# # DIR <- "/vast/projects/Spatial/tam/Differential_splicing/simulation_3vs3_mouse_genome/RNA-seq/DEU_mix/unbalanced_3_75_3/"
# # REF <- "/stornext/General/data/academic/lab_chen/tam/ref_genome/Mus_musculus/Gencode/"
# # targetp <- "/stornext/General/data/academic/lab_chen/tam/Differential_splicing/simulation_3vs3_mouse_genome/RNA-seq/config/design.tsv"
# 
# path_to_sj <- paste0(REF, "gencode.vM32.primary_assembly.annotation.uniqSJInfo.tsv")
# num_of_occur_sj <- paste0(REF, "gencode.vM32.primary_assembly.annotation.uniqSJInfo.numOfOccur.tsv")
# 
# annotated_sj <- read.table(path_to_sj, header=TRUE)
# numOfOccur_sj <- read.table(num_of_occur_sj, header=TRUE)
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
#   ### 02. Run limma
#   message('Start running limma ...')
#   bplapply(simID, FUN = function(i){
#     
#     message(paste0('Processing: ', i, " ..."))
#     OUTPUT <- paste0(DIR, "STAR_output/limma_fit_RL_75/", i, "/")
#     readCount_output <- paste0(DIR, "STAR_output/featureCounts_fit_RL_75/minMQS_255/", i, "/")
#     dir.create(OUTPUT, showWarnings = FALSE, recursive = TRUE)
#     
#     ### Read count quantification
#     message('Loading read counts at exon, junction ...')
#     # Exon-level read counts with double-counted junctions
#     exon_count_f <- read.table(paste0(readCount_output, "exon_count.tsv"), header=TRUE)
#     exon_count <- data.frame(exon_count_f[,1:ncol(exon_count_f)], row.names=NULL, check.names = FALSE)
#     exon_annot <- read.table(paste0(readCount_output, "exon_count_annotation.tsv"), header=TRUE)
#     
#     # Internal exon read counts
#     internal_exon_count_f <- read.table(paste0(readCount_output, "internal_exon_count.tsv"), header=TRUE)
#     internal_exon_count <- data.frame(internal_exon_count_f[,1:ncol(internal_exon_count_f)], row.names=NULL, check.names = FALSE)
#     internal_exon_annot <- read.table(paste0(readCount_output, "internal_exon_count_annotation.tsv"), header=TRUE)
#     
#     # Junction read counts
#     junc_count_f <- read.table(paste0(readCount_output, "junction_count.tsv"), header=TRUE)
#     junc_count_f$juncID <- paste(junc_count_f$Site1_chr, junc_count_f$Site1_location, junc_count_f$Site2_location, sep="_")
#     m_numOfOccur <- match(annotated_sj$juncID, numOfOccur_sj$juncID)
#     annotated_sj$numOfOccur <- numOfOccur_sj$numOfOccur[m_numOfOccur]
#     annotated_sj_numOfOccur_1 <- annotated_sj[annotated_sj$numOfOccur==1,]
#     m1 <- match(junc_count_f$juncID, annotated_sj_numOfOccur_1$juncID)
#     junc_count_f$PrimaryGene <- ifelse(!is.na(m1), annotated_sj_numOfOccur_1$geneID[m1], junc_count_f$PrimaryGene)
#     m2 <- match(junc_count_f$juncID, annotated_sj$juncID)
#     junc_count_f$annotated <- ifelse(!is.na(m2), 1, 0)
#     junc_count <- data.frame(junc_count_f[, 9:(ncol(junc_count_f)-2)], row.names = NULL, check.names = FALSE)
#     junc_annot <- data.frame(
#       GeneID=junc_count_f$PrimaryGene,
#       Chr=junc_count_f$Site1_chr,
#       Start=junc_count_f$Site1_location,
#       End=junc_count_f$Site2_location
#     )
#     m <- match(junc_annot$GeneID, internal_exon_annot$GeneID)
#     Strand <- internal_exon_annot$Strand[m]
#     junc_annot <- cbind(junc_annot, Strand=Strand, Length=1, Region="Junction", annotated=junc_count_f$annotated)
#     
#     internal_exon_junc_count <- rbind(internal_exon_count, junc_count)
#     internal_exon_annot <- cbind(internal_exon_annot, Region="Exon", annotated=1)
#     internal_exon_junc_annot <- rbind(internal_exon_annot, junc_annot)
#     
#     ### Preliminary analysis
#     message('Preliminary analysis ...')
#     ### design model
#     target <- read.table(targetp, sep="\t", header=TRUE, stringsAsFactors=FALSE)
#     group <- factor(target$Group, levels=c("Group_1","Group_2"))
#     design <- model.matrix(~ group)
#     
#     ## Existing approach
#     y_current_ap <- DGEList(counts=exon_count, genes=exon_annot, group=group)
#     colnames(y_current_ap) <- gsub("[.].*$", "", colnames(y_current_ap))
#     keep <- filterByExpr(y_current_ap, group=group)
#     y_current_ap <- y_current_ap[keep, , keep.lib.sizes=FALSE]
#     y_current_ap <- normLibSizes(y_current_ap)
#     v_current_ap <- voom(y_current_ap, design, plot=FALSE)
#     print(y_current_ap)
#     v_current_ap_fit <- lmFit(v_current_ap, design)
#     
#     ## Junction approach (internal exon read + junction read)
#     y_new_ap <- DGEList(counts=internal_exon_junc_count, genes=internal_exon_junc_annot, group=group)
#     colnames(y_new_ap) <- gsub("[.].*$", "", colnames(y_new_ap))
#     keep <- filterByExpr(y_new_ap, group=group)
#     y_new_ap <- y_new_ap[keep, , keep.lib.sizes=FALSE]
#     y_new_ap <- normLibSizes(y_new_ap)
#     v_new_ap <- voom(y_new_ap, design, plot=FALSE)
#     print(y_new_ap)
#     v_new_ap_fit <- lmFit(v_new_ap, design)
#     
#     ### Differential alternative splicing analysis
#     message('DEU analysis ...')
#     
#     message('Save diffSplice results ...')
#     ## Existing approach
#     ds_v_current_ap <- diffSplice(v_current_ap_fit, geneid="GeneID", robust=TRUE, exonid="Start")
#     v_current_ap_simes <- topSplice(ds_v_current_ap, test="simes", n=Inf)
#     v_current_ap_gene <- topSplice(ds_v_current_ap, test="F", n=Inf)
#     v_current_ap_exon <- topSplice(ds_v_current_ap, test="t", n=Inf)
#     write.table(v_current_ap_simes, file.path(OUTPUT, "DEU_current_approach_simes_minMQS_255.allGenes.tsv"), row.names=FALSE, sep="\t", quote=FALSE)
#     write.table(v_current_ap_gene, file.path(OUTPUT, "DEU_current_approach_gene_minMQS_255.allGenes.tsv"), row.names=FALSE, sep="\t", quote=FALSE)
#     write.table(v_current_ap_exon, file.path(OUTPUT, "DEU_current_approach_exon_minMQS_255.allGenes.tsv"), row.names=FALSE, sep="\t", quote=FALSE)
#     
#     ## Junction approach
#     ds_v_new_ap <- diffSplice(v_new_ap_fit, geneid="GeneID", robust=TRUE, exonid="Start")
#     v_new_ap_simes <- topSplice(ds_v_new_ap, test="simes", n=Inf)
#     v_new_ap_gene <- topSplice(ds_v_new_ap, test="F", n=Inf)
#     v_new_ap_exon <- topSplice(ds_v_new_ap, test="t", n=Inf)
#     write.table(v_new_ap_simes, file.path(OUTPUT, "DEU_new_approach_simes_minMQS_255_SJ.allGenes.tsv"), row.names=FALSE, sep="\t", quote=FALSE)
#     write.table(v_new_ap_gene, file.path(OUTPUT, "DEU_new_approach_gene_minMQS_255_SJ.allGenes.tsv"), row.names=FALSE, sep="\t", quote=FALSE)
#     write.table(v_new_ap_exon, file.path(OUTPUT, "DEU_new_approach_exon_minMQS_255_SJ.allGenes.tsv"), row.names=FALSE, sep="\t", quote=FALSE)
#     
#     message('Save significant DEUs ...')
#     ## Existing approach
#     v_current_ap_simes <- topSplice(ds_v_current_ap, test="simes", n=NULL, FDR=0.05)
#     v_current_ap_gene <- topSplice(ds_v_current_ap, test="F", n=NULL, FDR=0.05)
#     write.table(v_current_ap_simes, file.path(OUTPUT, "DEU_current_approach_simes_minMQS_255.tsv"), row.names=FALSE, sep="\t", quote=FALSE)
#     write.table(v_current_ap_gene, file.path(OUTPUT, "DEU_current_approach_gene_minMQS_255.tsv"), row.names=FALSE, sep="\t", quote=FALSE)
#     
#     ## Junction approach
#     v_new_ap_simes <- topSplice(ds_v_new_ap, test="simes", n=NULL, FDR=0.05)
#     v_new_ap_gene <- topSplice(ds_v_new_ap, test="F", n=NULL, FDR=0.05)
#     write.table(v_new_ap_simes, file.path(OUTPUT, "DEU_new_approach_simes_minMQS_255_SJ.tsv"), row.names=FALSE, sep="\t", quote=FALSE)
#     write.table(v_new_ap_gene, file.path(OUTPUT, "DEU_new_approach_gene_minMQS_255_SJ.tsv"), row.names=FALSE, sep="\t", quote=FALSE)
#     
#   }, BPPARAM = BPPARAM)
#   
# } else if (mode == "case_study") {
#   
#   target <- read.table(targetp, sep="\t", header=TRUE, stringsAsFactors=FALSE)
#   group <- factor(target$Group)
#   
#   OUTPUT <- paste0(DIR, "STAR_output/limma/")
#   dir.create(OUTPUT, recursive = TRUE, showWarnings = F)
#   readCount_output <- paste0(DIR, "STAR_output/featureCounts/")
#   
#   gene_info_p <- paste0(REF, "gene_info_GRCm39_M32_ensembl109.tsv")
#   gene_info <- read.table(gene_info_p, header=T, na.strings="." )
#   
#   exon_count_f <- read.table(paste0(readCount_output, "exon_count.tsv"), header=TRUE)
#   exon_count <- data.frame(exon_count_f[,1:ncol(exon_count_f)], row.names=NULL, check.names = FALSE)
#   colnames(exon_count) <- gsub("[.].*$", "", colnames(exon_count))
#   exon_annot <- read.table(paste0(readCount_output, "exon_count_annotation.tsv"), header=TRUE)
#   
#   ## Junction approach
#   # internal exon counts and annotation
#   internal_exon_count_f <- read.table(paste0(readCount_output, "internal_exon_count.tsv"), header=TRUE)
#   internal_exon_count <- data.frame(internal_exon_count_f[,1:ncol(internal_exon_count_f)], row.names=NULL, check.names = FALSE)
#   colnames(internal_exon_count) <- gsub("[.].*$", "", colnames(internal_exon_count))
#   internal_exon_annot <- read.table(paste0(readCount_output, "internal_exon_count_annotation.tsv"), header=TRUE)
#   
#   # Junction counts and its annotation
#   junc_count_f <- read.table(paste0(readCount_output, "junction_count.tsv"), header=TRUE)
#   junc_count_f$juncID <- paste(junc_count_f$Site1_chr, junc_count_f$Site1_location, junc_count_f$Site2_location, sep="_")
#   
#   # match juncID with number of occurences
#   m_numOfOccur <- match(annotated_sj$juncID, numOfOccur_sj$juncID)
#   annotated_sj$numOfOccur <- numOfOccur_sj$numOfOccur[m_numOfOccur]
#   
#   # only consider annotated sj that has numOfOccur equals to 1
#   annotated_sj_numOfOccur_1 <- annotated_sj[annotated_sj$numOfOccur==1,]
#   m1 <- match(junc_count_f$juncID, annotated_sj_numOfOccur_1$juncID)
#   junc_count_f$PrimaryGene <- ifelse(!is.na(m1), annotated_sj_numOfOccur_1$geneID[m1], junc_count_f$PrimaryGene)
#   m2 <- match(junc_count_f$juncID, annotated_sj$juncID)
#   junc_count_f$annotated <- ifelse(!is.na(m2), 1, 0)
#   junc_count <- data.frame(junc_count_f[, 9:(ncol(junc_count_f)-2)], row.names = NULL, check.names = FALSE)
#   colnames(junc_count) <- gsub("[.].*$", "", colnames(junc_count))
#   junc_annot <- data.frame(
#     GeneID=junc_count_f$PrimaryGene,
#     Chr=junc_count_f$Site1_chr,
#     Start=junc_count_f$Site1_location,
#     End=junc_count_f$Site2_location
#   )
#   m <- match(junc_annot$GeneID, internal_exon_annot$GeneID)
#   Strand <- internal_exon_annot$Strand[m]
#   junc_annot <- cbind(junc_annot, Strand=Strand, Length=1, Region="Junction", annotated=junc_count_f$annotated)
#   
#   # When we have internal exon count table and junction count table, we combine them into a internal exon and junction count table, which will be used to analyze differential splicing for the junction approach.
#   internal_exon_junc_count <- rbind(internal_exon_count, junc_count)
#   internal_exon_annot <- cbind(internal_exon_annot, Region="Exon", annotated=1)
#   internal_exon_junc_annot <- rbind(internal_exon_annot, junc_annot)
#   
#   ### Preliminary analysis
#   
#   # We firstly set up a design matrix for the differential analysis.
#   design <- model.matrix(~ 0 + group)
#   colnames(design) <- gsub("group", "", colnames(design))
#   
#   contr <- makeContrasts(
#     BvsL = Basal - LP,
#     BvsM = Basal - ML,
#     LvsM = LP - ML, 
#     levels=design)
#   
#   ## Existing approach
#   # We create DGEList object to store general exon counts and its annotation.
#   y_existing_ap <- DGEList(counts=exon_count, genes=exon_annot, group=group)
#   
#   # We add gene symbol to the exon annotation. Any missing gene symbol which marked as "." will be removed.
#   m <- match(y_existing_ap$genes$GeneID, gene_info$GeneID)
#   y_existing_ap$genes$Symbol <- gene_info$Symbol[m]
#   keep <- !is.na(y_existing_ap$genes$Symbol)
#   y_existing_ap <- y_existing_ap[keep, , keep.lib.sizes=FALSE]
#   
#   # Filtering
#   keep <- filterByExpr(y_existing_ap, group=group)
#   y_existing_ap <- y_existing_ap[keep, , keep.lib.sizes=FALSE]
#   
#   # Normalization
#   y_existing_ap <- normLibSizes(y_existing_ap)
#   
#   # Fitting model using voom
#   v_existing_ap <- voom(y_existing_ap, design, plot=TRUE)
#   v_existing_ap_fit <- lmFit(v_existing_ap, design)
#   v_existing_ap_cfit <- contrasts.fit(v_existing_ap_fit, contr)
#   
#   ## Junction Approach
#   # Like the existing approach, we also create the DEGList object to store the internal exonic and junction read counts and annotation.
#   y_junction_ap <- DGEList(counts=internal_exon_junc_count, genes=internal_exon_junc_annot, group=group)
#   
#   # We also add the gene symbol to the internal exonic and junction annotation.
#   m <- match(y_junction_ap$genes$GeneID, gene_info$GeneID)
#   y_junction_ap$genes$Symbol <- gene_info$Symbol[m]
#   keep <- !is.na(y_junction_ap$genes$Symbol)
#   y_junction_ap <- y_junction_ap[keep, , keep.lib.sizes=FALSE]
#   
#   # Filtering by expression 
#   keep <- filterByExpr(y_junction_ap, group=group)
#   y_junction_ap <- y_junction_ap[keep, , keep.lib.sizes=FALSE]
#   
#   # Normalizing lib sizes
#   y_junction_ap <- normLibSizes(y_junction_ap)
#   
#   # Fitting linear model using voom
#   v_junction_ap <- voom(y_junction_ap, design, plot=TRUE)
#   v_junction_ap_fit <- lmFit(v_junction_ap, design)
#   v_junction_ap_cfit <- contrasts.fit(v_junction_ap_fit, contr)
#   
#   ### DEJU analysis
#   gp_type <- c("Basal_LP", "Basal_ML", "LP_ML")
#   
#   ## Existing approach
#   sp_existing_ap <- diffSplice(v_existing_ap_cfit, geneid="GeneID", robust=TRUE, exonid="Start")
#   existing_ap_simes <- list()
#   existing_ap_gene <- list()
#   existing_ap_exon <- list()
#   coef_n <- 1
#   
#   for (gp in gp_type) {
#     existing_ap_simes[[gp]] <- topSplice(sp_existing_ap, coef=coef_n, test="simes", number=Inf)
#     existing_ap_gene[[gp]] <- topSplice(sp_existing_ap, coef=coef_n, test="F", number=Inf)
#     existing_ap_exon[[gp]] <- topSplice(sp_existing_ap, coef=coef_n, test="t", number=Inf)
#     dir.create(paste0(OUTPUT, gp), recursive = TRUE, showWarnings = F)
#     write.csv(existing_ap_simes[[gp]],
#               file.path(paste0(OUTPUT, gp), paste0("DEU_existing_approach_simes_test_", gp, ".allGenes.csv")),
#               row.names=FALSE)
#     write.csv(existing_ap_gene[[gp]],
#               file.path(paste0(OUTPUT, gp), paste0("DEU_existing_approach_gene_test_", gp, ".allGenes.csv")),
#               row.names=FALSE)
#     write.csv(existing_ap_exon[[gp]],
#               file.path(paste0(OUTPUT, gp), paste0("DEU_existing_approach_exon_test_", gp, ".allGenes.csv")),
#               row.names=FALSE)
#     coef_n <- coef_n + 1
#   }
#   
#   ## Junction approach
#   sp_junction_ap <- diffSplice(v_junction_ap_cfit, geneid="GeneID", robust=TRUE, exonid="Start")
#   junction_ap_simes <- list()
#   junction_ap_gene <- list()
#   junction_ap_exon <- list()
#   coef_n <- 1
#   
#   for (gp in gp_type) {
#     junction_ap_simes[[gp]] <- topSplice(sp_junction_ap, coef=coef_n, test="simes", number=Inf)
#     junction_ap_gene[[gp]] <- topSplice(sp_junction_ap, coef=coef_n, test="F", number=Inf)
#     junction_ap_exon[[gp]] <- topSplice(sp_junction_ap, coef=coef_n, test="t", number=Inf)
#     write.csv(junction_ap_simes[[gp]],
#               file.path(paste0(OUTPUT, gp), paste0("DEU_junction_approach_simes_test_", gp, "_with_SJ.allGenes.csv")),
#               row.names=FALSE)
#     write.csv(junction_ap_gene[[gp]],
#               file.path(paste0(OUTPUT, gp), paste0("DEU_junction_approach_gene_test_", gp, "_with_SJ.allGenes.csv")),
#               row.names=FALSE)
#     write.csv(junction_ap_exon[[gp]],
#               file.path(paste0(OUTPUT, gp), paste0("DEU_junction_approach_exon_test_", gp, "_with_SJ.allGenes.csv")),
#               row.names=FALSE)
#     coef_n <- coef_n + 1
#   }
#   
#   ### Significant DEUs
#   # Setting FDR cutoff 0.05, we extract the significant DEU detected by both approaches and save significant DEU genes.
#   
#   existing_ap_simes_0.05 <- list()
#   existing_ap_gene_0.05 <- list()
#   junction_ap_simes_0.05 <- list()
#   junction_ap_gene_0.05 <- list()
#   
#   for (gp in gp_type) {
#     
#     ### Existing method
#     existing_ap_simes_0.05[[gp]] <- existing_ap_simes[[gp]][existing_ap_simes[[gp]]$FDR<=0.05,]
#     write.table(existing_ap_simes_0.05[[gp]],
#                 file.path(paste0(OUTPUT, gp), paste0("DEU_existing_approach_simes_test_", gp, ".sigGenes.tsv")),
#                 sep="\t", quote=FALSE, row.names=FALSE)
#     write.csv(existing_ap_simes_0.05[[gp]],
#               file.path(paste0(OUTPUT, gp), paste0("DEU_existing_approach_simes_test_", gp, ".sigGenes.csv")),
#               row.names=FALSE)
#     
#     existing_ap_gene_0.05[[gp]] <- existing_ap_gene[[gp]][existing_ap_gene[[gp]]$FDR<=0.05,]
#     write.table(existing_ap_gene_0.05[[gp]],
#                 file.path(paste0(OUTPUT, gp), paste0("DEU_existing_approach_gene_test_", gp, ".sigGenes.tsv")),
#                 sep="\t", quote=FALSE, row.names=FALSE)
#     write.csv(existing_ap_gene_0.05[[gp]],
#               file.path(paste0(OUTPUT, gp), paste0("DEU_existing_approach_gene_test_", gp, ".sigGenes.csv")),
#               row.names=FALSE)
#     
#     ### Junction method
#     junction_ap_simes_0.05[[gp]] <- junction_ap_simes[[gp]][junction_ap_simes[[gp]]$FDR<=0.05,]
#     write.table(junction_ap_simes_0.05[[gp]], 
#                 file.path(paste0(OUTPUT, gp), paste0("DEU_junction_approach_simes_test_", gp, "_with_SJ.sigGenes.tsv")), 
#                 sep="\t", quote=FALSE, row.names=FALSE)
#     write.csv(junction_ap_simes_0.05[[gp]], 
#               file.path(paste0(OUTPUT, gp), paste0("DEU_junction_approach_simes_test_", gp, "_with_SJ.sigGenes.csv")), 
#               row.names=FALSE)
#     
#     junction_ap_gene_0.05[[gp]] <- junction_ap_gene[[gp]][junction_ap_gene[[gp]]$FDR<=0.05,]
#     write.table(junction_ap_gene_0.05[[gp]], 
#                 file.path(paste0(OUTPUT, gp), paste0("DEU_junction_approach_gene_test_", gp, "_with_SJ.sigGenes.tsv")), 
#                 sep="\t", quote=FALSE, row.names=FALSE)
#     write.csv(junction_ap_gene_0.05[[gp]], 
#               file.path(paste0(OUTPUT, gp), paste0("DEU_junction_approach_gene_test_", gp, "_with_SJ.sigGenes.csv")), 
#               row.names=FALSE)
#     
#   }
#   
# } else {
# 
#   print("Invalid mode specified. Choose 'simulation' or 'case_study'!")
#   
# }
# 
