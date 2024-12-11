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

DEU_analysis <- function(r_count, annot, group, design, mode, REF) {
  
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
                            mat$group, mat$design, mode, REF)

    message("Combining internal exon and junction read count matrix ...")  
    IE_J <- process_IE_J_count(featureCounts_o, SJ_database)
    
    message("Starting DEJU analysis for DEJU-edgeR ...")
    DEJU_fit <- DEU_analysis(IE_J$IE_J_count,
                             IE_J$IE_J_annot,
                             mat$group, mat$design, mode, REF)

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
                          mat$group, mat$design, mode, REF)
 
  message("Combining internal exon and junction read count matrix ...")  
  IE_J <- process_IE_J_count(featureCounts_o, SJ_database)

  message("Start DEJU analysis for DEU-edgeR ...")
  DEJU_fit <- DEU_analysis(IE_J$IE_J_count,
                           IE_J$IE_J_annot,
                           mat$group, mat$design, mode, REF)
  
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
