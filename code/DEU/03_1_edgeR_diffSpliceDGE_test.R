OK <- requireNamespace("devtools", quietly = TRUE)
if (!OK) {
  stop("devtools package required but is not installed (or can't be loaded)")
}
library(Rsubread)
library(edgeR)
library(BiocParallel)

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

DEU_analysis <- function(r_count, annot, group, design) {
  
  message("Constructing DGElist object ...")  
  y <- DGEList(counts=r_count, genes=annot, group=group)
  colnames(y) <- gsub("[.].*$", "", colnames(y))
  
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

DEU_res <- function(fit, design, p) {
  
  message("Running diffSpliceDGE ..")
  cmd <- paste("makeContrasts(", p, ",levels=design)", sep="")
  contr <- eval(parse(text=cmd))
  print(contr)
  
  sp <- diffSpliceDGE(fit, contrast=contr, geneid="GeneID", exonid="Start")
  DEU_simes <- topSpliceDGE(sp, test="Simes", number=Inf)
  DEU_F <- topSpliceDGE(sp, test="gene", number=Inf)
  DEU_exon <- topSpliceDGE(sp, test="exon", number=Inf)
  
}

# Getting parameters
args <- R.utils::commandArgs(asValues = TRUE)
print(args)

DIR <- as.character(args[['DIR']])
REF <- as.character(args[['REF']])
targetp <- as.character(args[['target']])
pair <- as.character(args[['pair']])
i <- "S1"

FASTA <- paste0(REF, "GRCm39.primary_assembly.genome.fa.gz")
flat_exon <- paste0(REF, "gencode.vM32.annotation.flattened.exon.saf")
SJ <- paste0(REF, "gencode.vM32.primary_assembly.annotation.SJdatabase.tsv")
SJ_database <- read.table(SJ, header=TRUE)

message("Check parameters ...")
print(DIR)
print(REF)
print(targetp)
print(FASTA)
print(flat_exon)
print(SJ)
print(pair)
    
message(paste0('Processing: ', i, " ..."))

message("Creating output folders to store read counts and DEU results ...")
featureCounts_o <- paste0(DIR, "featureCounts/", i, "/")
dir.create(featureCounts_o, showWarnings = FALSE, recursive = TRUE)
DEU_analysis_o <- paste0(DIR, "edgeR_diffSpliceDGE/", i, "/")
dir.create(DEU_analysis_o, showWarnings = FALSE, recursive = TRUE)

message("Constructing design matrix for contrast ...")
mat <- construct_model_matrix(targetp)

message("Constructing a matrix for exon read counts with duplicated junction count ...")
E <- process_E_count_withDupJ(featureCounts_o)

message("Starting DEU analysis for DEU-edgeR ...")
DEU_fit <- DEU_analysis(E$E_count,
                        E$E_annot,
                        mat$group, mat$design)

message("Combining internal exon and junction read count matrix ...")  
IE_J <- process_IE_J_count(featureCounts_o, SJ_database)

message("Starting DEJU analysis for DEJU-edgeR ...")
DEJU_fit <- DEU_analysis(IE_J$IE_J_count,
                         IE_J$IE_J_annot,
                         mat$group, mat$design)

message("Running diffSpliceDGE for DEU-edgeR ...")
res <- DEU_res(DEU_fit, mat$design, pair)  

message("Running diffSpliceDGE for DEJU-edgeR ...")
res <- DEU_res(DEJU_fit, mat$design, pair)  
  