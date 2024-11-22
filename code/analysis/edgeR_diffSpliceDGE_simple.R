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
              'IE_J_annot'=IE_J_annot,
              'J_count'=J_count,
              'J_annot'=J_annot,
              'IE_count'=IE_count,
              'IE_annot'=IE_annot))
}

norm_libSize <- function(r_count, annot, group, design, mode, REF) {
  
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
  
  return(y)
}

DEU_analysis <- function(y, design, p) {
  
  message("Fitting GLM-QL model for design matrix ...")
  fit <- glmQLFit(y, design, robust=TRUE)
  
  message("Running diffSpliceDGE ..")
  cmd <- paste("makeContrasts(", p, ",levels=design)", sep="")
  contr <- eval(parse(text=cmd))
  print(contr)
  sp <- diffSpliceDGE(fit, contrast=contr, geneid="GeneID", exonid="Start")
  
  return(list('fit'=fit, 'sp'=sp))
}
