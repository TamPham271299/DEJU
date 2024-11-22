OK <- requireNamespace("devtools", quietly = TRUE)
if (!OK) {
  stop("devtools package required but is not installed (or can't be loaded)")
}
library(Rsubread)
library(edgeR)
library(BiocParallel)

# Getting parameters
args <- R.utils::commandArgs(asValues = TRUE)
print(args)

DIR <- as.character(args[['DIR']])
REF <- as.character(args[['REF']])
targetp <- as.character(args[['target']])

path_to_sj <- paste0(REF, "gencode.vM32.primary_assembly.annotation.uniqSJInfo.tsv")
num_of_occur_sj <- paste0(REF, "gencode.vM32.primary_assembly.annotation.uniqSJInfo.numOfOccur.tsv")
i <- "S1"

### 0. Check parameter
message("Check valid parameters")
print(DIR)
print(REF)
print(path_to_sj)
print(num_of_occur_sj)
print(targetp)


### 02. Run edgeR
message('Start running edgeR ...')
  
message(paste0('Processing: ', i, " ..."))
OUTPUT <- paste0(DIR, "STAR_output/edgeR_fit_RL_75_test/", i, "/")
dir.create(OUTPUT, showWarnings = FALSE, recursive = TRUE)
readCount_output <- paste0(DIR, "STAR_output/featureCounts_fit_RL_75/minMQS_255/", i, "/")

### 021. Load target
target <- read.table(targetp, sep="\t", header=TRUE, stringsAsFactors=FALSE)
group <- factor(target$Group, levels=c("Group_1", "Group_2"))
design <- model.matrix(~ group)

### 022. Loading read counts
message('Loading read counts ...')

# Exon-level read counts with double-counted junctions
exon_count_f <- read.table(paste0(readCount_output, "exon_count.tsv"), header=TRUE)
exon_count <- data.frame(exon_count_f[,1:ncol(exon_count_f)], row.names=NULL, check.names = FALSE)
exon_annot <- read.table(paste0(readCount_output, "exon_count_annotation.tsv"), header=TRUE)

# Internal exon read counts
internal_exon_count_f <- read.table(paste0(readCount_output, "internal_exon_count.tsv"), header=TRUE)
internal_exon_count <- data.frame(internal_exon_count_f[,1:ncol(internal_exon_count_f)], row.names=NULL, check.names = FALSE)
internal_exon_annot <- read.table(paste0(readCount_output, "internal_exon_count_annotation.tsv"), header=TRUE)

# Junction read counts
junc_count_f <- read.table(paste0(readCount_output, "junction_count.tsv"), header=TRUE)
junc_count_f$juncID <- paste(junc_count_f$Site1_chr, junc_count_f$Site1_location, junc_count_f$Site2_location, sep="_")
annotated_sj <- read.table(path_to_sj, header=TRUE)
numOfOccur_sj <- read.table(num_of_occur_sj, header=TRUE)
m_numOfOccur <- match(annotated_sj$juncID, numOfOccur_sj$juncID)
annotated_sj$numOfOccur <- numOfOccur_sj$numOfOccur[m_numOfOccur]
annotated_sj_numOfOccur_1 <- annotated_sj[annotated_sj$numOfOccur==1,]
m1 <- match(junc_count_f$juncID, annotated_sj_numOfOccur_1$juncID)
junc_count_f$PrimaryGene <- ifelse(!is.na(m1), annotated_sj_numOfOccur_1$geneID[m1], junc_count_f$PrimaryGene)
m2 <- match(junc_count_f$juncID, annotated_sj$juncID)
junc_count_f$annotated <- ifelse(!is.na(m2), 1, 0)
junc_count <- data.frame(junc_count_f[, 9:(ncol(junc_count_f)-2)], row.names = NULL, check.names = FALSE)
junc_annot <- data.frame(
  GeneID=junc_count_f$PrimaryGene,
  Chr=junc_count_f$Site1_chr,
  Start=junc_count_f$Site1_location,
  End=junc_count_f$Site2_location
)
m <- match(junc_annot$GeneID, internal_exon_annot$GeneID)
Strand <- internal_exon_annot$Strand[m]
junc_annot <- cbind(junc_annot, Strand=Strand, Length=1, Region="Junction", annotated=junc_count_f$annotated)

internal_exon_junc_count <- rbind(internal_exon_count, junc_count)
internal_exon_annot <- cbind(internal_exon_annot, Region="Exon", annotated=1)
internal_exon_junc_annot <- rbind(internal_exon_annot, junc_annot)

### 023. Preliminary analysis
message('Preliminary analysis ...')

## Existing approach
y_current_ap <- DGEList(counts=exon_count, genes=exon_annot, group=group)
colnames(y_current_ap) <- gsub("[.].*$", "", colnames(y_current_ap))
keep <- filterByExpr(y_current_ap, group=group)
y_current_ap <- y_current_ap[keep, , keep.lib.sizes=FALSE]
y_current_ap <- normLibSizes(y_current_ap)
y_current_ap <- estimateDisp(y_current_ap, design, robust=TRUE)
print(y_current_ap)
fit_current_ap <- glmQLFit(y_current_ap, design, robust=TRUE)

## Junction approach (internal exon read + junction read)
y_new_ap <- DGEList(counts=internal_exon_junc_count, genes=internal_exon_junc_annot, group=group)
colnames(y_new_ap) <- gsub("[.].*$", "", colnames(y_new_ap))
keep <- filterByExpr(y_new_ap, group=group)
y_new_ap <- y_new_ap[keep, , keep.lib.sizes=FALSE]
y_new_ap <- normLibSizes(y_new_ap)
y_new_ap <- estimateDisp(y_new_ap, design, robust=TRUE)
print(y_new_ap)
fit_new_ap <- glmQLFit(y_new_ap, design, robust=TRUE)

### 024. Differential alternative splicing analysis
message('DEU analysis ...')

## Existing approach
sp_current_ap <- diffSpliceDGE(fit_current_ap, coef=2, geneid="GeneID", exonid="Start")
current_ap_simes <- topSpliceDGE(sp_current_ap, test="Simes", n=NULL, FDR=0.05)
current_ap_gene <- topSpliceDGE(sp_current_ap, test="gene", n=NULL, FDR=0.05)
write.table(current_ap_simes, file.path(OUTPUT, "DEU_current_approach_simes_minMQS_255.tsv"), row.names=FALSE, sep="\t", quote=FALSE)
write.table(current_ap_gene, file.path(OUTPUT, "DEU_current_approach_gene_minMQS_255.tsv"), row.names=FALSE, sep="\t", quote=FALSE)

## Junction approach
sp_new_ap <- diffSpliceDGE(fit_new_ap, coef=2, geneid="GeneID", exonid="Start")
new_ap_simes <- topSpliceDGE(sp_new_ap, test="Simes", n=NULL, FDR=0.05)
new_ap_gene <- topSpliceDGE(sp_new_ap, test="gene", n=NULL, FDR=0.05)
write.table(new_ap_simes, file.path(OUTPUT, "DEU_new_approach_simes_minMQS_255_SJ.tsv"), row.names=FALSE, sep="\t", quote=FALSE)
write.table(new_ap_gene, file.path(OUTPUT, "DEU_new_approach_gene_minMQS_255_SJ.tsv"), row.names=FALSE, sep="\t", quote=FALSE)
