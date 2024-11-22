library(dplyr)
library(Rsubread)

dir.create("../../data/simulation/customized_transcriptome",
           recursive = T,
           showWarnings = F)
setwd("../../data/simulation/customized_transcriptome")

message('Merging and flattening exon regions using flattenGTF ...')
GTF <- "../../../annotation/gencode.vM32.PC.annotation.gtf.gz"
SAF <- flattenGTF(GTF, GTF.featureType="exon", GTF.attrType="gene_id", method="merge")

message('Converting 1-based GTF into 0-based BED ...')
BED <- SAF %>%
  group_by(GeneID) %>%
  # calculate blockSizes and blockStarts of exons
  mutate(Start = Start - 1, blockSizes = End - Start, blockStarts = Start - min(Start)) %>% 
  ungroup()

message('Removing duplicated genes with flattened exons ...')
ds_f <- paste0("../../../annotation/duplicated_sequences.tsv") 
ds <- read.table(ds_f, header=FALSE)
BED <- BED[!(BED$GeneID %in% ds$V1),]

message('Extracting genes with more than 2 exons ...')
genes <- table(BED$GeneID) 
BED <- BED[BED$GeneID %in% names(genes[genes>2]), ]

message('Randomly selecting 5000 genes ...')
set.seed(123)
genes <- sample(unique(BED$GeneID), 5000, replace=FALSE)

message('Randomly spliting 5000 genes into 4 equal parts ...')
set.seed(123)
ES <- sample(genes, 1250, replace=FALSE)

message(' For ES pattern ...')
message('Generating full-length isoform ...')
ES_t1 <- BED[BED$GeneID %in% ES, ]
ES_t1 <- cbind(ES_t1, index=1:nrow(ES_t1))
ES_t1 <- ES_t1 %>%
  group_by(GeneID) %>%
  mutate(ExonID=1:n()) %>%
  ungroup()

message('Generating 2nd isoform by randomly skip 1 exon ...')
set.seed(123)
exon_ES <- ES_t1 %>%
  group_by(GeneID) %>%
  sample_n(1) %>%
  ungroup()
ES_t2 <- ES_t1[!(ES_t1$index %in% exon_ES$index), ]
# recalculate blockStarts
ES_t2 <- ES_t2 %>%
  group_by(GeneID) %>%
  mutate(blockStarts = Start - min(Start)) %>%
  ungroup()

message('For MXE pattern ...')
set.seed(123)
MXE <- sample(genes[!(genes %in% ES)], 1250, replace=FALSE)

MXE_genes <- BED[BED$GeneID %in% MXE, ]
MXE_genes <- cbind(MXE_genes, index=1:nrow(MXE_genes))
MXE_genes <- MXE_genes %>%
  group_by(GeneID) %>%
  mutate(ExonID=1:n()) %>%
  ungroup()

message('Generating 2nd isoform by randomly skip 1 exon ...')
set.seed(123)
exon_MXE_t1 <- MXE_genes %>%
  group_by(GeneID) %>%
  sample_n(1) %>%
  ungroup()
MXE_t1 <- MXE_genes[!(MXE_genes$index %in% exon_MXE_t1$index), ]
# recalculate blockStarts
MXE_t1 <- MXE_t1 %>%
  group_by(GeneID) %>%
  mutate(blockStarts = Start - min(Start)) %>%
  ungroup()

message('Generating 2nd isoform by randomly skip 1 exon that mutually exclusive ...')
set.seed(123)
exon_MXE_t2 <- MXE_t1 %>%
  group_by(GeneID) %>%
  sample_n(1) %>%
  ungroup()
MXE_t2 <- MXE_genes[!(MXE_genes$index %in% exon_MXE_t2$index), ]
# recalculate blockStarts
MXE_t2 <- MXE_t2 %>%
  group_by(GeneID) %>%
  mutate(blockStarts = Start - min(Start)) %>%
  ungroup()

message('For ASS pattern ...')
ASS <- sample(genes[!(genes %in% union(ES, MXE))], 1250, replace=FALSE)

message('Generating full-length isoform ...')
ASS_t1 <- BED[BED$GeneID %in% ASS, ]
ASS_t1 <- cbind(ASS_t1, index=1:nrow(ASS_t1))
ASS_t1 <- ASS_t1 %>%
  group_by(GeneID) %>%
  mutate(ExonID=1:n()) %>%
  ungroup()

message('Generating 2nd isoform by randomly creating alternative splice site ...')
# Define a function to randomly select one exon and create alternative splice site
select_and_modify_exon <- function(gene_data) {
  if (nrow(gene_data)==3){
    selected_row <- 2
  } else {
    selected_row <- sample(2:(nrow(gene_data)-1), 1)  # Select one random row except first, last row
  }
  selected_exon <- gene_data[selected_row, ]
  random_fraction <- sample(c(0.3, 0.4, 0.5, 0.6, 0.7), 1)  # Random fraction
  delta <- selected_exon$blockSizes * random_fraction
  selected_exon$Start <- round(selected_exon$Start + delta)
  # Recalculate blockSizes and blockStarts
  selected_exon$blockSizes <- selected_exon$End - selected_exon$Start
  selected_exon$blockStarts <- selected_exon$Start - gene_data[1,]$Start
  gene_data[gene_data$End == selected_exon$End, ] <- selected_exon
  return(list(gene_data, selected_exon))
}

# Create an empty data frame to store the modified exons
ASS_t2 <- data.frame()
exon_ASS <- data.frame()
# Loop through unique genes and select/modify one exon per gene
unique_genes <- unique(ASS_t1$GeneID)
for (gene in unique_genes) {
  set.seed(123)
  gene_data <- ASS_t1[ASS_t1$GeneID == gene, ]
  gene_exon_list <- select_and_modify_exon(gene_data)
  modified_gene <- gene_exon_list[[1]]
  selected_exon <- gene_exon_list[[2]]
  ASS_t2 <- rbind(ASS_t2, modified_gene)
  exon_ASS <- rbind(exon_ASS, selected_exon)
}

message('For RI pattern ...')
RI <- genes[!(genes %in% c(ES, MXE, ASS))]

message('Generating full-length isoform ...')
RI_t1 <- BED[BED$GeneID %in% RI, ]
RI_t1 <- cbind(RI_t1, index=1:nrow(RI_t1))
RI_t1 <- RI_t1 %>%
  group_by(GeneID) %>%
  mutate(ExonID=1:n()) %>%
  ungroup()

message('Generating 2nd isoform with 1 intron-retained exon ...')
set.seed(123)
# Define a function to randomly select one exon (except the last exon) and merge with the downstream adjacent exon
select_and_modify_exon <- function(gene_data) {
  selected_row <- sample(1:(nrow(gene_data)-1), 1)
  selected_exon_1 <- gene_data[selected_row, ]
  selected_exon_2 <- gene_data[selected_row + 1, ]
  selected_exon_1$End <- selected_exon_2$End
  # Recalculate blockSizes
  selected_exon_1$blockSizes <- selected_exon_1$End - selected_exon_1$Start
  gene_data[gene_data$Start == selected_exon_1$Start, ] <- selected_exon_1
  gene_data <- gene_data[-(selected_row + 1), ]
  return(list(gene_data, selected_exon_1))
}

# Create an empty data frame to store the modified exons
RI_t2 <- data.frame()
exon_RI <- data.frame()
# Loop through unique genes and select/modify one exon per gene
unique_genes <- unique(RI_t1$GeneID)
for (gene in unique_genes) {
  gene_data <- RI_t1[RI_t1$GeneID == gene, ]
  gene_exon_list <- select_and_modify_exon(gene_data)
  modified_gene <- gene_exon_list[[1]]
  selected_exon <- gene_exon_list[[2]]
  RI_t2 <- rbind(RI_t2, modified_gene)
  exon_RI <- rbind(exon_RI, selected_exon)
}

message('Save metadata ...')
write.table(exon_ES, "ES_exon.info.tsv", sep="\t", quote=FALSE, row.names=FALSE)
write.table(ES_t1, "ES_t1.pre.bed12.tsv", sep="\t", quote=FALSE, row.names=FALSE)
write.table(ES_t2, "ES_t2.pre.bed12.tsv", sep="\t", quote=FALSE, row.names=FALSE)

write.table(exon_MXE_t1, "MXE_exon_t1.info.tsv", sep="\t", quote=FALSE, row.names=FALSE)
write.table(MXE_t1, "MXE_t1.pre.bed12.tsv", sep="\t", quote=FALSE, row.names=FALSE)
write.table(exon_MXE_t2, "MXE_exon_t2.info.tsv", sep="\t", quote=FALSE, row.names=FALSE)
write.table(MXE_t2, "MXE_t2.pre.bed12.tsv", sep="\t", quote=FALSE, row.names=FALSE)

write.table(exon_ASS, "ASS_exon.info.tsv", sep="\t", quote=FALSE, row.names=FALSE)
write.table(ASS_t1, "ASS_t1.pre.bed12.tsv", sep="\t", quote=FALSE, row.names=FALSE)
write.table(ASS_t2, "ASS_t2.pre.bed12.tsv", sep="\t", quote=FALSE, row.names=FALSE)

write.table(exon_RI, "RI_exon.info.tsv", sep="\t", quote=FALSE, row.names=FALSE)
write.table(RI_t1, "RI_t1.pre.bed12.tsv", sep="\t", quote=FALSE, row.names=FALSE)
write.table(RI_t2, "RI_t2.pre.bed12.tsv", sep="\t", quote=FALSE, row.names=FALSE)

AS_info <- data.frame(gene=c(ES, MXE, ASS, RI), 
                      type=c(rep("ES", length(ES)), 
                             rep("MXE", length(MXE)), 
                             rep("ASS", length(ASS)), 
                             rep("RI", length(RI))))
write.table(AS_info, "DEU_genes.info.tsv", sep="\t", quote=FALSE, row.names=FALSE)
