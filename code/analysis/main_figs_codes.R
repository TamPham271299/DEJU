library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr)
library(gridExtra)
library(edgeR)
library(Gviz)
library(VennDiagram)
library(UpSetR)
library(stringr)

options(ucscChromosomeNames=FALSE)

# setwd("/vast/projects/Spatial/tam/Differential_splicing/github/code/analysis/")
fig <- "../../figures/main/"
dir.create(fig, recursive = T, showWarnings = F)

################################################################################
###### Codes to produce results in Figure 2
################################################################################
simDir <- "../../data/simulation/"

n <- 3
ls <- "balanced"
fc <- 3
rl <- 75
s <- paste(ls, n, rl, fc, sep="_")

methods <- c("DEJU-edgeR-simes", 
             "DEU-edgeR-simes", 
             "DEJU-limma-simes", 
             "DEU-limma-simes", 
             "DEXSeq", 
             "JunctionSeq")

### Relable method names
methods.lb <-  gsub("-simes", "", methods)
res <- read.csv(paste0(simDir, s, "/performance_analysis/benchmarking_result_sp.mean.csv"), header=TRUE)
res <- res[grep(paste(methods, collapse="|"), res$Method),]
labels <- setNames(methods.lb, methods)

sp <- c("ES", "MXE", "ASS")

fig2 <- list()
for (i in sp) {
  df <- res[grep(i, res$Splicing_pattern), 
                  c("Method", "TP.mean", "FP.mean", "FDR.mean")]
  df <- reshape2::melt(df, id.vars = c("Method", "FDR.mean"))
  df$variable <- factor(df$variable, levels=c("FP.mean", "TP.mean"))
  df$FDR.mean <- ifelse(df$variable == "FP.mean", df$FDR.mean, NA)
  df$Method <- factor(df$Method, levels=methods)
  
  fig2[[i]] <- ggplot(df, aes(fill=variable, y=value, x=Method)) +
    geom_bar(stat="identity", color = "black") +
    scale_x_discrete(labels = labels) +
    labs(y="DEU genes", x=NULL, title = i) +
    ylim(0, 350) +
    scale_fill_manual(values = c("FP.mean" = "red", "TP.mean" = "grey")) +
    theme_test() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color="black"),
          axis.text.y = element_text(color="black"),
          plot.title = element_text(size = 11, hjust = 0.5, vjust = 2),
          axis.title.x = element_text(size=11)) +
    geom_text(aes(label=round(FDR.mean,2)), 
              position = position_stack(vjust = 1), 
              vjust = -0.5)
}

df <- res[grep("RI", res$Splicing_pattern), 
                c("Method", "TP.mean", "FP.mean", "FDR.mean")]
df <- reshape2::melt(df, id.vars = c("Method", "FDR.mean"))
df$variable <- factor(df$variable, levels=c("FP.mean", "TP.mean"))
df$FDR.mean <- ifelse(df$variable == "FP.mean", df$FDR.mean, NA)
df$Method <- factor(df$Method, levels=methods)

fig2[["RI"]] <- ggplot(df, aes(fill=variable, y=value, x=Method)) + 
  geom_bar(stat="identity", color = "black") +
  scale_x_discrete(labels = labels) +
  labs(y="DEU genes", x=NULL, title = "IR") +
  ylim(0, 350) +
  scale_fill_manual(values = c("FP.mean" = "red", "TP.mean" = "grey"), 
                    labels = c("FP.mean" = "FP", "TP.mean" = "TP")) +
  theme_test() +
  theme(legend.position = c(0.95, 0.95),
        legend.justification = c(1, 1),
        legend.title = element_blank(),
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text(size = 9),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color="black"),
        axis.text.y = element_text(color="black"),
        plot.title = element_text(size = 11, hjust = 0.5, vjust = 2),
        axis.title.x = element_text(size=11)) +
  geom_text(aes(label=round(FDR.mean,2)), position = position_stack(vjust = 1), vjust = -0.5)

pdf(paste0(fig, "fig_2.pdf"), height = 5.5, width = 6)
grid.arrange(grobs=fig2, nrow=2, ncol=2, widths=c(2,2))
dev.off()

################################################################################
###### Codes to produce results in Figure 3
################################################################################
simDir <- "../../data/simulation/"

methods <- c("DEJU-edgeR", "DEU-edgeR", 
             "DEJU-limma", "DEU-limma", 
             "DEXSeq", "JunctionSeq")
ngenes <- 5000
nmethods <- length(methods)
nsim <- 20

### Reference gene list
map_f <- "../../data/simulation/customized_transcriptome/DEU_genes.info.tsv"
map <- read.table(map_f, header=TRUE)

### Scenarios
nlibs <- c(3,5,10)
libSize <- c("balanced", "unbalanced")
fc <- 3
rl <- 75

fd_l <- list()
fdmax_l <- list()
nd_l <- list()

for (ls in libSize) {
  for (n in nlibs) {
    s <- paste(ls, n, rl, fc, sep = "_")
    DIR <- paste0(simDir, s, "/")
    print(DIR)
    
    fd <- ranking <- matrix(0, nrow=ngenes, ncol=nmethods)
    colnames(fd) <- colnames(ranking) <- methods
    nd <- rep(0, nmethods)
    names(nd) <- methods    
    
    # BEGIN SIM
    for (simID in paste0("S", 1:nsim)) {
      print(paste0("SIM = ", simID))
      
      simReads_o <- paste0(DIR, "simReads/", simID, "/")
      DEU <- read.table(paste0(simReads_o, "DEU.tsv"), header=TRUE)
      ref.DEUs <- unique(gsub("\\.(skip_exon|splice_site_at_exon|retained_intron_).*", "", DEU$GeneID))
      map$status <- ifelse(map$gene %in% ref.DEUs, 1, 0)
      o <- order(map$status, decreasing = TRUE)
      map <- map[o, ]
      
      # DEJU-edgeR (simes)
      res_f <- paste0(DIR, "edgeR_diffSpliceDGE/", 
                      simID, "/DEJU_simes_test.allGenes.tsv")
      res <- read.table(res_f, header = TRUE)
      FDR <- res$FDR
      i <- 1
      res.1 <- merge(map, res, by.x="gene", by.y="GeneID", all.x=TRUE)
      status <- res.1[order(res.1$P.Value), "status"]
      ranking[1:ngenes,i] <- status
      nd[i] <- nd[i] + sum(FDR<0.05)
      
      # DEU-edgeR (simes)
      res_f <- paste0(DIR, "edgeR_diffSpliceDGE/", 
                      simID, "/DEU_simes_test.allGenes.tsv")
      res <- read.table(res_f, header = TRUE)
      FDR <- res$FDR
      i <- 2
      res.1 <- merge(map, res, by.x="gene", by.y="GeneID", all.x=TRUE)
      status <- res.1[order(res.1$P.Value), "status"]
      ranking[1:ngenes,i] <- status
      nd[i] <- nd[i] + sum(FDR<0.05)
      
      # DEJU-limma (simes)
      res_f <- paste0(DIR, "limma_diffSplice/", 
                      simID, "/DEJU_simes_test.allGenes.tsv")
      res <- read.table(res_f, header = TRUE)
      FDR <- res$FDR
      i <- 3
      res.1 <- merge(map, res, by.x="gene", by.y="GeneID", all.x=TRUE)
      status <- res.1[order(res.1$P.Value), "status"]
      ranking[1:ngenes,i] <- status
      nd[i] <- nd[i] + sum(FDR<0.05)
      
      # DEU-limma (simes)
      res_f <- paste0(DIR, "limma_diffSplice/", 
                      simID, "/DEU_simes_test.allGenes.tsv")
      res <- read.table(res_f, header = TRUE)
      FDR <- res$FDR
      i <- 4
      res.1 <- merge(map, res, by.x="gene", by.y="GeneID", all.x=TRUE)
      status <- res.1[order(res.1$P.Value), "status"]
      ranking[1:ngenes,i] <- status
      nd[i] <- nd[i] + sum(FDR<0.05)
      
      # DEXSeq
      res_f <- paste0(DIR, "DEXSeq/", 
                      simID, "/DEUs.geneWise.allGenes.tsv")
      res <- read.table(res_f, header = TRUE)
      FDR <- res$padj.gene
      i <- 5
      res.1 <- merge(map, res, by.x="gene", by.y="groupID", all.x=TRUE)
      status <- res.1[order(res.1$padj.gene), "status"]
      ranking[1:ngenes,i] <- status
      nd[i] <- nd[i] + sum(FDR<0.05)
      
      # JunctionSeq
      res_f <- paste0(DIR, "JunctionSeq/",
                      simID, "/allGenes.results.txt.gz")
      res <- read.table(res_f, header = TRUE)
      res <- res[, c("geneID", "geneWisePadj")]
      res <- unique(res[complete.cases(res), ])
      FDR <- res$geneWisePadj
      i <- 6
      res.1 <- merge(map, res, by.x="gene", by.y="geneID", all.x=TRUE)
      status <- res.1[order(res.1$geneWisePadj), "status"]
      ranking[1:ngenes,i] <- status
      nd[i] <- nd[i] + sum(FDR<0.05)
      
      # DEJU-edgeR (F)
      res_f <- paste0(DIR, "edgeR_diffSpliceDGE/", 
                      simID, "/DEJU_F_test.allGenes.tsv")
      res <- read.table(res_f, header = TRUE)
      FDR <- res$FDR
      i <- 7
      res.1 <- merge(map, res, by.x="gene", by.y="GeneID", all.x=TRUE)
      status <- res.1[order(res.1$P.Value), "status"]
      ranking[1:ngenes,i] <- status
      nd[i] <- nd[i] + sum(FDR<0.05)
      
      # DEU-edgeR (F)
      res_f <- paste0(DIR, "edgeR_diffSpliceDGE/", 
                      simID, "/DEU_F_test.allGenes.tsv")
      res <- read.table(res_f, header = TRUE)
      FDR <- res$FDR
      i <- 8
      res.1 <- merge(map, res, by.x="gene", by.y="GeneID", all.x=TRUE)
      status <- res.1[order(res.1$P.Value), "status"]
      ranking[1:ngenes,i] <- status
      nd[i] <- nd[i] + sum(FDR<0.05)
      
      # DEJU-limma (F)
      res_f <- paste0(DIR, "limma_diffSplice/", 
                      simID, "/DEJU_F_test.allGenes.tsv")
      res <- read.table(res_f, header = TRUE)
      FDR <- res$FDR
      i <- 9
      res.1 <- merge(map, res, by.x="gene", by.y="GeneID", all.x=TRUE)
      status <- res.1[order(res.1$P.Value), "status"]
      ranking[1:ngenes,i] <- status
      nd[i] <- nd[i] + sum(FDR<0.05)
      
      # DEU-limma (F)
      res_f <- paste0(DIR, "limma_diffSplice/", 
                      simID, "/DEU_F_test.allGenes.tsv")
      res <- read.table(res_f, header = TRUE)
      FDR <- res$FDR
      i <- 10
      res.1 <- merge(map, res, by.x="gene", by.y="GeneID", all.x=TRUE)
      status <- res.1[order(res.1$P.Value), "status"]
      ranking[1:ngenes,i] <- status
      nd[i] <- nd[i] + sum(FDR<0.05)
      
      # Calculate cumulative fd
      fd <- fd + apply(1-abs(ranking),2,cumsum)
    }
    fd_l[[s]] <- fd/nsim
    fdmax_l[[s]] <- apply(fd_l[[s]],1,max)
    nd_l[[s]] <- nd/nsim
  }
}

i <- 1:1000
fd.col <- c("green","green",
            "blue", "blue", 
            "violet", "black")
fd.type <- c(1,2,1,2,1,1)

pdf(paste0(fig, "fig_3.pdf"), height =7, width = 7)
par(mfrow=c(2,2))

s <- "balanced_3_75_3"
plot(i, fdmax_l[[s]][i], type="n", xlab="Genes chosen", ylab="False discoveries", main="Balanced, n=3")
lines(i,fd_l[[s]][i,1],col="green",lwd=2)
lines(i,fd_l[[s]][i,2],col="green",lwd=2, lty=2)
lines(i,fd_l[[s]][i,3],col="blue",lwd=2)
lines(i,fd_l[[s]][i,4],col="blue",lwd=2, lty=2)
lines(i,fd_l[[s]][i,5],col="violet",lwd=2)
lines(i,fd_l[[s]][i,6],col="black",lwd=2)
order <- 1:6
legend("topleft", legend=methods[order], lwd=2, col=fd.col[order], lty=fd.type[order], bty="n", cex=1)

s <- "unbalanced_3_75_3"
plot(i, fdmax_l[[s]][i], type="n", xlab="Genes chosen", ylab="False discoveries", main="Unbalanced, n=3")
lines(i,fd_l[[s]][i,1],col="green",lwd=2)
lines(i,fd_l[[s]][i,2],col="green",lwd=2, lty=2)
lines(i,fd_l[[s]][i,3],col="blue",lwd=2)
lines(i,fd_l[[s]][i,4],col="blue",lwd=2, lty=2)
lines(i,fd_l[[s]][i,5],col="violet",lwd=2)
lines(i,fd_l[[s]][i,6],col="black",lwd=2)
order <- 1:6

s <- "balanced_5_75_3"
plot(i, fdmax_l[[s]][i], type="n", xlab="Genes chosen", ylab="False discoveries", main="Balanced, n=5")
lines(i,fd_l[[s]][i,1],col="green",lwd=2)
lines(i,fd_l[[s]][i,2],col="green",lwd=2, lty=2)
lines(i,fd_l[[s]][i,3],col="blue",lwd=2)
lines(i,fd_l[[s]][i,4],col="blue",lwd=2, lty=2)
lines(i,fd_l[[s]][i,5],col="violet",lwd=2)
lines(i,fd_l[[s]][i,6],col="black",lwd=2)
order <- 1:6

s <- "unbalanced_5_75_3"
plot(i, fdmax_l[[s]][i], type="n", xlab="Genes chosen", ylab="False discoveries", main="Unbalanced, n=5")
lines(i,fd_l[[s]][i,1],col="green",lwd=2)
lines(i,fd_l[[s]][i,2],col="green",lwd=2, lty=2)
lines(i,fd_l[[s]][i,3],col="blue",lwd=2)
lines(i,fd_l[[s]][i,4],col="blue",lwd=2, lty=2)
lines(i,fd_l[[s]][i,5],col="violet",lwd=2)
lines(i,fd_l[[s]][i,6],col="black",lwd=2)
order <- 1:6

dev.off()

################################################################################
###### Codes to produce results in Figure 4A and 4B
################################################################################
setwd("/vast/projects/Spatial/tam/Differential_splicing/github/code/analysis/")
source("edgeR_diffSpliceDGE_simple.R")
source("plotJunc3_diffSpliceDGE.R")

# DIR <- "/vast/projects/MM/tam/Differential_splicing/milevskiy_2023_GSE227748/"
DIR <- "../../data/case_study/GSE227748/"
p <- "LP-ML"
# targetp <- paste0(DIR, "target/target.tsv")
targetp <- paste0(DIR, "target/target.", p, ".tsv")
featureCounts_o <- paste0(DIR, "featureCounts/")
REF <- "../../annotation/"
SJ <- paste0(REF, "gencode.vM32.primary_assembly.annotation.SJdatabase.tsv")
SJ_database <- read.table(SJ, header=TRUE)

message("Constructing design matrix for contrast ...")
mat <- construct_model_matrix(targetp)

message("Constructing a matrix for exon read counts with duplicated junction count ...")
E <- process_E_count_withDupJ(featureCounts_o)

message("Combining internal exon and junction read count matrix ...")  
IE_J <- process_IE_J_count(featureCounts_o, SJ_database)

message("Analyze libSize for DEU-edgeR ...")
E_y <- norm_libSize(E$E_count,
                    E$E_annot,
                    mat$group, mat$design,
                    "case_study", REF)

message("Analyze libSize for DEJU-edgeR ...")
IE_J_y <- norm_libSize(IE_J$IE_J_count,
                       IE_J$IE_J_annot,
                       mat$group, mat$design,
                       "case_study", REF)

message("Start DEU analysis for DEU-edgeR ...")
DEU <- DEU_analysis(E_y, mat$design, p)

message("Start DEJU analysis for DEJU-edgeR ...")
DEJU <- DEU_analysis(IE_J_y, mat$design, p)

### Figure 4A
message(paste0("Make VennDiagram for ", p, " ..."))
sigDEUgenes <- list()
# DEJU-edgeR-simes
res_f <- paste0(DIR, "edgeR_diffSpliceDGE/", 
                p, "/DEJU_simes_test.sigGenes_0.05.tsv")
res <- read.table(res_f, header = TRUE)
i <- 1
sigDEUgenes[[i]] <- res$GeneID

# DEU-edgeR-simes
res_f <- paste0(DIR, "edgeR_diffSpliceDGE/", 
                p, "/DEU_simes_test.sigGenes_0.05.tsv")
res <- read.table(res_f, header = TRUE)
i <- 2
sigDEUgenes[[i]] <- res$GeneID

# DEJU-edgeR-F
res_f <- paste0(DIR, "edgeR_diffSpliceDGE/", 
                p, "/DEJU_F_test.sigGenes_0.05.tsv")
res <- read.table(res_f, header = TRUE)
i <- 3
sigDEUgenes[[i]] <- res$GeneID

# DEU-edgeR-F
res_f <- paste0(DIR, "edgeR_diffSpliceDGE/", 
                p, "/DEU_F_test.sigGenes_0.05.tsv")
res <- read.table(res_f, header = TRUE)
i <- 4
sigDEUgenes[[i]] <- res$GeneID

l_simes <- list("DEJU-edgeR" = sigDEUgenes[[1]],
                "DEU-edgeR" = sigDEUgenes[[2]])

l_F <- list("DEJU-edgeR" = sigDEUgenes[[3]],
            "DEU-edgeR" = sigDEUgenes[[4]])

colors <- c("blue", "grey")
v_simes <- venn.diagram(l_simes, filename = NULL, disable.logging = TRUE, 
             main=NULL, print.mode=c("raw", "percent"),
             fill = colors, lwd=1)
v_F <- venn.diagram(l_F, filename = NULL, disable.logging = TRUE, 
                        main=NULL, print.mode=c("raw", "percent"),
                        fill = colors, lwd=1)
pdf(paste0(fig, "fig_4A.pdf"), height = 5, width = 2.5)
grid.arrange(v_simes, v_F, nrow=2, ncol=1)
dev.off()

# Fgfr1 gene
g <- "Fgfr1"
message(paste0("Make plotJunc plot for ", g, " ..."))
pdf(paste0(fig, "fig_4B_1.", g, ".pdf"), height = 8.5, width = 8.5)
par(mfrow=c(2,1))
plotJunc(DEU$sp, geneid=g, genecol="Symbol", annotation=IE_J$IE_annot)
plotJunc(DEJU$sp, geneid=g, genecol="Symbol", annotation=IE_J$IE_annot)
dev.off()

message(paste0("Make sashimi plot for ", g, " ..."))
geneSymbol <- "Fgfr1"
geneID <- "ENSMUSG00000031565.19"
chr <- "chr8"
g_from <- 26002669
g_to <- 26066734
d <- 1000

message("Load alignmnent track from bam files generated using make_bam_for_visualization.sh ...")
samples <- c("LP_rep1", "LP_rep2", "LP_rep3", "ML_rep1", "ML_rep2", "ML_rep3")
OUTPUT <- paste0(DIR, "Gviz/", geneSymbol, "_", geneID, "_", d, "/")
PE_bam_files <- paste0(OUTPUT, samples, ".", geneSymbol, ".bam")

alTrack <- list()

for (idx in 1:3) {
  s <- samples[idx]
  alTrack[[s]] <- AlignmentsTrack(PE_bam_files[idx], 
                                  isPaired = TRUE, 
                                  fill.coverage="blue", 
                                  col.sashimi="blue", 
                                  fill="white", 
                                  fontsize=8, cex.axis=1, col="blue")
}

for (idx in 4:6) {
  s <- samples[idx]
  alTrack[[s]] <- AlignmentsTrack(PE_bam_files[idx], 
                                  isPaired = TRUE, 
                                  fill.coverage="green", 
                                  col.sashimi="green", 
                                  fill="white", 
                                  fontsize=8, cex.axis=1, col="green")
}

message("Load transcript annotation track ...")
knownGenes <- UcscTrack(genome = "mm39", chromosome = chr, 
                        track = "All GENCODE VM32", table="wgEncodeGencodeCompVM32", 
                        from = g_from, to = g_to,
                        trackType = "GeneRegionTrack", 
                        rstarts = "exonStarts", rends = "exonEnds", 
                        gene = "name", symbol = "name", 
                        transcript = "name", strand = "strand", 
                        fill = "black", name = "ALL GENCODE VM32", col="black")

message("Load genome axis track ...")
idxTrack <- IdeogramTrack(genome="mm39", chromosome=chr)
axTrack <- GenomeAxisTrack()

message("Combine all tracks ...")
pdf(paste0(fig, "fig_4B_2.pdf"), height = 7, width = 4)
plotTracks(c(idxTrack, axTrack, knownGenes, alTrack),
           from = 26052007, to = 26059369, chromosome = chr, type = c("coverage", "sashimi"), 
           sashimiNumbers=TRUE,
           showTitle = FALSE,
           sashimiScore=9,
           lwd.sashimiMax=2,
           sizes = c(0.1,0.2,0.3,rep(0.3, 6)),
           lwd.title=2,
           lwd.border.title=4,
           background.title="white",
           col.axis="black",
           fontcolor="black",
           col.line="black")
dev.off()

svg(paste0(fig, "fig_4B_2.svg"), height = 7, width = 4)
plotTracks(c(idxTrack, axTrack, knownGenes, alTrack),
           from = 26052007, to = 26059369, chromosome = chr, type = c("coverage", "sashimi"), 
           sashimiNumbers=TRUE,
           showTitle = FALSE,
           sashimiScore=9,
           lwd.sashimiMax=2,
           sizes = c(0.1,0.2,0.3,rep(0.3, 6)),
           lwd.title=2,
           lwd.border.title=4,
           background.title="white",
           col.axis="black",
           fontcolor="black",
           col.line="black")
dev.off()

### Figure 4C
sigDEUgenes <- list()

# DEJU-edgeR
res_f <- paste0(DIR, "edgeR_diffSpliceDGE/", 
                p, "/DEJU_simes_test.sigGenes_0.05.tsv")
res <- read.table(res_f, header = TRUE)
i <- 1
sigDEUgenes[[i]] <- res$GeneID

# DEU-edgeR
res_f <- paste0(DIR, "edgeR_diffSpliceDGE/", 
                p, "/DEU_simes_test.sigGenes_0.05.tsv")
res <- read.table(res_f, header = TRUE)
i <- 2
sigDEUgenes[[i]] <- res$GeneID

# DEJU-limma
res_f <- paste0(DIR, "limma_diffSplice/", 
                p, "/DEJU_simes_test.sigGenes_0.05.tsv")
res <- read.table(res_f, header = TRUE)
i <- 3
sigDEUgenes[[i]] <- res$GeneID

# DEU-limma
res_f <- paste0(DIR, "limma_diffSplice/", 
                p, "/DEU_simes_test.sigGenes_0.05.tsv")
res <- read.table(res_f, header = TRUE)
i <- 4
sigDEUgenes[[i]] <- res$GeneID

# DEXSeq
res_f <- paste0(DIR, "DEXSeq/", 
                p, "/DEUs.geneWise.sigGenes.0.05.withGeneSymbol.tsv")
res <- read.table(res_f, header = TRUE)
i <- 5
sigDEUgenes[[i]] <- res$groupID

# JunctionSeq
res_f <- paste0(DIR, "JunctionSeq/",
                p, "/sigGenes_0.05.genewiseResults.flt.withGeneSymbol.txt")
res <- read.table(res_f, header = TRUE)
i <- 6
sigDEUgenes[[i]] <- res$geneID

# DEJU-edgeR
res_f <- paste0(DIR, "edgeR_diffSpliceDGE/", 
                p, "/DEJU_F_test.sigGenes_0.05.tsv")
res <- read.table(res_f, header = TRUE)
i <- 7
sigDEUgenes[[i]] <- res$GeneID

# DEU-edgeR
res_f <- paste0(DIR, "edgeR_diffSpliceDGE/", 
                p, "/DEU_F_test.sigGenes_0.05.tsv")
res <- read.table(res_f, header = TRUE)
i <- 8
sigDEUgenes[[i]] <- res$GeneID

set_name <- c("DEJU-edgeR (simes)", "DEU-edgeR (simes)", 
              "DEJU-limma", "DEU-limma",
              "DEXSeq", "JunctionSeq",
              "DEJU-edgeR (F)", "DEU-edgeR (F)")
names(sigDEUgenes) <- set_name

pdf(paste0(fig, "fig_4C.pdf"), height = 5, width = 6)
upset(fromList(sigDEUgenes), order.by = "freq", nsets = length(sigDEUgenes), nintersects = 10)
dev.off()

################################################################################
###### Codes to produce results in Table 1
################################################################################
simDir <- "../../data/simulation/"

nlibs <- c(3, 5, 10)
ls <- "balanced"
fc <- 3
rl <- 75
methods <- c("DEJU-edgeR-simes",
             "DEJU-limma-simes")
methods.lb <-  gsub("-simes", "", methods)
labels <- setNames(methods.lb, methods)
sp <- c("ES", "MXE", "ASS", "RI")
# num_row <- sum(length(sp)+length(nlibs)+length(methods))

tab1 <- data.frame(matrix(ncol=5, nrow=0))
colnames(tab1) <- c("Pattern", "n", "Method", "FDR", "Power")
idx <- 1
for (i in sp) {
  for (n in nlibs) {
    for (m in methods) {
      s <- paste(ls, n, rl, fc, sep="_")
      res <- read.csv(paste0(simDir, s, "/performance_analysis/benchmarking_result_sp.mean.csv"), header=TRUE)
      res <- res[res$Method==m & res$Splicing_pattern==i,]
      tab1[idx,] <- c(i, n, m, round(res$FDR.mean,3), round(res$Sensitivity.mean,3)) 
      idx <- idx + 1
    }
  }
}
tab1$Method <- factor(tab1$Method, levels=methods)
tab1$n <- factor(tab1$n, levels=c(3,5,10))
tab1$Pattern <- factor(tab1$Pattern, levels=c("ES", "MXE", "ASS", "RI"))
DEJU.edgeR <- tab1[tab1$Method=="DEJU-edgeR-simes", -3]
DEJU.limma <- tab1[tab1$Method=="DEJU-limma-simes", -3]

tab1 <- merge(DEJU.edgeR, DEJU.limma, by=c("Pattern", "n"), suffixes=c(".DEJU.edgeR", ".DEJU.limma"))
tab1 <- tab1[order(tab1$Pattern, tab1$n), ]
write.csv(tab1, paste0(fig, "table_1.csv"), row.names=FALSE)

################################################################################
###### Codes to produce results in Table 2
################################################################################
simDir <- "../../data/simulation/"

nlibs <- c(3, 5, 10)
ls <- c("balanced", "unbalanced")
fc <- 3
rl <- 75
tools <- c("edgeR_diffSpliceDGE", "limma_diffSplice", "DEXSeq", "JunctionSeq")
idx <- 1:length(tools)

extract_info <- function(f, p) {
  lines <- readLines(f)
  matched_lines <- lines[grep(p, lines)]
  return(matched_lines)
}

tab2 <- data.frame(matrix(ncol=5, nrow=0))
colnames(tab2) <- c("libSize", "n", "Method", "Walltime", "Memory")

for (l in ls) {
  for (n in nlibs) {
    s <- paste(l, n, rl, fc, sep="_")
    print(s)
    logs <- paste0(simDir, s, "/log/", "03_", idx, "_DS_", tools, "_test.log")
    for (i in idx) {
      
      # Extract time-wall (TW)
      TW <- extract_info(logs[i], "wall clock")
      TW <- gsub("\tElapsed \\(wall clock\\) time \\(h:mm:ss or m:ss\\): ", "", TW)
      # Extract peak memory usage (PMU)
      PMU <- extract_info(logs[i], "Maximum resident set size")
      PMU <- gsub("\tMaximum resident set size \\(kbytes\\): ", "", PMU)
      df <- data.frame(libSize=l, n=n, Method=tools[i], Walltime=TW, Memory=PMU)
      tab2 <- rbind(tab2, df)
    }
  }
}

time_convert <- function(TW) {
  TW <- strsplit(TW, ":")[[1]]
  if (length(TW)==2){
    min <- as.numeric(TW[1])
    sec <- as.numeric(TW[2])
    TW <- round(min + sec/60,2)
  } else if (length(TW)==3){
    hr <- as.numeric(TW[1])
    min <- as.numeric(TW[2])
    sec <- as.numeric(TW[3])
    TW <- round(hr*60 + min + sec/60,2)
  }
  return(TW)  
}

tab2$Walltime <- sapply(tab2$Walltime, time_convert)
tab2$Memory <- round(as.numeric(tab2$Memory)/(1024*1024),2)
tab2$mode <- ifelse(tab2$Method=="edgeR_diffSpliceDGE" | tab2$Method=="limma_diffSplice",
                    "serial", NA)

tab2[which(is.na(tab2$mode)), "mode"] <- rep(c("serial", "parallel"), length(which(is.na(tab2$mode)))/2)
write.csv(tab2, paste0(fig, "table_2.csv"), row.names=FALSE)

##############################################################################################################
# library(ggplot2)
# library(reshape2)
# library(plyr)
# library(dplyr)
# library(gridExtra)
# library(edgeR)
# library(Gviz)
# library(VennDiagram)
# library(UpSetR)
# library(stringr)

# options(ucscChromosomeNames=FALSE)

# # setwd("/vast/projects/Spatial/tam/Differential_splicing/github/code/analysis/")
# fig <- "../../figures/main/"
# dir.create(fig, recursive = T, showWarnings = F)

# ################################################################################
# ###### Codes to produce results in Figure 2
# ################################################################################
# simDir <- "../../data/simulation/"

# n <- 3
# ls <- "balanced"
# fc <- 3
# rl <- 75
# s <- paste(ls, n, rl, fc, sep="_")

# methods <- c("DEJU-edgeR-simes", 
#              "DEU-edgeR-simes", 
#              "DEJU-limma-simes", 
#              "DEU-limma-simes", 
#              "DEXSeq", 
#              "JunctionSeq")

# ### Relable method names
# methods.lb <-  gsub("-simes", "", methods)
# res <- read.csv(paste0(simDir, s, "/performance_analysis/benchmarking_result_sp.mean.csv"), header=TRUE)
# res <- res[grep(paste(methods, collapse="|"), res$Method),]
# labels <- setNames(methods.lb, methods)

# sp <- c("ES", "MXE", "ASS")

# fig2 <- list()
# for (i in sp) {
#   df <- res[grep(i, res$Splicing_pattern), 
#                   c("Method", "TP.mean", "FP.mean", "FDR.mean")]
#   df <- reshape2::melt(df, id.vars = c("Method", "FDR.mean"))
#   df$variable <- factor(df$variable, levels=c("FP.mean", "TP.mean"))
#   df$FDR.mean <- ifelse(df$variable == "FP.mean", df$FDR.mean, NA)
#   df$Method <- factor(df$Method, levels=methods)
  
#   fig2[[i]] <- ggplot(df, aes(fill=variable, y=value, x=Method)) +
#     geom_bar(stat="identity", color = "black") +
#     scale_x_discrete(labels = labels) +
#     labs(y="DEU genes", x=NULL, title = i) +
#     ylim(0, 350) +
#     scale_fill_manual(values = c("FP.mean" = "red", "TP.mean" = "grey")) +
#     theme_test() +
#     theme(legend.position = "none",
#           axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color="black"),
#           axis.text.y = element_text(color="black"),
#           plot.title = element_text(size = 11, hjust = 0.5, vjust = 2),
#           axis.title.x = element_text(size=11)) +
#     geom_text(aes(label=round(FDR.mean,2)), 
#               position = position_stack(vjust = 1), 
#               vjust = -0.5)
# }

# df <- res[grep("RI", res$Splicing_pattern), 
#                 c("Method", "TP.mean", "FP.mean", "FDR.mean")]
# df <- reshape2::melt(df, id.vars = c("Method", "FDR.mean"))
# df$variable <- factor(df$variable, levels=c("FP.mean", "TP.mean"))
# df$FDR.mean <- ifelse(df$variable == "FP.mean", df$FDR.mean, NA)
# df$Method <- factor(df$Method, levels=methods)

# fig2[["RI"]] <- ggplot(df, aes(fill=variable, y=value, x=Method)) + 
#   geom_bar(stat="identity", color = "black") +
#   scale_x_discrete(labels = labels) +
#   labs(y="DEU genes", x=NULL, title = "IR") +
#   ylim(0, 350) +
#   scale_fill_manual(values = c("FP.mean" = "red", "TP.mean" = "grey"), 
#                     labels = c("FP.mean" = "FP", "TP.mean" = "TP")) +
#   theme_test() +
#   theme(legend.position = c(0.95, 0.95),
#         legend.justification = c(1, 1),
#         legend.title = element_blank(),
#         legend.key.size = unit(0.4, "cm"),
#         legend.text = element_text(size = 9),
#         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color="black"),
#         axis.text.y = element_text(color="black"),
#         plot.title = element_text(size = 11, hjust = 0.5, vjust = 2),
#         axis.title.x = element_text(size=11)) +
#   geom_text(aes(label=round(FDR.mean,2)), position = position_stack(vjust = 1), vjust = -0.5)

# pdf(paste0(fig, "fig_2.pdf"), height = 5.5, width = 6)
# grid.arrange(grobs=fig2, nrow=2, ncol=2, widths=c(2,2))
# dev.off()

# ################################################################################
# ###### Codes to produce results in Figure 3
# ################################################################################
# simDir <- "../../data/simulation/"

# methods <- c("DEJU-edgeR", "DEU-edgeR", 
#              "DEJU-limma", "DEU-limma", 
#              "DEXSeq", "JunctionSeq")
# ngenes <- 5000
# nmethods <- length(methods)
# nsim <- 20

# ### Reference gene list
# map_f <- "../../data/simulation/customized_transcriptome/DEU_genes.info.tsv"
# map <- read.table(map_f, header=TRUE)

# ### Scenarios
# nlibs <- c(3,5,10)
# libSize <- c("balanced", "unbalanced")
# fc <- 3
# rl <- 75

# fd_l <- list()
# fdmax_l <- list()
# nd_l <- list()

# for (ls in libSize) {
#   for (n in nlibs) {
#     s <- paste(ls, n, rl, fc, sep = "_")
#     DIR <- paste0(simDir, s, "/")
#     print(DIR)
    
#     fd <- ranking <- matrix(0, nrow=ngenes, ncol=nmethods)
#     colnames(fd) <- colnames(ranking) <- methods
#     nd <- rep(0, nmethods)
#     names(nd) <- methods    
    
#     # BEGIN SIM
#     for (simID in paste0("S", 1:nsim)) {
#       print(paste0("SIM = ", simID))
      
#       simReads_o <- paste0(DIR, "simReads/", simID, "/")
#       DEU <- read.table(paste0(simReads_o, "DEU.tsv"), header=TRUE)
#       ref.DEUs <- unique(gsub("\\.(skip_exon|splice_site_at_exon|retained_intron_).*", "", DEU$GeneID))
#       map$status <- ifelse(map$gene %in% ref.DEUs, 1, 0)
#       o <- order(map$status, decreasing = TRUE)
#       map <- map[o, ]
      
#       # DEJU-edgeR (simes)
#       res_f <- paste0(DIR, "edgeR_diffSpliceDGE/", 
#                       simID, "/DEJU_simes_test.allGenes.tsv")
#       res <- read.table(res_f, header = TRUE)
#       FDR <- res$FDR
#       i <- 1
#       res.1 <- merge(map, res, by.x="gene", by.y="GeneID", all.x=TRUE)
#       status <- res.1[order(res.1$P.Value), "status"]
#       ranking[1:ngenes,i] <- status
#       nd[i] <- nd[i] + sum(FDR<0.05)
      
#       # DEU-edgeR (simes)
#       res_f <- paste0(DIR, "edgeR_diffSpliceDGE/", 
#                       simID, "/DEU_simes_test.allGenes.tsv")
#       res <- read.table(res_f, header = TRUE)
#       FDR <- res$FDR
#       i <- 2
#       res.1 <- merge(map, res, by.x="gene", by.y="GeneID", all.x=TRUE)
#       status <- res.1[order(res.1$P.Value), "status"]
#       ranking[1:ngenes,i] <- status
#       nd[i] <- nd[i] + sum(FDR<0.05)
      
#       # DEJU-limma (simes)
#       res_f <- paste0(DIR, "limma_diffSplice/", 
#                       simID, "/DEJU_simes_test.allGenes.tsv")
#       res <- read.table(res_f, header = TRUE)
#       FDR <- res$FDR
#       i <- 3
#       res.1 <- merge(map, res, by.x="gene", by.y="GeneID", all.x=TRUE)
#       status <- res.1[order(res.1$P.Value), "status"]
#       ranking[1:ngenes,i] <- status
#       nd[i] <- nd[i] + sum(FDR<0.05)
      
#       # DEU-limma (simes)
#       res_f <- paste0(DIR, "limma_diffSplice/", 
#                       simID, "/DEU_simes_test.allGenes.tsv")
#       res <- read.table(res_f, header = TRUE)
#       FDR <- res$FDR
#       i <- 4
#       res.1 <- merge(map, res, by.x="gene", by.y="GeneID", all.x=TRUE)
#       status <- res.1[order(res.1$P.Value), "status"]
#       ranking[1:ngenes,i] <- status
#       nd[i] <- nd[i] + sum(FDR<0.05)
      
#       # DEXSeq
#       res_f <- paste0(DIR, "DEXSeq/", 
#                       simID, "/DEUs.geneWise.allGenes.tsv")
#       res <- read.table(res_f, header = TRUE)
#       FDR <- res$padj.gene
#       i <- 5
#       res.1 <- merge(map, res, by.x="gene", by.y="groupID", all.x=TRUE)
#       status <- res.1[order(res.1$padj.gene), "status"]
#       ranking[1:ngenes,i] <- status
#       nd[i] <- nd[i] + sum(FDR<0.05)
      
#       # JunctionSeq
#       res_f <- paste0(DIR, "JunctionSeq/",
#                       simID, "/allGenes.results.txt.gz")
#       res <- read.table(res_f, header = TRUE)
#       res <- res[, c("geneID", "geneWisePadj")]
#       res <- unique(res[complete.cases(res), ])
#       FDR <- res$geneWisePadj
#       i <- 6
#       res.1 <- merge(map, res, by.x="gene", by.y="geneID", all.x=TRUE)
#       status <- res.1[order(res.1$geneWisePadj), "status"]
#       ranking[1:ngenes,i] <- status
#       nd[i] <- nd[i] + sum(FDR<0.05)
      
#       # DEJU-edgeR (F)
#       res_f <- paste0(DIR, "edgeR_diffSpliceDGE/", 
#                       simID, "/DEJU_F_test.allGenes.tsv")
#       res <- read.table(res_f, header = TRUE)
#       FDR <- res$FDR
#       i <- 7
#       res.1 <- merge(map, res, by.x="gene", by.y="GeneID", all.x=TRUE)
#       status <- res.1[order(res.1$P.Value), "status"]
#       ranking[1:ngenes,i] <- status
#       nd[i] <- nd[i] + sum(FDR<0.05)
      
#       # DEU-edgeR (F)
#       res_f <- paste0(DIR, "edgeR_diffSpliceDGE/", 
#                       simID, "/DEU_F_test.allGenes.tsv")
#       res <- read.table(res_f, header = TRUE)
#       FDR <- res$FDR
#       i <- 8
#       res.1 <- merge(map, res, by.x="gene", by.y="GeneID", all.x=TRUE)
#       status <- res.1[order(res.1$P.Value), "status"]
#       ranking[1:ngenes,i] <- status
#       nd[i] <- nd[i] + sum(FDR<0.05)
      
#       # DEJU-limma (F)
#       res_f <- paste0(DIR, "limma_diffSplice/", 
#                       simID, "/DEJU_F_test.allGenes.tsv")
#       res <- read.table(res_f, header = TRUE)
#       FDR <- res$FDR
#       i <- 9
#       res.1 <- merge(map, res, by.x="gene", by.y="GeneID", all.x=TRUE)
#       status <- res.1[order(res.1$P.Value), "status"]
#       ranking[1:ngenes,i] <- status
#       nd[i] <- nd[i] + sum(FDR<0.05)
      
#       # DEU-limma (F)
#       res_f <- paste0(DIR, "limma_diffSplice/", 
#                       simID, "/DEU_F_test.allGenes.tsv")
#       res <- read.table(res_f, header = TRUE)
#       FDR <- res$FDR
#       i <- 10
#       res.1 <- merge(map, res, by.x="gene", by.y="GeneID", all.x=TRUE)
#       status <- res.1[order(res.1$P.Value), "status"]
#       ranking[1:ngenes,i] <- status
#       nd[i] <- nd[i] + sum(FDR<0.05)
      
#       # Calculate cumulative fd
#       fd <- fd + apply(1-abs(ranking),2,cumsum)
#     }
#     fd_l[[s]] <- fd/nsim
#     fdmax_l[[s]] <- apply(fd_l[[s]],1,max)
#     nd_l[[s]] <- nd/nsim
#   }
# }

# i <- 1:1000
# fd.col <- c("green","green",
#             "blue", "blue", 
#             "violet", "black")
# fd.type <- c(1,2,1,2,1,1)

# pdf(paste0(fig, "fig_3.pdf"), height =7, width = 7)
# par(mfrow=c(2,2))

# s <- "balanced_3_75_3"
# plot(i, fdmax_l[[s]][i], type="n", xlab="Genes chosen", ylab="False discoveries", main="Balanced, n=3")
# lines(i,fd_l[[s]][i,1],col="green",lwd=2)
# lines(i,fd_l[[s]][i,2],col="green",lwd=2, lty=2)
# lines(i,fd_l[[s]][i,3],col="blue",lwd=2)
# lines(i,fd_l[[s]][i,4],col="blue",lwd=2, lty=2)
# lines(i,fd_l[[s]][i,5],col="violet",lwd=2)
# lines(i,fd_l[[s]][i,6],col="black",lwd=2)
# order <- 1:6
# legend("topleft", legend=methods[order], lwd=2, col=fd.col[order], lty=fd.type[order], bty="n", cex=1)

# s <- "unbalanced_3_75_3"
# plot(i, fdmax_l[[s]][i], type="n", xlab="Genes chosen", ylab="False discoveries", main="Unbalanced, n=3")
# lines(i,fd_l[[s]][i,1],col="green",lwd=2)
# lines(i,fd_l[[s]][i,2],col="green",lwd=2, lty=2)
# lines(i,fd_l[[s]][i,3],col="blue",lwd=2)
# lines(i,fd_l[[s]][i,4],col="blue",lwd=2, lty=2)
# lines(i,fd_l[[s]][i,5],col="violet",lwd=2)
# lines(i,fd_l[[s]][i,6],col="black",lwd=2)
# order <- 1:6

# s <- "balanced_5_75_3"
# plot(i, fdmax_l[[s]][i], type="n", xlab="Genes chosen", ylab="False discoveries", main="Balanced, n=5")
# lines(i,fd_l[[s]][i,1],col="green",lwd=2)
# lines(i,fd_l[[s]][i,2],col="green",lwd=2, lty=2)
# lines(i,fd_l[[s]][i,3],col="blue",lwd=2)
# lines(i,fd_l[[s]][i,4],col="blue",lwd=2, lty=2)
# lines(i,fd_l[[s]][i,5],col="violet",lwd=2)
# lines(i,fd_l[[s]][i,6],col="black",lwd=2)
# order <- 1:6

# s <- "unbalanced_5_75_3"
# plot(i, fdmax_l[[s]][i], type="n", xlab="Genes chosen", ylab="False discoveries", main="Unbalanced, n=5")
# lines(i,fd_l[[s]][i,1],col="green",lwd=2)
# lines(i,fd_l[[s]][i,2],col="green",lwd=2, lty=2)
# lines(i,fd_l[[s]][i,3],col="blue",lwd=2)
# lines(i,fd_l[[s]][i,4],col="blue",lwd=2, lty=2)
# lines(i,fd_l[[s]][i,5],col="violet",lwd=2)
# lines(i,fd_l[[s]][i,6],col="black",lwd=2)
# order <- 1:6

# dev.off()

# ################################################################################
# ###### Codes to produce results in Figure 4A and 4B
# ################################################################################
# setwd("/vast/projects/Spatial/tam/Differential_splicing/github/code/analysis/")
# source("edgeR_diffSpliceDGE_simple.R")
# source("plotJunc3_diffSpliceDGE.R")

# # DIR <- "/vast/projects/MM/tam/Differential_splicing/milevskiy_2023_GSE227748/"
# DIR <- "../../data/case_study/GSE227748/"
# p <- "LP-ML"
# # targetp <- paste0(DIR, "target/target.tsv")
# targetp <- paste0(DIR, "target/target.", p, ".tsv")
# featureCounts_o <- paste0(DIR, "featureCounts/")
# REF <- "../../annotation/"
# SJ <- paste0(REF, "gencode.vM32.primary_assembly.annotation.SJdatabase.tsv")
# SJ_database <- read.table(SJ, header=TRUE)

# message("Constructing design matrix for contrast ...")
# mat <- construct_model_matrix(targetp)

# message("Constructing a matrix for exon read counts with duplicated junction count ...")
# E <- process_E_count_withDupJ(featureCounts_o)

# message("Combining internal exon and junction read count matrix ...")  
# IE_J <- process_IE_J_count(featureCounts_o, SJ_database)

# message("Analyze libSize for DEU-edgeR ...")
# E_y <- norm_libSize(E$E_count,
#                     E$E_annot,
#                     mat$group, mat$design,
#                     "case_study", REF)

# message("Analyze libSize for DEJU-edgeR ...")
# IE_J_y <- norm_libSize(IE_J$IE_J_count,
#                        IE_J$IE_J_annot,
#                        mat$group, mat$design,
#                        "case_study", REF)

# message("Start DEU analysis for DEU-edgeR ...")
# DEU <- DEU_analysis(E_y, mat$design, p)

# message("Start DEJU analysis for DEJU-edgeR ...")
# DEJU <- DEU_analysis(IE_J_y, mat$design, p)

# ### Figure 4A
# message(paste0("Make VennDiagram for ", p, " ..."))
# sigDEUgenes <- list()
# # DEJU-edgeR-simes
# res_f <- paste0(DIR, "edgeR_diffSpliceDGE/", 
#                 p, "/DEJU_simes_test.sigGenes_0.05.tsv")
# res <- read.table(res_f, header = TRUE)
# i <- 1
# sigDEUgenes[[i]] <- res$GeneID

# # DEU-edgeR-simes
# res_f <- paste0(DIR, "edgeR_diffSpliceDGE/", 
#                 p, "/DEU_simes_test.sigGenes_0.05.tsv")
# res <- read.table(res_f, header = TRUE)
# i <- 2
# sigDEUgenes[[i]] <- res$GeneID

# # DEJU-edgeR-F
# res_f <- paste0(DIR, "edgeR_diffSpliceDGE/", 
#                 p, "/DEJU_F_test.sigGenes_0.05.tsv")
# res <- read.table(res_f, header = TRUE)
# i <- 3
# sigDEUgenes[[i]] <- res$GeneID

# # DEU-edgeR-F
# res_f <- paste0(DIR, "edgeR_diffSpliceDGE/", 
#                 p, "/DEU_F_test.sigGenes_0.05.tsv")
# res <- read.table(res_f, header = TRUE)
# i <- 4
# sigDEUgenes[[i]] <- res$GeneID

# l_simes <- list("DEJU-edgeR" = sigDEUgenes[[1]],
#                 "DEU-edgeR" = sigDEUgenes[[2]])

# l_F <- list("DEJU-edgeR" = sigDEUgenes[[3]],
#             "DEU-edgeR" = sigDEUgenes[[4]])

# colors <- c("blue", "grey")
# v_simes <- venn.diagram(l_simes, filename = NULL, disable.logging = TRUE, 
#              main=NULL, print.mode=c("raw", "percent"),
#              fill = colors, lwd=1)
# v_F <- venn.diagram(l_F, filename = NULL, disable.logging = TRUE, 
#                         main=NULL, print.mode=c("raw", "percent"),
#                         fill = colors, lwd=1)
# pdf(paste0(fig, "fig_4A.pdf"), height = 5, width = 2.5)
# grid.arrange(v_simes, v_F, nrow=2, ncol=1)
# dev.off()

# # Fgfr1 gene
# g <- "Fgfr1"
# message(paste0("Make plotJunc plot for ", g, " ..."))
# pdf(paste0(fig, "fig_4B_1.", g, ".pdf"), height = 8.5, width = 8.5)
# par(mfrow=c(2,1))
# plotJunc(DEU$sp, geneid=g, genecol="Symbol", annotation=IE_J$IE_annot)
# plotJunc(DEJU$sp, geneid=g, genecol="Symbol", annotation=IE_J$IE_annot)
# dev.off()

# message(paste0("Make sashimi plot for ", g, " ..."))
# geneSymbol <- "Fgfr1"
# geneID <- "ENSMUSG00000031565.19"
# chr <- "chr8"
# g_from <- 26002669
# g_to <- 26066734
# d <- 1000

# message("Load alignmnent track from bam files generated using make_bam_for_visualization.sh ...")
# samples <- c("LP_rep1", "LP_rep2", "LP_rep3", "ML_rep1", "ML_rep2", "ML_rep3")
# OUTPUT <- paste0(DIR, "Gviz/", geneSymbol, "_", geneID, "_", d, "/")
# PE_bam_files <- paste0(OUTPUT, samples, ".", geneSymbol, ".bam")

# alTrack <- list()

# for (idx in 1:3) {
#   s <- samples[idx]
#   alTrack[[s]] <- AlignmentsTrack(PE_bam_files[idx], 
#                                   isPaired = TRUE, 
#                                   fill.coverage="blue", 
#                                   col.sashimi="blue", 
#                                   fill="white", 
#                                   fontsize=8, cex.axis=1, col="blue")
# }

# for (idx in 4:6) {
#   s <- samples[idx]
#   alTrack[[s]] <- AlignmentsTrack(PE_bam_files[idx], 
#                                   isPaired = TRUE, 
#                                   fill.coverage="green", 
#                                   col.sashimi="green", 
#                                   fill="white", 
#                                   fontsize=8, cex.axis=1, col="green")
# }

# message("Load transcript annotation track ...")
# knownGenes <- UcscTrack(genome = "mm39", chromosome = chr, 
#                         track = "All GENCODE VM32", table="wgEncodeGencodeCompVM32", 
#                         from = g_from, to = g_to,
#                         trackType = "GeneRegionTrack", 
#                         rstarts = "exonStarts", rends = "exonEnds", 
#                         gene = "name", symbol = "name", 
#                         transcript = "name", strand = "strand", 
#                         fill = "black", name = "ALL GENCODE VM32", col="black")

# message("Load genome axis track ...")
# idxTrack <- IdeogramTrack(genome="mm39", chromosome=chr)
# axTrack <- GenomeAxisTrack()

# message("Combine all tracks ...")
# pdf(paste0(fig, "fig_4B_2.pdf"), height = 7, width = 4)
# plotTracks(c(idxTrack, axTrack, knownGenes, alTrack),
#            from = 26052007, to = 26059369, chromosome = chr, type = c("coverage", "sashimi"), 
#            sashimiNumbers=TRUE,
#            showTitle = FALSE,
#            sashimiScore=9,
#            lwd.sashimiMax=2,
#            sizes = c(0.1,0.2,0.3,rep(0.3, 6)),
#            lwd.title=2,
#            lwd.border.title=4,
#            background.title="white",
#            col.axis="black",
#            fontcolor="black",
#            col.line="black")
# dev.off()

# svg(paste0(fig, "fig_4B_2.svg"), height = 7, width = 4)
# plotTracks(c(idxTrack, axTrack, knownGenes, alTrack),
#            from = 26052007, to = 26059369, chromosome = chr, type = c("coverage", "sashimi"), 
#            sashimiNumbers=TRUE,
#            showTitle = FALSE,
#            sashimiScore=9,
#            lwd.sashimiMax=2,
#            sizes = c(0.1,0.2,0.3,rep(0.3, 6)),
#            lwd.title=2,
#            lwd.border.title=4,
#            background.title="white",
#            col.axis="black",
#            fontcolor="black",
#            col.line="black")
# dev.off()

# ### Figure 4C
# sigDEUgenes <- list()

# # DEJU-edgeR
# res_f <- paste0(DIR, "edgeR_diffSpliceDGE/", 
#                 p, "/DEJU_simes_test.sigGenes_0.05.tsv")
# res <- read.table(res_f, header = TRUE)
# i <- 1
# sigDEUgenes[[i]] <- res$GeneID

# # DEU-edgeR
# res_f <- paste0(DIR, "edgeR_diffSpliceDGE/", 
#                 p, "/DEU_simes_test.sigGenes_0.05.tsv")
# res <- read.table(res_f, header = TRUE)
# i <- 2
# sigDEUgenes[[i]] <- res$GeneID

# # DEJU-limma
# res_f <- paste0(DIR, "limma_diffSplice/", 
#                 p, "/DEJU_simes_test.sigGenes_0.05.tsv")
# res <- read.table(res_f, header = TRUE)
# i <- 3
# sigDEUgenes[[i]] <- res$GeneID

# # DEU-limma
# res_f <- paste0(DIR, "limma_diffSplice/", 
#                 p, "/DEU_simes_test.sigGenes_0.05.tsv")
# res <- read.table(res_f, header = TRUE)
# i <- 4
# sigDEUgenes[[i]] <- res$GeneID

# # DEXSeq
# res_f <- paste0(DIR, "DEXSeq/", 
#                 p, "/DEUs.geneWise.sigGenes.0.05.withGeneSymbol.tsv")
# res <- read.table(res_f, header = TRUE)
# i <- 5
# sigDEUgenes[[i]] <- res$groupID

# # JunctionSeq
# res_f <- paste0(DIR, "JunctionSeq/",
#                 p, "/sigGenes_0.05.genewiseResults.flt.withGeneSymbol.txt")
# res <- read.table(res_f, header = TRUE)
# i <- 6
# sigDEUgenes[[i]] <- res$geneID

# # DEJU-edgeR
# res_f <- paste0(DIR, "edgeR_diffSpliceDGE/", 
#                 p, "/DEJU_F_test.sigGenes_0.05.tsv")
# res <- read.table(res_f, header = TRUE)
# i <- 7
# sigDEUgenes[[i]] <- res$GeneID

# # DEU-edgeR
# res_f <- paste0(DIR, "edgeR_diffSpliceDGE/", 
#                 p, "/DEU_F_test.sigGenes_0.05.tsv")
# res <- read.table(res_f, header = TRUE)
# i <- 8
# sigDEUgenes[[i]] <- res$GeneID

# set_name <- c("DEJU-edgeR (simes)", "DEU-edgeR (simes)", 
#               "DEJU-limma", "DEU-limma",
#               "DEXSeq", "JunctionSeq",
#               "DEJU-edgeR (F)", "DEU-edgeR (F)")
# names(sigDEUgenes) <- set_name

# pdf(paste0(fig, "fig_4C.pdf"), height = 5, width = 6)
# upset(fromList(sigDEUgenes), order.by = "freq", nsets = length(sigDEUgenes), nintersects = 10)
# dev.off()

# ################################################################################
# ###### Codes to produce results in Tables 1
# ################################################################################
# simDir <- "../../data/simulation/"

# nlibs <- c(3, 5, 10)
# ls <- "balanced"
# fc <- 3
# rl <- 75
# methods <- c("DEJU-edgeR-simes",
#              "DEJU-limma-simes")
# methods.lb <-  gsub("-simes", "", methods)
# labels <- setNames(methods.lb, methods)
# sp <- c("ES", "MXE", "ASS", "RI")
# # num_row <- sum(length(sp)+length(nlibs)+length(methods))

# tab1 <- data.frame(matrix(ncol=5, nrow=0))
# colnames(tab1) <- c("Pattern", "n", "Method", "FDR", "Power")
# idx <- 1
# for (i in sp) {
#   for (n in nlibs) {
#     for (m in methods) {
#       s <- paste(ls, n, rl, fc, sep="_")
#       res <- read.csv(paste0(simDir, s, "/performance_analysis/benchmarking_result_sp.mean.csv"), header=TRUE)
#       res <- res[res$Method==m & res$Splicing_pattern==i,]
#       tab1[idx,] <- c(i, n, m, round(res$FDR.mean,3), round(res$Sensitivity.mean,3)) 
#       idx <- idx + 1
#     }
#   }
# }
# tab1$Method <- factor(tab1$Method, levels=methods)
# tab1$n <- factor(tab1$n, levels=c(3,5,10))
# tab1$Pattern <- factor(tab1$Pattern, levels=c("ES", "MXE", "ASS", "RI"))
# DEJU.edgeR <- tab1[tab1$Method=="DEJU-edgeR-simes", -3]
# DEJU.limma <- tab1[tab1$Method=="DEJU-limma-simes", -3]

# tab1 <- merge(DEJU.edgeR, DEJU.limma, by=c("Pattern", "n"), suffixes=c(".DEJU.edgeR", ".DEJU.limma"))
# tab1 <- tab1[order(tab1$Pattern, tab1$n), ]
# write.csv(tab1, paste0(fig, "table_1.csv"), row.names=FALSE)

# ################################################################################
# ###### Codes to produce results in Figure 4A and 4B
# ################################################################################
# simDir <- "../../data/simulation/"

# nlibs <- c(3, 5, 10)
# ls <- c("balanced", "unbalanced")
# fc <- 3
# rl <- 75
# tools <- c("edgeR_diffSpliceDGE", "limma_diffSplice", "DEXSeq", "JunctionSeq")
# idx <- 1:length(tools)

# extract_info <- function(f, p) {
#   lines <- readLines(f)
#   matched_lines <- lines[grep(p, lines)]
#   return(matched_lines)
# }

# tab2 <- data.frame(matrix(ncol=5, nrow=0))
# colnames(tab2) <- c("libSize", "n", "Method", "Walltime", "Memory")

# for (l in ls) {
#   for (n in nlibs) {
#     s <- paste(l, n, rl, fc, sep="_")
#     print(s)
#     logs <- paste0(simDir, s, "/log/", "03_", idx, "_DS_", tools, "_test.log")
#     for (i in idx) {
      
#       # Extract time-wall (TW)
#       TW <- extract_info(logs[i], "wall clock")
#       TW <- gsub("\tElapsed \\(wall clock\\) time \\(h:mm:ss or m:ss\\): ", "", TW)
#       # Extract peak memory usage (PMU)
#       PMU <- extract_info(logs[i], "Maximum resident set size")
#       PMU <- gsub("\tMaximum resident set size \\(kbytes\\): ", "", PMU)
#       df <- data.frame(libSize=l, n=n, Method=tools[i], Walltime=TW, Memory=PMU)
#       tab2 <- rbind(tab2, df)
#     }
#   }
# }

# time_convert <- function(TW) {
#   TW <- strsplit(TW, ":")[[1]]
#   if (length(TW)==2){
#     min <- as.numeric(TW[1])
#     sec <- as.numeric(TW[2])
#     TW <- round(min + sec/60,2)
#   } else if (length(TW)==3){
#     hr <- as.numeric(TW[1])
#     min <- as.numeric(TW[2])
#     sec <- as.numeric(TW[3])
#     TW <- round(hr*60 + min + sec/60,2)
#   }
#   return(TW)  
# }

# tab2$Walltime <- sapply(tab2$Walltime, time_convert)
# tab2$Memory <- round(as.numeric(tab2$Memory)/(1024*1024),2)
# tab2$mode <- ifelse(tab2$Method=="edgeR_diffSpliceDGE" | tab2$Method=="limma_diffSplice",
#                     "serial", NA)

# tab2[which(is.na(tab2$mode)), "mode"] <- rep(c("serial", "parallel"), length(which(is.na(tab2$mode)))/2)
# write.csv(tab2, paste0(fig, "table_2.csv"), row.names=FALSE)
