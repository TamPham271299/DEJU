library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr)
library(gridExtra)
library(edgeR)
library(Gviz)
library(VennDiagram)
library(JunctionSeq)
library(DEXSeq)

options(ucscChromosomeNames=FALSE)

# setwd("/vast/projects/Spatial/tam/Differential_splicing/github/code/analysis/")
fig <- "../../figures/supp/"

dir.create(fig, recursive = T, showWarnings = F)

################################################################################
###### Codes to produce results in Figure S1
################################################################################
simDir <- "../../data/simulation/"

nlibs <- c(3, 5, 10)
libSize <- c("balanced", "unbalanced")
fc <- 3
rl <- 75
scenario <- c()
for (ls in libSize) {
  scenario <- c(scenario, paste(ls, nlibs, rl, fc, sep = "_"))
}

methods <- c("DEJU-edgeR-simes", 
             "DEU-edgeR-simes", 
             "DEJU-limma-simes", 
             "DEU-limma-simes", 
             "DEXSeq", 
             "JunctionSeq")
methods.lb <-  gsub("-simes", "", methods)
labels <- setNames(methods.lb, methods)

sp <- c("ES", "MXE", "ASS", "RI")
sp.lb <- c("ES", "MXE", "ASS", "IR")
sp.labels <- setNames(sp.lb, sp)

scenario.lb <- c("Balanced, n=3", "Balanced, n=5", "Balanced, n=10", 
                 "Unbalanced, n=3", "Unbalanced, n=5", "Unbalanced, n=10")
scenario.labels <- setNames(scenario.lb, scenario)

p1 <- list()
for (s in scenario[2]) {
  res <- read.csv(paste0(simDir, s, "/performance_analysis/benchmarking_result_sp.mean.csv"), header=TRUE)
  res <- res[grep(paste(methods, collapse="|"), res$Method),]
  
  df <- res[, c("Method", "Splicing_pattern", "TP.mean", "FP.mean", "FDR.mean")]
  df <- reshape2::melt(df, id.vars = c("Method", "FDR.mean", "Splicing_pattern"))
  df$variable <- factor(df$variable, levels=c("FP.mean", "TP.mean"))
  df$FDR.mean <- ifelse(df$variable == "FP.mean", df$FDR.mean, NA)
  df$Method <- factor(df$Method, levels=methods)
  df$Splicing_pattern <- factor(df$Splicing_pattern, levels=sp)
  
  p1[[s]] <- ggplot(df, aes(fill=variable, y=value, x=Method)) +
    geom_bar(stat="identity", color = "black") +
    facet_wrap(~ Splicing_pattern, ncol=4, labeller = labeller(.default = sp.labels)) +    
    scale_x_discrete(labels = labels) +
    labs(y="DEU genes", x=NULL, title = scenario.labels[[s]]) +
    ylim(0, 350) +
    scale_fill_manual(values = c("FP.mean" = "red", "TP.mean" = "grey"), 
                      labels=c("FP.mean" = "FP", "TP.mean" = "TP")) +
    theme_test() +
    theme(legend.position = c(0.99, 0.99),
          legend.justification = c(1, 1),
          legend.title = element_blank(),
          legend.key.size = unit(0.4, "cm"),
          legend.text = element_text(size = 9),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color="black"),
          axis.text.y = element_text(color="black"),
          plot.title = element_text(size = 11, hjust = 0.5, vjust = 2),
          axis.title.x = element_text(size=11),
          strip.background = element_blank()) +
    geom_text(aes(label=round(FDR.mean,2)), 
              position = position_stack(vjust = 1), 
              vjust = -0.5)
}

for (s in scenario[-c(1,2)]) {
  res <- read.csv(paste0(simDir, s, "/performance_analysis/benchmarking_result_sp.mean.csv"), header=TRUE)
  res <- res[grep(paste(methods, collapse="|"), res$Method),]
  
  df <- res[, c("Method", "Splicing_pattern", "TP.mean", "FP.mean", "FDR.mean")]
  df <- reshape2::melt(df, id.vars = c("Method", "Splicing_pattern", "FDR.mean"))
  df$variable <- factor(df$variable, levels=c("FP.mean", "TP.mean"))
  df$FDR.mean <- ifelse(df$variable == "FP.mean", df$FDR.mean, NA)
  df$Method <- factor(df$Method, levels=methods)
  df$Splicing_pattern <- factor(df$Splicing_pattern, levels=sp)
  
  p1[[s]] <- ggplot(df, aes(fill=variable, y=value, x=Method)) +
    geom_bar(stat="identity", color = "black") +
    facet_wrap(~ Splicing_pattern, ncol=4, labeller = labeller(.default = sp.labels)) +    
    scale_x_discrete(labels = labels) +
    labs(y="DEU genes", x=NULL, title = scenario.labels[[s]]) +
    ylim(0, 350) +
    scale_fill_manual(values = c("FP.mean" = "red", "TP.mean" = "grey")) +
    theme_test() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color="black"),
          axis.text.y = element_text(color="black"),
          plot.title = element_text(size = 11, hjust = 0.5, vjust = 2),
          axis.title.x = element_text(size=11),
          strip.text = element_blank()) +
    geom_text(aes(label=round(FDR.mean,2)), 
              position = position_stack(vjust = 1), 
              vjust = -0.5)  
}

pdf(paste0(fig, "fig_S1.pdf"), height = 12.5, width = 9)
grid.arrange(grobs=p1, 
             nrow=5, ncol=1, heights = c(2.2, rep(2, 4)))
dev.off()

################################################################################
###### Codes to produce results in Figure S2
################################################################################
simDir <- "../../data/simulation/"

methods <- c("DEJU-edgeR", "DEU-edgeR", "DEJU-limma", "DEU-limma", "DEXSeq", "JunctionSeq")
ngenes <- 5000
nmethods <- length(methods)
nsim <- 20

fd <- ranking <- matrix(0, nrow=ngenes, ncol=nmethods)
colnames(fd) <- colnames(ranking) <- methods
nd <- rep(0, nmethods)
names(nd) <- methods

### Reference gene list
map_f <- "../../data/simulation/customized_transcriptome/DEU_genes.info.tsv"
map <- read.table(map_f, header=TRUE)

### Scenarios
nlibs <- 10
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
    
    # BEGIN SIM
    for (simID in paste0("S", 1:nsim)) {
      print(paste0("SIM = ", simID))
      
      simReads_o <- paste0(DIR, "simReads/", simID, "/")
      DEU <- read.table(paste0(simReads_o, "DEU.tsv"), header=TRUE)
      ref.DEUs <- unique(gsub("\\.(skip_exon|splice_site_at_exon|retained_intron_).*", "", DEU$GeneID))
      map$status <- ifelse(map$gene %in% ref.DEUs, 1, 0)
      o <- order(map$status, decreasing = TRUE)
      map <- map[o, ]
      
      # DEJU-edgeR
      res_f <- paste0(DIR, "edgeR_diffSpliceDGE/", 
                      simID, "/DEJU_simes_test.allGenes.tsv")
      res <- read.table(res_f, header = TRUE)
      FDR <- res$FDR
      i <- 1
      res.1 <- merge(map, res, by.x="gene", by.y="GeneID", all.x=TRUE)
      status <- res.1[order(res.1$P.Value), "status"]
      ranking[1:ngenes,i] <- status
      nd[i] <- nd[i] + sum(FDR<0.05)
      
      # DEU-edgeR
      res_f <- paste0(DIR, "edgeR_diffSpliceDGE/", 
                      simID, "/DEU_simes_test.allGenes.tsv")
      res <- read.table(res_f, header = TRUE)
      FDR <- res$FDR
      i <- 2
      res.1 <- merge(map, res, by.x="gene", by.y="GeneID", all.x=TRUE)
      status <- res.1[order(res.1$P.Value), "status"]
      ranking[1:ngenes,i] <- status
      nd[i] <- nd[i] + sum(FDR<0.05)
      
      # DEJU-limma
      res_f <- paste0(DIR, "limma_diffSplice/", 
                      simID, "/DEJU_simes_test.allGenes.tsv")
      res <- read.table(res_f, header = TRUE)
      FDR <- res$FDR
      i <- 3
      res.1 <- merge(map, res, by.x="gene", by.y="GeneID", all.x=TRUE)
      status <- res.1[order(res.1$P.Value), "status"]
      ranking[1:ngenes,i] <- status
      nd[i] <- nd[i] + sum(FDR<0.05)
      
      # DEU-limma
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
      
      # Calculate cumulative fd
      fd <- fd + apply(1-abs(ranking),2,cumsum)
    }
    fd_l[[s]] <- fd/nsim
    fdmax_l[[s]] <- apply(fd_l[[s]],1,max)
    nd_l[[s]] <- nd/nsim
  }
}

i <- 1:1000
fd.col <- c("green","green","blue", "blue", "violet", "black")
fd.type <- c(1,2,1,2,1,1)

pdf(paste0(fig, "fig_S2.pdf"), height = 4, width = 8)
par(mfrow=c(1,2))

s <- "balanced_10_75_3"
plot(i, fdmax_l[[s]][i], ylim=c(1,700), type="n", xlab="Genes chosen", ylab="False discoveries", main="Balanced, n=10")
lines(i,fd_l[[s]][i,1],col="green",lwd=2)
lines(i,fd_l[[s]][i,2],col="green",lwd=2, lty=2)
lines(i,fd_l[[s]][i,3],col="blue",lwd=2)
lines(i,fd_l[[s]][i,4],col="blue",lwd=2, lty=2)
lines(i,fd_l[[s]][i,5],col="violet",lwd=2)
lines(i,fd_l[[s]][i,6],col="black",lwd=2)
order <- 1:6
legend("topleft", legend=methods[order], lwd=2, col=fd.col[order], lty=fd.type[order], bty="n", cex=1)

s <- "unbalanced_10_75_3"
plot(i, fdmax_l[[s]][i], ylim=c(1,700), type="n", xlab="Genes chosen", ylab="False discoveries", main="Unbalanced, n=10")
lines(i,fd_l[[s]][i,1],col="green",lwd=2)
lines(i,fd_l[[s]][i,2],col="green",lwd=2, lty=2)
lines(i,fd_l[[s]][i,3],col="blue",lwd=2)
lines(i,fd_l[[s]][i,4],col="blue",lwd=2, lty=2)
lines(i,fd_l[[s]][i,5],col="violet",lwd=2)
lines(i,fd_l[[s]][i,6],col="black",lwd=2)
order <- 1:6

dev.off()

################################################################################
###### Codes to produce results in Figure S3
################################################################################
simDir <- "../../data/simulation/"

nlibs <- c(3, 5, 10)
libSize <- c("balanced", "unbalanced")
fc <- 3
rl <- 75
scenario <- c()
for (ls in libSize) {
  scenario <- c(scenario, paste(ls, nlibs, rl, fc, sep = "_"))
}

methods <- c("DEJU-edgeR-simes", 
             "DEU-edgeR-simes", 
             "DEJU-limma-simes", 
             "DEU-limma-simes", 
             "DEXSeq", 
             "JunctionSeq")
methods.lb <-  gsub("-simes", "", methods)
labels <- setNames(methods.lb, methods)

sp <- c("ES", "MXE", "ASS", "RI")
sp.lb <- c("ES", "MXE", "ASS", "IR")
sp.labels <- setNames(sp.lb, sp)

scenario.lb <- c("Balanced, n=3", "Balanced, n=5", "Balanced, n=10", 
                 "Unbalanced, n=3", "Unbalanced, n=5", "Unbalanced, n=10")
scenario.labels <- setNames(scenario.lb, scenario)

p1 <- list()
for (s in scenario[1]) {
  res <- read.csv(paste0(simDir, s, "/performance_analysis/benchmarking_result_sp.mean.csv"), header=TRUE)
  res <- res[grep(paste(methods, collapse="|"), res$Method),]
  
  df <- res[, c("Method", "Splicing_pattern", "Sensitivity.mean", "Sensitivity.sd")]
  df$Method <- factor(df$Method, levels=methods)
  df$Splicing_pattern <- factor(df$Splicing_pattern, levels=sp)
  
  p1[[s]] <- ggplot(df) +
    geom_bar(aes(x=Method, y=Sensitivity.mean), stat="identity", fill="lightgray") +
    geom_errorbar(aes(x=Method, ymin=Sensitivity.mean-Sensitivity.sd,
                      ymax=Sensitivity.mean+Sensitivity.sd),
                  color="black", width=0.3, alpha=0.9, size=0.2) +
    facet_wrap(~ Splicing_pattern, ncol=4, labeller = labeller(.default = sp.labels)) +
    scale_x_discrete(labels = labels) +
    scale_y_continuous(limits = c(0, 1.2), breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1)) +
    labs(y="Power", x=NULL, title = scenario.labels[[s]]) +
    theme_test() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color="black"),
          axis.text.y = element_text(color="black"),
          plot.title = element_text(size = 10, hjust = 0.5, vjust = 2),
          strip.background = element_blank()) +
    geom_text(aes(x=Method, y=Sensitivity.mean, label=round(Sensitivity.mean,2)), vjust=-1.5, size=2.3)
}

for (s in scenario[-1]) {
  res <- read.csv(paste0(simDir, s, "/performance_analysis/benchmarking_result_sp.mean.csv"), header=TRUE)
  res <- res[grep(paste(methods, collapse="|"), res$Method),]
  
  df <- res[, c("Method", "Splicing_pattern", "Sensitivity.mean", "Sensitivity.sd")]
  df$Method <- factor(df$Method, levels=methods)
  df$Splicing_pattern <- factor(df$Splicing_pattern, levels=sp)
  
  p1[[s]] <- ggplot(df) +
    geom_bar(aes(x=Method, y=Sensitivity.mean), stat="identity", fill="lightgray") +
    geom_errorbar(aes(x=Method, ymin=Sensitivity.mean-Sensitivity.sd,
                      ymax=Sensitivity.mean+Sensitivity.sd),
                  color="black", width=0.3, alpha=0.9, size=0.2) +
    facet_wrap(~ Splicing_pattern, ncol=4, labeller = labeller(.default = sp.labels)) +
    scale_x_discrete(labels = labels) +
    scale_y_continuous(limits = c(0, 1.2), breaks = c(0, 0.5, 1)) +
    labs(y="Power", x=NULL, title = scenario.labels[[s]]) +
    theme_test() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color="black"),
          axis.text.y = element_text(color="black"),
          plot.title = element_text(size = 10, hjust = 0.5, vjust = 2),
          strip.text = element_blank()) +
    geom_text(aes(x=Method, y=Sensitivity.mean, label=round(Sensitivity.mean,2)), vjust=-1.5, size=2.3)

}

pdf(paste0(fig, "fig_S3.pdf"), height = 13, width = 9)
grid.arrange(grobs=p1, 
             nrow=6, ncol=1, heights = c(2.2, rep(2, 5)))
dev.off()

################################################################################
###### Codes to produce results in Figure S4
################################################################################
simDir <- "../../data/simulation/"

nlibs <- c(3, 5, 10)
libSize <- c("balanced", "unbalanced")
fc <- 3
rl <- 75
scenario <- c()
for (ls in libSize) {
  scenario <- c(scenario, paste(ls, nlibs, rl, fc, sep = "_"))
}

methods <- c("DEJU-edgeR-simes", 
             "DEU-edgeR-simes", 
             "DEJU-limma-simes", 
             "DEU-limma-simes", 
             "DEXSeq", 
             "JunctionSeq")
methods.lb <-  gsub("-simes", "", methods)
labels <- setNames(methods.lb, methods)

sp <- c("ES", "MXE", "ASS", "RI")
sp.lb <- c("ES", "MXE", "ASS", "IR")
sp.labels <- setNames(sp.lb, sp)

scenario.lb <- c("Balanced, n=3", "Balanced, n=5", "Balanced, n=10", 
                 "Unbalanced, n=3", "Unbalanced, n=5", "Unbalanced, n=10")
scenario.labels <- setNames(scenario.lb, scenario)

colors <- c("green", "green", "blue", "blue", "violet", "black")
ltys <- c(1,2,1,2,1,1)

df.full <- data.frame()
for (s in scenario) {
  res <- read.csv(paste0(simDir, s, "/performance_analysis/benchmarking_result_sp.mean.csv"), header=TRUE)
  res <- res[grep(paste(methods, collapse="|"), res$Method),]
  
  df <- res[, c("Method", "Splicing_pattern", "Sensitivity.mean", "Sensitivity.sd")]
  ls <- strsplit(s, "_")[[1]][1]
  n <- strsplit(s, "_")[[1]][2]
  df$libSize <- ls
  df$n <- paste0(n, "vs", n)
  df$Method <- factor(df$Method, levels=methods)
  df$Splicing_pattern <- factor(df$Splicing_pattern, levels=sp)
  df$n <- factor(df$n, levels=c("3vs3", "5vs5", "10vs10"))
  df.full <- rbind(df.full, df)
}

p1 <- ggplot(df.full, aes(x=n, y=Sensitivity.mean, group=Method, color=Method, linetype=Method)) +
  geom_line(linewidth=0.5) +
  geom_errorbar(aes(ymin=Sensitivity.mean-Sensitivity.sd, ymax=Sensitivity.mean+Sensitivity.sd, color=Method),
                width=0.1, size=0.4) +
  facet_grid(libSize ~ Splicing_pattern,
             labeller = labeller(Splicing_pattern = sp.labels,
                                 libSize = c("balanced" = "Balanced",
                                             "unbalanced" = "Unbalanced"))) +
  theme_test() +
  scale_color_manual(values = setNames(colors, methods), 
                     labels=labels) + 
  scale_linetype_manual(values = setNames(ltys, methods),
                        labels=labels) +
  labs(y="Power", x=NULL) +
  theme(legend.position = "bottom",        # Move legend to bottom
        legend.box = "horizontal",         # Arrange legend items horizontally
        plot.margin = margin(10, 10, 40, 10),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        plot.title = element_text(size = 10, hjust = 0.5, vjust = 2),
        strip.background = element_blank())

pdf(paste0(fig, "fig_S4.pdf"), height = 7, width = 9)
p1
dev.off()

################################################################################
###### Codes to produce results in Figure S5-S9 and Table S3
################################################################################
# setwd("/vast/projects/Spatial/tam/Differential_splicing/github/code/analysis/")
source("edgeR_diffSpliceDGE_simple.R")
source("plotJunc3_1.R")

# DIR <- "/vast/projects/MM/tam/Differential_splicing/milevskiy_2023_GSE227748/"
DIR <- "../../data/case_study/GSE227748/"
p <- "LP-ML"
targetp <- paste0(DIR, "target/target.", p, ".tsv")
featureCounts_o <- paste0(DIR, "featureCounts/")
REF <- "../../annotation/"
SJ <- paste0(REF, "gencode.vM32.primary_assembly.annotation.SJdatabase.tsv")
SJ_database <- read.table(SJ, header=TRUE)
g <- gsub("-", "|", p)

message("Constructing design matrix for contrast ...")
mat <- construct_model_matrix(targetp)

message("Constructing a matrix for exon read counts with duplicated junction count ...")
E <- process_E_count_withDupJ(featureCounts_o)
E$E_count <- E$E_count[, grep(g, colnames(E$E_count))]

message("Combining internal exon and junction read count matrix ...")  
IE_J <- process_IE_J_count(featureCounts_o, SJ_database)
IE_J$IE_J_count <- IE_J$IE_J_count[, grep(g, colnames(IE_J$IE_J_count))]
IE_J$J_count <- IE_J$J_count[, grep(g, colnames(IE_J$J_count))]
IE_J$IE_count <- IE_J$IE_count[, grep(g, colnames(IE_J$IE_count))]

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

message("Analyze libSize for junction reads ...")
J_y <- norm_libSize(IE_J$J_count,
                    IE_J$J_annot,
                    mat$group, mat$design,
                    "case_study", REF)

message("Analyze libSize for internal exon reads ...")
IE_y <- norm_libSize(IE_J$IE_count,
                     IE_J$IE_annot,
                     mat$group, mat$design,
                     "case_study", REF)

message("Start DEU analysis for DEU-edgeR ...")
DEU <- DEU_analysis(E_y, mat$design, p)

message("Start DEJU analysis for DEJU-edgeR ...")
DEJU <- DEU_analysis(IE_J_y, mat$design, p)

### Figure S5
message(paste0("Generate BCV, QL, MDS plot for DEU-edgeR and DEJU-edgeR for ", p, " ..."))
points <- c(1,2)
colors <- c("blue", "green")

# DEU-edgeR
png(file.path(fig, "fig_S5_1.png"), height = 2.5*300, width = 8*300, res=300)
par(mfrow=c(1,3))
plotMDS(E_y, col=colors[group])
plotBCV(E_y)
plotQLDisp(DEU$fit)
dev.off()

# DEJU-edgeR
png(file.path(fig, "fig_S5_2.png"), height = 2.5*300, width = 8*300, res=300)
par(mfrow=c(1,3))
plotMDS(IE_J_y, col=colors[group])
plotBCV(IE_J_y)
plotQLDisp(DEJU$fit)
dev.off()

### Table S3
message("Calculate libSize for internal exon, junction reads, DEU-edgeR, and DEJU-edgeR ...")
libSize <- data.frame(Sample=row.names(J_y$samples), 
                      Group=J_y$samples$group, 
                      internal_exon=IE_y$samples$lib.size,
                      Junction=J_y$samples$lib.size,
                      `DEU-edgeR`=E_y$samples$lib.size,
                      `DEJU-edgeR`=IE_J_y$samples$lib.size)
colnames(libSize) <- c("Sample", "Group", "Internal exon reads", "Junction reads", "DEU-edgeR", "DEJU-edgeR")
write.csv(libSize, paste0(fig, "table_S3.csv"), row.names=FALSE)

### Figure S6
# Numb gene
g <- "Numb"
message(paste0("Make plotJunc plot for ", g, " ..."))
pdf(paste0(fig, "fig_S6_1.pdf"), height = 8.5, width = 8.5)
par(mfrow=c(2,1))
plotJunc(DEU$sp, geneid=g, genecol="Symbol", annotation=IE_J$IE_annot)
plotJunc(DEJU$sp, geneid=g, genecol="Symbol", annotation=IE_J$IE_annot)
dev.off()

message(paste0("Make sashimi plot for ", g, " ..."))
g <- "Numb"
geneID <- "ENSMUSG00000021224.16"
chr <- "chr12"
g_from <- 83839807
g_to <- 83969708
d <- 1000

message("Load alignmnent track from bam files generated using make_bam_for_visualization.sh ...")
samples <- c("LP_rep1", "LP_rep2", "LP_rep3", "ML_rep1", "ML_rep2", "ML_rep3")
OUTPUT <- paste0(DIR, "Gviz/", g, "_", geneID, "_", d, "/")
PE_bam_files <- paste0(OUTPUT, samples, ".", g, ".bam")

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

pdf(paste0(fig, "fig_S6_2.exon8-14.", g, ".pdf"), width = 4, height = 7)
plotTracks(c(idxTrack, axTrack, knownGenes, alTrack),
           from = 83840892, to = 83858214, chromosome = chr, type = c("coverage", "sashimi"), 
           sashimiNumbers=TRUE,
           showTitle = FALSE,
           sashimiScore=20,
           lwd.sashimiMax=2,
           sizes = c(0.1,0.2,0.3,rep(0.3, 6)),
           lwd.title=2,
           lwd.border.title=4,
           background.title="white",
           col.axis="black",
           fontcolor="black",
           col.line="black")
dev.off()

# svg(paste0(fig, "fig_S6_2_exon8-14.", g, ".svg"), width = 4, height = 7)
# plotTracks(c(idxTrack, axTrack, knownGenes, alTrack),
#            from = 83840892, to = 83858214, chromosome = chr, type = c("coverage", "sashimi"), 
#            sashimiNumbers=TRUE,
#            showTitle = FALSE,
#            sashimiScore=20,
#            lwd.sashimiMax=2,
#            sizes = c(0.1,0.2,0.3,rep(0.3, 6)),
#            lwd.title=2,
#            lwd.border.title=4,
#            background.title="white",
#            col.axis="black",
#            fontcolor="black",
#            col.line="black")
# dev.off()

pdf(paste0(fig, "fig_S6_2_exon1-6.", g, ".pdf"), width = 4, height = 7)
plotTracks(c(idxTrack, axTrack, knownGenes, alTrack),
           from = 83878280, to = 83968734, chromosome = chr, type = c("coverage", "sashimi"), 
           sashimiNumbers=TRUE,
           showTitle = FALSE,
           sashimiScore=5,
           lwd.sashimiMax=2,
           sizes = c(0.1,0.2,0.3,rep(0.3, 6)),
           lwd.title=2,
           lwd.border.title=4,
           background.title="white",
           col.axis="black",
           fontcolor="black",
           col.line="black", ylim=c(0,100))
dev.off()

# svg(paste0(fig, "fig_S6_2_exon1-6.", g, ".svg"), width = 4, height = 7)
# plotTracks(c(idxTrack, axTrack, knownGenes, alTrack),
#            from = 83878280, to = 83968734, chromosome = chr, type = c("coverage", "sashimi"), 
#            sashimiNumbers=TRUE,
#            showTitle = FALSE,
#            sashimiScore=5,
#            lwd.sashimiMax=2,
#            sizes = c(0.1,0.2,0.3,rep(0.3, 6)),
#            lwd.title=2,
#            lwd.border.title=4,
#            background.title="white",
#            col.axis="black",
#            fontcolor="black",
#            col.line="black", ylim=c(0,100))
# dev.off()

# Dusp16 gene
g <- "Dusp16"
message(paste0("Make plotJunc plot for ", g, " ..."))
pdf(paste0(fig, "fig_S6_1.", g, ".pdf"), height = 8.5, width = 8.5)
par(mfrow=c(2,1))
plotJunc(DEU$sp, geneid=g, genecol="Symbol", annotation=IE_J$IE_annot)
plotJunc(DEJU$sp, geneid=g, genecol="Symbol", annotation=IE_J$IE_annot)
dev.off()

message(paste0("Make sashimi plot for ", g, " ..."))
g <- "Dusp16"
geneID <- "ENSMUSG00000030203.18"
chr <- "chr6"
g_from <- 134691430
g_to <- 134770588
d <- 1000

message("Load alignmnent track from bam files generated using make_bam_for_visualization.sh ...")
samples <- c("LP_rep1", "LP_rep2", "LP_rep3", "ML_rep1", "ML_rep2", "ML_rep3")
OUTPUT <- paste0(DIR, "Gviz/", g, "_", geneID, "_", d, "/")
PE_bam_files <- paste0(OUTPUT, samples, ".", g, ".bam")

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

pdf(paste0(fig, "fig_S6_2.", g, ".pdf"), width = 4, height = 7)
plotTracks(c(idxTrack, axTrack, knownGenes, alTrack),
           from = 134737452, to = 134789796, chromosome = chr, type = c("coverage", "sashimi"), 
           sashimiNumbers=TRUE,
           showTitle = FALSE,
           sashimiScore=15,
           lwd.sashimiMax=2,
           sizes = c(0.1,0.2,0.3,rep(0.3, 6)),
           lwd.title=2,
           lwd.border.title=4,
           background.title="white",
           col.axis="black",
           fontcolor="black",
           col.line="black",
           ylim=c(0,200))
dev.off()

# svg(paste0(fig, "fig_S6_2.", g, ".svg"), width = 4, height = 7)
# plotTracks(c(idxTrack, axTrack, knownGenes, alTrack),
#            from = 134737452, to = 134789796, chromosome = chr, type = c("coverage", "sashimi"), 
#            sashimiNumbers=TRUE,
#            showTitle = FALSE,
#            sashimiScore=15,
#            lwd.sashimiMax=2,
#            sizes = c(0.1,0.2,0.3,rep(0.3, 6)),
#            lwd.title=2,
#            lwd.border.title=4,
#            background.title="white",
#            col.axis="black",
#            fontcolor="black",
#            col.line="black",
#            ylim=c(0,200))
# dev.off()

### Figure S9 (JunctionSeq plot)
jscs <- readRDS(paste0(DIR, "JunctionSeq/", p, "/jscs.rds"))

# Dusp16 gene
g <- "Dusp16"
geneID <- "ENSMUSG00000030203.18"
buildAllPlotsForGene(geneID, jscs, outfile.prefix = paste0(fig, "/fig_S9.", g, "."))

### Figure S8 (DEXSeq plots)
dxr <- readRDS(paste0(DIR, "DEXSeq/", p, "/dxr_exon.rds"))

# Dusp16 gene
g <- "Dusp16"
geneID <- "ENSMUSG00000030203.18"
pdf(paste0(fig, "fig_S8.", g, ".pdf"), width = 12, height = 7)
plotDEXSeq(dxr, geneID, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2)
dev.off()

# Fgfr1 gene
g <- "Fgfr1"
geneID <- "ENSMUSG00000031565.19"
pdf(paste0(fig, "fig_S8.", g, ".pdf"), width = 12, height = 7)
plotDEXSeq(dxr, geneID, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2)
dev.off()

# Numb gene
g <- "Numb"
geneID <- "ENSMUSG00000021224.16"
pdf(paste0(fig, "fig_S8.", g, ".pdf"), width = 12, height = 7)
plotDEXSeq(dxr, geneID, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2)
dev.off()
