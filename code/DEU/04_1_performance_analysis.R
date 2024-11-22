OK <- requireNamespace("devtools", quietly = TRUE)
if (!OK) {
  stop("devtools package required but is not installed (or can't be loaded)")
}
library(BiocParallel)

# Getting parameters
args <- R.utils::commandArgs(asValues = TRUE)
print(args)

DIR <- as.character(args[['DIR']])
DEU_map <- as.character(args[['DEU_map']])
fdr_cutoff <- as.numeric(args[['fdr_cutoff']])
noOfSim <- as.integer(args[['noOfSim']])
seed <- as.integer(args[['seed']])
workers <- as.integer(args[['workers']])

# DIR <- "/vast/projects/Spatial/tam/Differential_splicing/github/data/simulation/unbalanced_3_75_3/"
# DEU_map <- "/vast/projects/Spatial/tam/Differential_splicing/github/data/simulation/customized_transcriptome/DEU_genes.info.tsv"
# noOfSim <- 20
# seed <- 2024
# workers <- 90
# fdr_cutoff <- 0.05

simID <- paste0("S", 1:noOfSim)

### 00. Check parameters
print(DIR)
print(DEU_map)
print(fdr_cutoff)
print(noOfSim)
print(simID)
print(seed)
print(workers)

methods <- c("DEU-edgeR-simes", "DEU-edgeR-F",
             "DEJU-edgeR-simes", "DEJU-edgeR-F", 
             "DEU-limma-simes", "DEU-limma-F",
             "DEJU-limma-simes", "DEJU-limma-F",
             "DEXSeq", "JunctionSeq")

AP_type <- c("ES", "MXE", "ASS", "RI")

### 01. Register BPPARAM
message('Register BPPARAM')
BPPARAM <- MulticoreParam(workers = workers, progressbar = TRUE, RNGseed = seed)
register(BPPARAM = BPPARAM)

### 02. Performance analysis
message('Start running JunctionSeq ...')
bplapply(simID, FUN = function(i){
  
  message(paste0('Processing: ', i, " ..."))
  analysis_o <- paste0(DIR, "performance_analysis/", i, "/")
  dir.create(analysis_o, showWarnings = F, recursive = T)
  
  print(analysis_o)
  
  ### 021. Load DEU results
  message('Loading DEU results ...')

  pred_DEUs <-list()
  
  # DEU-edgeR-simes
  k <- 1
  f <- paste0(DIR, "edgeR_diffSpliceDGE/", i, "/DEU_simes_test.sigGenes_0.05.tsv")
  df <- read.table(f, header = TRUE)
  pred_DEUs[[k]] <- df$GeneID

  # DEU-edgeR-F
  k <- 2
  f <- paste0(DIR, "edgeR_diffSpliceDGE/", i, "/DEU_F_test.sigGenes_0.05.tsv")
  df <- read.table(f, header = TRUE)
  pred_DEUs[[k]] <- df$GeneID
  
  # DEJU-edgeR-simes
  k <- 3
  f <- paste0(DIR, "edgeR_diffSpliceDGE/", i, "/DEJU_simes_test.sigGenes_0.05.tsv")
  df <- read.table(f, header = TRUE)
  pred_DEUs[[k]] <- df$GeneID
  
  # DEJU-edgeR-gene
  k <- 4
  f <- paste0(DIR, "edgeR_diffSpliceDGE/", i, "/DEJU_F_test.sigGenes_0.05.tsv")
  df <- read.table(f, header = TRUE)
  pred_DEUs[[k]] <- df$GeneID
  
  # DEU-limma-simes
  k <- 5
  f <- paste0(DIR, "limma_diffSplice/", i, "/DEU_simes_test.sigGenes_0.05.tsv")
  df <- read.table(f, header = TRUE)
  pred_DEUs[[k]] <- df$GeneID

  # DEU-limma-F
  k <- 6
  f <- paste0(DIR, "limma_diffSplice/", i, "/DEU_F_test.sigGenes_0.05.tsv")
  df <- read.table(f, header = TRUE)
  pred_DEUs[[k]] <- df$GeneID
  
  # DEJU-limma-simes
  k <- 7
  f <- paste0(DIR, "limma_diffSplice/", i, "/DEJU_simes_test.sigGenes_0.05.tsv")
  df <- read.table(f, header = TRUE)
  pred_DEUs[[k]] <- df$GeneID
  
  # DEJU-limma-F
  k <- 8
  f <- paste0(DIR, "limma_diffSplice/", i, "/DEJU_F_test.sigGenes_0.05.tsv")
  df <- read.table(f, header = TRUE)
  pred_DEUs[[k]] <- df$GeneID
  
  # DEXSeq
  k <- 9
  f <- paste0(DIR, "DEXSeq/", i, "/DEUs.geneWise.allGenes.tsv")
  df <- read.table(f, header = TRUE)
  pred_DEUs[[k]] <- df[df$padj.gene <= fdr_cutoff, "groupID"]
  # genes <- df[df$padj.gene <= fdr_cutoff, "groupID"]
  # pred_DEUs[[k]] <- genes[-grep("\\+", genes)]
  
  # JunctionSeq
  k <- 10
  f <- paste0(DIR, "JunctionSeq/", i, "/sigGenes.genewiseResults.txt.gz")
  df <- read.table(f, header = TRUE)
  pred_DEUs[[k]] <- df$geneID

  message('Load ground truth ...')
  
  simReads_o <- paste0(DIR, "simReads/", i, "/")
  DEU <- read.table(paste0(simReads_o, "DEU.tsv"), header=TRUE)
  ref_DEUs <- unique(gsub("\\.(skip_exon|splice_site_at_exon|retained_intron_).*", "", DEU$GeneID))
  map <- read.table(DEU_map, header=TRUE)
  
  ### Benchmarking results relative to splicing patterns
  res_sp <- data.frame(matrix(ncol=11, nrow=0))
  colnames(res_sp) <- c("simID", "method", "splicing_pattern", 
                        "TP", "FP", "FN", "FDR", "FNR", 
                        "Precision", "Sensitivity", "undefined")
  
  r <- 1
  for (t in AP_type){
    g <- ref_DEUs[ref_DEUs %in% map[map$type==t, "gene"]]
    
    ### Extract predicted DEU genes relative to splicing patterns 
    pred_DEUs_sp <- list()
    
    for (k in 1:10) {
      pred_DEUs_sp <- pred_DEUs[[k]][pred_DEUs[[k]] %in% map[map$type==t, "gene"]]
      # genes in both reference list and predicted list
      tp <- length(intersect(g, pred_DEUs_sp))
      # genes in predicted list but not in reference list
      fp <- length(setdiff(pred_DEUs_sp, g))
      # genes in reference list but not in predicted list
      fn <- length(setdiff(g, pred_DEUs_sp))
      # undefined genes (genes not simulated but detected)
      undefined <- length(setdiff(pred_DEUs_sp, map$gene))
      FDR <- fp / (tp + fp)
      FNR <- fn / (tp + fn)
      Pre <- tp / (tp + fp)
      Sen <- tp / (tp + fn)
      res_sp[r, ] <- c(i, methods[k], t, tp, fp, fn, FDR, FNR, Pre, Sen, undefined)
      r <- r + 1      
    }
  }
  
  ### General benchmarking results
  res_all <- data.frame(matrix(ncol=9, nrow=0))
  colnames(res_all) <- c("simID", "Method", 
                         "TP", "FP", "FN", "FDR", "FNR", 
                         "Precision", "Sensitivity")
  
  r <- 1
  for (k in 1:10) {
    # genes in both reference list and predicted list
    tp <- length(intersect(ref_DEUs, pred_DEUs[[k]]))
    # genes in predicted list but not in reference list
    fp <- length(setdiff(pred_DEUs[[k]], ref_DEUs))
    # genes in reference list but not in predicted list
    fn <- length(setdiff(ref_DEUs, pred_DEUs[[k]]))
    FDR <- fp / (tp + fp)
    FNR <- fn / (tp + fn)
    Pre <- tp / (tp + fp)
    Sen <- tp / (tp + fn)
    res_all[r, ] <- c(i, methods[k], tp, fp, fn, FDR, FNR, Pre, Sen)
    r <- r+1
  }
  
  ### Save final results
  write.csv(res_sp, 
              paste0(analysis_o, "/benchmarking_result_sp.csv"), 
              row.names = F)
  write.csv(res_all, 
              paste0(analysis_o, "/benchmarking_result.csv"), 
              row.names = F)
  
}, BPPARAM = BPPARAM)

message("Calculate mean and sd of FDR, FNR, Precision, and Sensitivity averaged for 20 simulations ...")

ave_res_sp <- data.frame(matrix(nrow = 0, ncol=16))
colnames(ave_res_sp) <- c("Method", "Splicing_pattern", 
                          "TP.mean", "TP.sd", 
                          "FP.mean", "FP.sd", 
                          "FN.mean", "FN.sd", 
                          "FDR.mean", "FDR.sd", 
                          "FNR.mean", "FNR.sd",
                          "Precision.mean", "Precision.sd",
                          "Sensitivity.mean", "Sensitivity.sd")

ave_res_all <- data.frame(matrix(nrow = 0, ncol=15))
colnames(ave_res_all) <- c("Method",
                           "TP.mean", "TP.sd", 
                           "FP.mean", "FP.sd", 
                           "FN.mean", "FN.sd", 
                           "FDR.mean", "FDR.sd", 
                           "FNR.mean", "FNR.sd",
                           "Precision.mean", "Precision.sd",
                           "Sensitivity.mean", "Sensitivity.sd")

res_sp <- list()
r_sp <- 1

res_all <- list()
r_all <- 1

message("Read bechmarking results from each simulation ...")
for (i in 1:noOfSim) {
  f <- paste0(DIR, "performance_analysis/S", i, "/benchmarking_result_sp.csv")
  res_sp[[i]] <- read.csv(f, header=TRUE)
  
  f <- paste0(DIR, "performance_analysis/S", i, "/benchmarking_result.csv")
  res_all[[i]] <- read.csv(f, header=TRUE)
}

for (m in methods) {
  
  message("Averaged for each splicing patterns ...")
  for (p in AP_type) {
    TP.sp <- rep(NaN, noOfSim)
    FP.sp <- rep(NaN, noOfSim)
    FN.sp <- rep(NaN, noOfSim)
    FDR.sp <- rep(NaN, noOfSim)
    FNR.sp <- rep(NaN, noOfSim)
    Pre.sp <- rep(NaN, noOfSim)
    Sen.sp <- rep(NaN, noOfSim)
    
    for (i in 1:noOfSim) {
      TP.sp[i] <- res_sp[[i]][(res_sp[[i]]$method==m & res_sp[[i]]$splicing_pattern==p), "TP"]
      FP.sp[i] <- res_sp[[i]][(res_sp[[i]]$method==m & res_sp[[i]]$splicing_pattern==p), "FP"]
      FN.sp[i] <- res_sp[[i]][(res_sp[[i]]$method==m & res_sp[[i]]$splicing_pattern==p), "FN"]
      FDR.sp[i] <- res_sp[[i]][(res_sp[[i]]$method==m & res_sp[[i]]$splicing_pattern==p), "FDR"]
      FNR.sp[i] <- res_sp[[i]][(res_sp[[i]]$method==m & res_sp[[i]]$splicing_pattern==p), "FNR"]
      Pre.sp[i] <- res_sp[[i]][(res_sp[[i]]$method==m & res_sp[[i]]$splicing_pattern==p), "Precision"]
      Sen.sp[i] <- res_sp[[i]][(res_sp[[i]]$method==m & res_sp[[i]]$splicing_pattern==p), "Sensitivity"]   
    }
    
    ave_res_sp[r_sp,] <- c(m, p,
                                  mean(TP.sp, na.rm=TRUE), sd(TP.sp, na.rm=TRUE), 
                                  mean(FP.sp, na.rm=TRUE), sd(FP.sp, na.rm=TRUE), 
                                  mean(FN.sp, na.rm=TRUE), sd(FN.sp, na.rm=TRUE),
                                  mean(FDR.sp, na.rm=TRUE), sd(FDR.sp, na.rm=TRUE), 
                                  mean(FNR.sp, na.rm=TRUE), sd(FNR.sp, na.rm=TRUE), 
                                  mean(Pre.sp, na.rm=TRUE), sd(Pre.sp, na.rm=TRUE), 
                                  mean(Sen.sp, na.rm=TRUE), sd(Sen.sp, na.rm=TRUE))
    r_sp <- r_sp + 1
  }
  
  TP.all <- rep(NaN, noOfSim)
  FP.all <- rep(NaN, noOfSim)
  FN.all <- rep(NaN, noOfSim)
  FDR.all <- rep(NaN, noOfSim)
  FNR.all <- rep(NaN, noOfSim)
  Pre.all <- rep(NaN, noOfSim)
  Sen.all <- rep(NaN, noOfSim)
  
  message("Averaged regardless of splicing patterns ...")
  for (i in 1:noOfSim) {
    TP.all[i] <- res_all[[i]][(res_all[[i]]$Method==m), "TP"]
    FP.all[i] <- res_all[[i]][(res_all[[i]]$Method==m), "FP"]
    FN.all[i] <- res_all[[i]][(res_all[[i]]$Method==m), "FN"]
    FDR.all[i] <- res_all[[i]][(res_all[[i]]$Method==m), "FDR"]
    FNR.all[i] <- res_all[[i]][(res_all[[i]]$Method==m), "FNR"]
    Pre.all[i] <- res_all[[i]][(res_all[[i]]$Method==m), "Precision"]
    Sen.all[i] <- res_all[[i]][(res_all[[i]]$Method==m), "Sensitivity"]
  }
  
  ave_res_all[r_all,] <- c(m,
                                  mean(TP.all, na.rm=TRUE), sd(TP.all, na.rm=TRUE),
                                  mean(FP.all, na.rm=TRUE), sd(FP.all, na.rm=TRUE),
                                  mean(FN.all, na.rm=TRUE), sd(FN.all, na.rm=TRUE),
                                  mean(FDR.all, na.rm=TRUE), sd(FDR.all, na.rm=TRUE),
                                  mean(FNR.all, na.rm=TRUE), sd(FNR.all, na.rm=TRUE),
                                  mean(Pre.all, na.rm=TRUE), sd(Pre.all, na.rm=TRUE),
                                  mean(Sen.all, na.rm=TRUE), sd(Sen.all, na.rm=TRUE))    
  r_all <- r_all + 1
}

write.csv(ave_res_sp, 
            paste0(DIR, "performance_analysis/benchmarking_result_sp.mean.csv"), 
            row.names = F)
write.csv(ave_res_all, 
            paste0(DIR, "performance_analysis/benchmarking_result.mean.csv"), 
            row.names = F)
