par.df <- data.frame(matrix(ncol=4))
colnames(par.df) <- c("scenario", "libs", "rlen", "fc")
i <- 1

dest.data <- "../../data/simulation/"

for (scenario in c("balanced", "unbalanced")) {
  for (libs in c(3, 5, 10)) {
    for (rlen in c(75)) {
      for (fc in c(3)) {
        par.df[i, ] <- c(scenario, libs, rlen, fc)
        i <- i + 1
        data.folder <- paste(scenario, libs, rlen, fc, sep="_")
        for (sub.folder in c("log",
                             "raw", "simReads",
                             "aligned_pass1", "SJ", "reindexed_genome", "aligned_pass2",
                             "featureCounts" , "edgeR_diffSpliceDGE", "limma_diffSplice",
                             "featureCounts_subread", "DEXSeq",
                             "QoRTs", "JunctionSeq")) {
          dir.create(paste0(dest.data, data.folder, "/", sub.folder), recursive = T, showWarnings = F)
        }
      }
    }
  }
}

write.table(par.df, "parameters_DEU.txt", row.names = F, sep="\t", quote = F)