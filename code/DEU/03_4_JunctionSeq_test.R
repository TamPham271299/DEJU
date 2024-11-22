OK <- requireNamespace("devtools", quietly = TRUE)
if (!OK) {
  stop("devtools package required but is not installed (or can't be loaded)")
}
library(JunctionSeq)

# Getting parameters
args <- R.utils::commandArgs(asValues = TRUE)
print(args)

DIR <- as.character(args[['DIR']])
targetp <- as.character(args[['target']])
ncores <- as.integer(args[['ncores']])

target <- read.table(targetp, header=TRUE, stringsAsFactors=FALSE)
design <- data.frame(condition = factor(target$group.ID))
i <- "S1"

### 0. Check parameter
message("Check valid parameters")
print(DIR)
print(targetp)
print(target)
print(design)
print(ncores)

message(paste0('Processing: ', i, " ..."))
count_files <- paste0(DIR, "STAR_output/QoRTs_fit_RL_75/", i, "/", target$sample.ID, "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz")
GFF <- paste0(DIR, "STAR_output/QoRTs_fit_RL_75/", i, "/withNovel.forJunctionSeq.gff.gz")
OUTPUT <- paste0(DIR, "STAR_output/JunctionSeq_fit_RL_75_test/", i, "/")
dir.create(OUTPUT, showWarnings = FALSE, recursive = TRUE)

print(count_files)
print(GFF)
print(OUTPUT)

### DEU analysis
message("readJunctionSeqCounts ...")
jscs <- readJunctionSeqCounts(countfiles = count_files,
                              samplenames = target$sample.ID,
                              design = design,
                              flat.gff.file = GFF)
message("Estimate size factors and dispersion ...")
jscs <- estimateJunctionSeqSizeFactors(jscs)
jscs <- estimateJunctionSeqDispersions(jscs, nCores = ncores)
jscs <- fitJunctionSeqDispersionFunction(jscs)
message("Test for DEU ...")
jscs <- testForDiffUsage(jscs, nCores = ncores, optimizeFilteringForAlpha = 0.05)
jscs <- estimateEffectSizes( jscs, nCores = ncores)
message("Write DEU results ...")
writeCompleteResults(jscs,
                     outfile.prefix = OUTPUT,
                     FDR.threshold = 0.05,
                     save.allGenes = TRUE,
                     save.sigGenes = TRUE,
                     save.bedTracks = FALSE)

