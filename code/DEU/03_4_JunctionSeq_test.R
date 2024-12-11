OK <- requireNamespace("devtools", quietly = TRUE)
if (!OK) {
  stop("devtools package required but is not installed (or can't be loaded)")
}
library(JunctionSeq)

read_target <- function(targetp, pair) {
  p <- strsplit(pair, "-")[[1]]
  target <- read.table(targetp, header=TRUE, stringsAsFactors=FALSE)
  design <- data.frame(condition = factor(target$group.ID, levels=c(p[2], p[1])))
  
  return(list("target"=target, "design"=design))
}

DEU_analysis <- function(count_files, target, design, GFF, ncores) {
  
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
  jscs <- testForDiffUsage(jscs, 
                           nCores = ncores, 
                           optimizeFilteringForAlpha = 0.05)
  jscs <- estimateEffectSizes(jscs, nCores = ncores)
}

# Getting parameters
args <- R.utils::commandArgs(asValues = TRUE)
print(args)

DIR <- as.character(args[['DIR']])
REF <- as.character(args[['REF']])
targetp <- as.character(args[['target']])
ncores <- as.integer(args[['ncores']])
pair <- as.character(args[['pair']])
i <- "S1"
target <- read_target(targetp, pair)

message("Check parameters ...")
print(DIR)
print(REF)
print(targetp)
print(target$target)
print(target$design)
print(pair)
print(ncores)

message(paste0('Processing: ', i, " ..."))
JunctionSeq_o <- paste0(DIR, "JunctionSeq/", i, "/")
dir.create(JunctionSeq_o, showWarnings = FALSE, recursive = TRUE)
QoRTs_o <- paste0(DIR, "QoRTs/", i, "/")

count_files <- paste0(QoRTs_o, target$target$sample.ID, 
                      "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz")
GFF <- paste0(QoRTs_o, "/withNovel.forJunctionSeq.gff.gz")

print(count_files)
print(GFF)

message("Run DEU analysis ...")
res <- DEU_analysis(count_files, target$target, target$design, GFF, ncores)
