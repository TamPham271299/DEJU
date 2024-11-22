OK <- requireNamespace("devtools", quietly = TRUE)
if (!OK) {
  stop("devtools package required but is not installed (or can't be loaded)")
}

library(Rsubread)
library(BiocParallel)

# Getting parameters
args <- R.utils::commandArgs(asValues = TRUE)
print(args)

### Parameters for projects
project_p <- as.character(args[['path_to_project']])

### Parameters for customized transcriptome  
fasta <- as.character(args[['transcriptome_fasta']])
metadatap <- as.character(args[['metadata_of_transcriptome']])

### simReads parameters
libs.per.group <- as.integer(args[['libs']])
scenario <- args[['scenario']]
read.length <- as.integer(args[['rlen']])
paired.end <- as.logical(args[['pe']])
simulate.sequencing.error <- as.logical(args[['simulate_sequencing_error']])
quality.reference <- if (args[['quality_reference']] == "NA") {
  NULL
} else {
  as.character(args[['quality_reference']])
}
nlibs <- sum(rep(libs.per.group,2))

if (scenario == 'balanced') {
  libSize <- rep(50e6, nlibs)
}

if (scenario == 'unbalanced') {
  libSize <- rep(rep(c(25e6,100e6), length.out = libs.per.group), 2)
}

### Parameters to simulate DEUs
FC <- as.numeric(args[['fc']])
numOfDEUs.ES <- as.integer(args[['number_of_DEUs_for_ES_pattern']])
numOfDEUs.MXE <- as.integer(args[['number_of_DEUs_for_MXE_pattern']])
numOfDEUs.ASS <- as.integer(args[['number_of_DEUs_for_ASS_pattern']])
numOfDEUs.RI <- as.integer(args[['number_of_DEUs_for_RI_pattern']])
numOfDEUs <- list()
numOfDEUs[["ES"]] <- numOfDEUs.ES
numOfDEUs[["MXE"]] <- numOfDEUs.MXE
numOfDEUs[["ASS"]] <- numOfDEUs.ASS
numOfDEUs[["RI"]] <- numOfDEUs.RI

### General parameters
simulation <- as.integer(args[['number_of_simulation']])
seed <- as.integer(args[['seed']])
workers <- as.integer(args[['workers']])

source("simulation.R")

simulateExperiment(project_p = project_p,
                   fasta = fasta,
                   metadatap = metadatap,
                   nlibs = nlibs,
                   libSize = libSize,
                   read.length = read.length,
                   paired.end = paired.end,
                   simulate.sequencing.error = simulate.sequencing.error,
                   quality.reference = quality.reference,
                   numOfDEUs = numOfDEUs,
                   FC = FC,
                   simulation = simulation,
                   workers = workers,
                   seed = seed)
