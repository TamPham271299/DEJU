
# library(Rsubread)
# library(BiocParallel)

### Step 1:
readFasta <- function(fasta, metadatap) {
  
  ### scan fasta
  contigs <- scanFasta(fasta)
  metadata <- read.table(metadatap, header=TRUE)
  
  if (any(!contigs$Unique)) {
    
    ### Removing genes with duplicated isoforms if any
    dup <- unique(gsub("\\.(skip_exon|splice_site_at_exon|retained_intron_).*", "",
                       contigs[contigs$Unique==FALSE, ]$TranscriptID))
    dup_idx <- c()
    for (i in dup) {
      idx <- grep(i, contigs$TranscriptID)
      dup_idx <- c(dup_idx, idx)
    }
    
    contigs2 <- contigs[-dup_idx, ]
    ngenes <- nrow(contigs) / 2
    
    ### Updating gene map without genes with duplicated isoforms
    metadata <- metadata[!(metadata$gene %in% dup),]
    metadata <- "../../data/simulation/customized_transcriptome/DEU_genes.info.updated.tsv"
    if (!file.exists(metadata)) {
      write.table(metadata, 
                  "../../data/simulation/customized_transcriptome/DEU_genes.info.updated.tsv", 
                  sep="\t", quote=FALSE, row.names=FALSE)
    }

  } else {
    
    contigs2 <- contigs
    ngenes <- nrow(contigs) / 2
  }
  
  return(list('contigs' = contigs, 'contigs2' = contigs2, 'ngenes' = ngenes, 'metadata' = metadata))
}

### Step 2:
simulateBaseline <- function(ngenes, nlibs) {
  
  ### Use Zipf's law for generating baseline abundance for the first isoform
  k <- -0.9
  r <- 1:ngenes
  a <- 3e-5
  b <- 3e-10
  baselineprop <- exp(k*log(r)-r*a-(r^2)*b)
  baselineprop <- baselineprop/sum(baselineprop)
  mu0_tr <- baselineprop * 1e6
  
  ### shuffling transcript-level abundance
  mu0_tr <- sample(mu0_tr)
  mu0 <- c()
  
  ### simulating baseline abundance for the second isoform
  for (mu in mu0_tr){
    mu0 <- c(mu0, mu, mu * sample(seq(0.8, 1.2, 0.01), 1))
  }
  
  ### constructing transcript-level baseline abundance matrix across replicates
  mu0_mat <- matrix(rep(mu0, nlibs), ngenes * 2, nlibs)
  
  return(mu0_mat)
}

### Step 3:
selectDEUGenes <- function(metadata, contigs, numOfDEUs) {
  
  ### Collecting geneID from customized transcriptome
  genes <- gsub("\\.(skip_exon|splice_site_at_exon|retained_intron_).*", "", contigs$TranscriptID)
  
  ### Randomly partitioning a gene set into nonDEU and DEU genes featuring ES, MXE, ASS, RI
  splicing_pattern <- c("ES", "MXE", "ASS", "RI")
  DEU_G1 <- c()
  DEU_G2 <- c()
  
  for (p in splicing_pattern){
    idx <- which(genes %in% metadata[metadata$type==p, "gene"])
    idx_DEU_G1 <- idx[sample(seq(1, length(idx), 2), numOfDEUs[[p]])]
    idx_DEU_G2 <- idx_DEU_G1 + 1
    
    ### a list of indexes of DEU genes 
    DEU_G1 <- c(DEU_G1, idx_DEU_G1)
    DEU_G2 <- c(DEU_G2, idx_DEU_G2)
  }
  
  return(list('DEU_G1' = DEU_G1, 'DEU_G2' = DEU_G2))
}

### Step 4:
simulateDEUs <- function(DEU_G1, DEU_G2, nlibs, mu0_mat, FC){
  
  ### Multiply baseline expression value of the first isoform of the first group with a FC
  for (i in sort(DEU_G1)){
    for (j in 1:(nlibs/2)){
      mu0_mat[i, j] <- mu0_mat[i, j] * FC
    }
  }
  
  ### Multiply baseline expression value of the second isoform of the second group with a FC
  for (i in sort(DEU_G2)){
    for (j in (nlibs/2+1):nlibs){
      mu0_mat[i, j] <- mu0_mat[i, j] * FC
    }
  }
  
  return(mu0_mat)
}

### Step 5:
simulateBioVar <- function(mu0_mat, BCV0, ngenes, df.BCV, nlibs) {
  
  ### linearize baseline matrix into baseline vector
  mu0 <- c()
  for (col in 1:dim(mu0_mat)[2]){
    mu0 <- c(mu0, mu0_mat[, col])
  }
  
  ### Underlying dispersion trend
  BCV0 <- BCV0 + 1/sqrt(mu0 + 5) # dispersion trend
  
  ### Simulate true transcript-wise dispersions
  message('Simulate true transcript-wise dispersions ...')
  disp <- BCV0^2 * df.BCV/rep(rchisq(ngenes*2, df = df.BCV), nlibs) # tagwise dispersion
  shape <- 1/disp
  
  ### Simulate transcript abundance following gamma distribution
  message('Simulate transcript abundance following gamma distribution ...')
  scale <- mu0/shape
  mu <- matrix(rgamma(ngenes*2*nlibs, shape = shape, scale = scale), ngenes*2, nlibs) 
  
  return(mu)
}

### Step 6:
simulateTPM <- function(mu, contigs) {
  
  tpm <- mu / contigs$Length
  tpm <- 1e6 * tpm / matrix(colSums(tpm), dim(tpm)[1], dim(tpm)[2], byrow=TRUE)
  
  return(tpm)
}

### Step 7:
simulateFASTQ <- function(fasta, metadatap, nlibs, libSize, read.length, paired.end, simulate.sequencing.error, quality.reference, numOfDEUs, FC, BCV0, df.BCV, workers, seed, simReads, OUTPUT) {
  
  ### Setting up parallel computing
  if(is.null(seed)){
    BPPARAM <- MulticoreParam(workers = workers, progressbar = TRUE)
  } else{
    BPPARAM <- MulticoreParam(workers = workers, progressbar = TRUE, RNGseed = seed)
  }
  register(BPPARAM = BPPARAM)
  
  message('Preparing fasta file for simulation...')
  rFasta <- readFasta(fasta, metadatap)
  
  message('Use Zipf law for generating transcript-level baseline abundance ...')
  mu0_mat <- simulateBaseline(rFasta$ngenes, nlibs)
  
  message('Randomly subsetting genes featuring ES, MXE, ASS, RI...')
  selectIdx <- selectDEUGenes(rFasta$metadata, rFasta$contigs2, numOfDEUs)
  
  message('Simulate DEUs by multiplying transcript-level baseline with a fold-change ...')
  mu0_mat <- simulateDEUs(selectIdx$DEU_G1, selectIdx$DEU_G2, nlibs, mu0_mat, FC)
  
  message('Simulate variation of transcript abundance across replicates')
  mu <- simulateBioVar(mu0_mat, BCV0, rFasta$ngenes, df.BCV, nlibs)
  
  message('Simulate TPM value for each transcript ...')
  tpm <- simulateTPM(mu, rFasta$contigs2)
  
  message('Match back to the transcripts in the original FASTA file ...')
  TPM <- matrix(0, nrow(rFasta$contigs), nlibs) # create an empty matrix
  m <- match(rFasta$contigs2$TranscriptID, rFasta$contigs$TranscriptID)
  TPM[m,] <- tpm
  
  if (simulate.sequencing.error){
    if (is.null(quality.reference)) {
      # Getting quality reference
      quality.source <- ifelse(read.length %in% c(75, 100), 'Rsubread' , 'rfun')
      quality.reference <- list.files(system.file(package = quality.source,'qualf'),
                                      paste0('-',read.length,'bp'), full.names = TRUE)
    }
  } 
  
  ### Save metadata
  message('Save metadata: TPM values, a list of genuine DEU genes, read count output from simReads ...')
  
  # TPM values
  TPM_full <- cbind(rFasta$contigs$TranscriptID, as.data.frame(TPM))
  colnames(TPM_full)[1:(nlibs+1)] <- c("TranscriptID", paste0(rep("Rep", nlibs), 1:nlibs))
  write.table(TPM_full, file.path(simReads, "TPM_full.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
  
  # A list of genuine DEU genes
  DEU <- data.frame(GeneID=rbind(rFasta$contigs2[sort(selectIdx$DEU_G1),], rFasta$contigs2[sort(selectIdx$DEU_G2),])$TranscriptID, group=rep(c("Group_1", "Group_2"), each = sum(unlist(numOfDEUs))))
  write.table(DEU, file.path(simReads, "DEU.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

  message('Record parameters ...')
  print(sort(selectIdx$DEU_G1)[1:10])
  print(mu0_mat[1:14, ])
  print(mu[1:14, ])
  print(TPM[1:14, ])
  print(tpm[1:14, ])
  print(TPM_full[1:14, ])
  print(libSize)
  print(paired.end)
  print(read.length)
  print(simulate.sequencing.error)
  print(quality.reference)
  print(workers)
  print(seed)
  print(numOfDEUs)
  print(FC)
  print(BCV0)
  print(df.BCV)
  print(OUTPUT)
  print(DEU[1:10,])
  
  ### Run simReads
  message('Generate FASTQ files ...')
  out <- bplapply(seq_len(nlibs), FUN = function(i){
    simReads(transcript.file = fasta, expression.levels = TPM[,i],
             output.prefix = paste0(OUTPUT, "/Rep", i), library.size = libSize[i],
             paired.end = paired.end,
             read.length = read.length,
             simulate.sequencing.error = simulate.sequencing.error,
             quality.reference = quality.reference)
  }, BPPARAM = BPPARAM)
  
  ### Simulated read counts for each transcript from simReads
  mat <- do.call(cbind, lapply(out, function(x){x$NReads}))
  colnames(mat) <- colnames(TPM_full[, -1])
  out <- cbind(rFasta$contigs[,c("TranscriptID", "Length")], mat)
  rownames(out) <- NULL
  write.table(out, file.path(simReads, "simulation_3vs3_full_table_result.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
  
  message('Simulation of sequencing reads completed!')
}

### Run simulation
simulateExperiment <- function(project_p, 
                               fasta = "../../data/simulation/customized_transcriptome/gencode.vM32.custom.transcriptome.fa",
                               metadatap = "../../data/simulation/customized_transcriptome/DEU_genes.info.tsv",
                               simulation = 1, 
                               workers = 1, 
                               seed = NULL,
                               nlibs = 6, 
                               libSize = rep(50e6, nlibs), 
                               read.length = 75, 
                               paired.end = TRUE, 
                               simulate.sequencing.error = FALSE, 
                               quality.reference = NULL,
                               numOfDEUs = list(ES = 250, MXE = 250, ASS = 250, RI = 250), 
                               FC = 3, 
                               BCV0 = 0.05, 
                               df.BCV = 10){
  
  
  ### Generate FASTQ files for each simulation
  set.seed(seed)
  for (i in 1:simulation){
    
    OUTPUT <- paste0(project_p, "raw/S", i)
    simReads <- paste0(project_p, "simReads/S", i)
    
    dir.create(OUTPUT, showWarnings = FALSE, recursive = TRUE)
    dir.create(simReads, showWarnings = FALSE, recursive = TRUE)
    
    simulateFASTQ(fasta = fasta, 
                  metadatap = metadatap, 
                  nlibs = nlibs, 
                  libSize = libSize, 
                  read.length = read.length, 
                  paired.end = paired.end, 
                  simulate.sequencing.error = simulate.sequencing.error, 
                  quality.reference = quality.reference,
                  numOfDEUs = numOfDEUs, 
                  FC = FC, 
                  BCV0 = BCV0, 
                  df.BCV = df.BCV, 
                  workers = workers, 
                  seed = seed,
                  simReads = simReads,
                  OUTPUT = OUTPUT)
    
  }
}
