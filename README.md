# Differential Exon-Junction Usage (DEJU)

This is the repository for the code used to perform analysis and generate figures for the following paper titled "Incorporating exon-exon junction reads enhances differential splicing detection". 

### Citation

If you are using code or pipelines from this repository, please consider citing our associated article:

Pham, M.T., Milevskiy, M.J.G., Visvader, J.E. et al. Incorporating exon–exon junction reads enhances differential splicing detection. BMC Bioinformatics 26, 193 (2025). [https://doi.org/10.1186/s12859-025-06210-4](https://doi.org/10.1186/s12859-025-06210-4)

### Introduction

The paper introduced a DEJU analysis workflow implementing a **STAR-Rsubread-edgeR-limma framework** to identify DEU genes indicative of differential splicing between experimental conditions in RNA-seq data. Here is a schematic presentation of our proposed DEJU workflow.

<p align="center">
  <img src="figures/fig_1_v2.png" alt="DEJU_workflow" width="400"/>
</p>

The paper also benchmarked the DEJU analysis methods implemented in edgeR and limma (DEJU-edgeR, DEJU-limma) against the existing DEU analysis methods (DEU-edgeR, DEU-limma) and other popular R-based DEU/DEJU tools (DEXSeq, JunctionSeq) based on simulated RNA-seq datasets. We also performed DEU analysis on RNA-sequencing experiments of [NCBI GEO database (GSE227748)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE227748) using our proposed method and other DEU analysis pipelines benchmarked in the study.

### Results and conclusion

By incorporating exon-exon junction reads, our DEJU methods demonstrates higher performance over other benchmarked methods in FDR control, statistical power, computational efficiency (turn-around time + memory usage), and flexibility in detecting a broad range of AS events, notably alternative splice sites and intron retention, making it the most suitable candidate for DEJU analysis in RNA-seq data.
In practical applications, our DEJU method effectively handles splicing alterations involving multiple known and novel transcripts, while supporting complex experimental designs. The effectiveness of the DEJU approach increases with sample size, making it particularly advantageous in large-scale studies (e.g. TCGA). However, it is not recommended for non-model organisms with the incomplete reference genome and does not provide transcript-level abundance estimates.

### Repository Structure

```plaintext
DEJU/
├── annotation/
│   ├── gene_info_GRCm39_M32_ensembl109.tsv    # Gene information
│   └── duplicated_sequences.tsv    # duplicate sequences that should be ignored for the simulation process
├── code/
│   ├── annotation_dl    # Script to prepare reference genome datasets of mm39 (GTF, FASTA, exon-junction SAF...) 
│   ├── alignment      # Script to map reads to reference genome mm39
│   ├── DEU            # Script to quantify exon-junction counts and detect DEU genes by 6 benchmarked methods
│   ├── simulation    # Main script to generate customized transcriptome, generate and analyze simulation data
│   ├── fastq_dl      # Script to download FASTQ files with GEO acession number
│   ├── case_study    # Main script to analyze case-study data
│   └── analysis      # Main script to perform overall analysis, create figures for the paper
├── tools/
│   └── Subread_to_DEXSeq/    # External scripts to generate exon counts for DEXSeq
├── figures/
│   └── fig_S1.png    # DEJU analysis worflow
└── README.md                  # Main documentation for the repository
```

### DEJU workflow tutorial

In this example, we provide basic steps of our DEJU-edgeR method to detect DEU genes from paired-end RNA-seq data for 2 groups with 2 biological replicates (e.g `sample1_G1`, `sample2_G1`, `sample1_G2`, `sample2_G2`).\
We can also apply this pipeline for single-end RNA-seq data.\
More replicates we have, higher sensitivity and specificity we get for the differential splicing detection result.

#### 0. Reference genome

First, we download genomic annotation `hg38.genome.gtf` and genomic sequence `hg38.genome.fasta` of the reference genome (e.g., from Gencode, UCSC database).\
**Note:** In this tutorial, we use genomic human annotation hg38 from the Gencode database.\
We can download latest version of FASTA and GTF of hg38 from [Gencode](https://www.gencodegenes.org/human/).\
To generate flattened and merged exon annotation, please visit [code/annotation_dl/GTF2SAF.R](code/annotation_dl/GTF2SAF.R) for more details.\
To generate junction database, please visit [code/annotation_dl/GTF2SJdatabase.R](code/annotation_dl/GTF2SJdatabase.R) for more details.

**Input:** `hg38.genome.gtf`, `hg38.genome.fasta`

**Output:** `hg38.flat_exon.saf`, `hg38.SJdatabase.tsv`

```bash
# Install Rsubread if not yet installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Rsubread")

library(Rsubread)

# convert GTF to SAF
SAF <- flattenGTF("hg38.genome.gtf", 
                  GTF.featureType="exon", 
                  GTF.attrType="gene_id", 
                  method="merge")

write.table(SAF, "hg38.flat_exon.saf", quote=F, row.names=F, sep="\t")

# Generate junction database for the reference hg38 genome
# install the following packages if not yet installed
install.packages("dplyr")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("rtracklayer", "GenomicRanges"))

library(dplyr)
library(rtracklayer)
library(GenomicRanges)

# Getting parameters
message('Extracting exon information from GTF ...')
GTF <- import("hg38.genome.gtf")
exon_dt <- subset(GTF, type == "exon")
exon_info <- data.frame(chr = seqnames(exon_dt),
                        start = start(exon_dt),
                        end = end(exon_dt),
                        strand = strand(exon_dt),
                        geneID = mcols(exon_dt)$gene_id,
                        transcriptID = mcols(exon_dt)$transcript_id)

message('Converting exon annotation SAF to SJ database ...')

SJ <- exon_info %>%
  group_by(transcriptID) %>%
  arrange(chr, start, end) %>%
  mutate(start = lead(start)) %>%
  slice(-n()) %>%
  ungroup()

SJ <- SJ[c("geneID", "chr", "end", "start", "strand")]
colnames(SJ)[3:4] <- c("start", "end")
SJ$juncID <- paste(SJ$chr, SJ$start, SJ$end, sep="_")

SJ.1 <- SJ %>%
  group_by(geneID) %>%
  distinct(juncID, .keep_all = TRUE) %>%
  ungroup() %>%
  count(juncID, name = "freq") %>%
  arrange(geneID, chr, start, end)

write.table(SJ.1, "hg38.SJdatabase.tsv", quote=F, row.names=F, sep="\t")
```
#### 1. Per-sample FASTQ pre-processing and exon-junction read mapping with 2-pass mode (STAR)

**Input:** `hg38.genome.gtf`, `hg38.genome.fasta`,\
`sample1_G1_R1.fastq.gz`, `sample1_G1_R2.fastq.gz`, `sample2_G1_R1.fastq.gz`,  `sample2_G1_R2.fastq.gz`, `sample1_G2_R1.fastq.gz`, `sample1_G2_R2.fastq.gz`, `sample2_G2_R1.fastq.gz`, `sample2_G2_R2.fastq.gz` inside `raw` folder

**Final Output:**\
`sample1_G1_R1.fastq.gz`, `sample1_G1_R2.fastq.gz`, `sample2_G1_R1.fastq.gz`,  `sample2_G1_R2.fastq.gz`, `sample1_G2_R1.fastq.gz`, `sample1_G2_R2.fastq.gz`, `sample2_G2_R1.fastq.gz`, `sample2_G2_R2.fastq.gz` inside `trimmed` folder
`sample1_G1.Aligned.sortedByCoord.out.bam`, `sample2_G1.Aligned.sortedByCoord.out.bam`, `sample1_G2.Aligned.sortedByCoord.out.bam`, `sample2_G2.Aligned.sortedByCoord.out.bam` inside `aligned_pass2` folder

We recommend run STAR with 2-pass mapping with re-generated genome + manual SJ filtering process like below to get the high mapping quality of exon-exon junction reads.
Otherwise, we can freely do mapping with any splice-aware aligners, e.g. STAR, HISAT2, Rsubread::subjuc

An example of STAR with 2-pass mapping mode, followed by FASTQC/trim_galore used in this study is descibed below. For more details, please visit [`code/alignment/`](code/alignment/).
This mapping step is referred to [https://www.reneshbedre.com/blog/star-aligner-twopass-mode.html](https://www.reneshbedre.com/blog/star-aligner-twopass-mode.html) and [STAR documentation](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf)

```bash
# Install FASTQC, trim_galore, STAR if not yet installed

# For example, two FASTQ files of paired end reads of sample1_G1 are stored in raw folder with the name "sample1_G1_R1.fastq.gz" and "sample1_G1_R2.fastq.gz"
# R1 and R2 stands for forward and reverse FASTQ reads, respectively

# Step 1: Run FASTQC + trim reads and adapters (if has not been trimmed)

fastqc --nogroup -t 16 --dir raw -o raw raw/sample1_G1_R1.fastq.gz raw/sample1_G1_R2.fastq.gz

trim_galore --cores 4 -q 30 --length 20 --paired --fastqc_args "--nogroup" --output_dir trimmed raw/sample1_G1_R1.fastq.gz raw/sample1_G1_R2.fastq.gz

# Then, we can rename trimmed fastq files (often with "*_val_*.fq.gz") into "sample1_G1_R1.fastq.gz" and "sample1_G1_R2.fastq.gz" (if paired-end reads)
# for example, run this:
for i in `find trimmed -name "*_val_*.fq.gz"`; do mv "$i" "${i/_val_?.fq.gz/.fastq.gz}"; done
# Then trimmed FASTQ files will be stored in trimmed folder with the name "sample1_G1_R1.fastq.gz" and "sample1_G1_R2.fastq.gz"

# Index the reference human genome (e.g. hg38) if not yet indexed
# Note: Only run once for 1 reference genome
sjdbOverhang=99 # optimal SJ overhang if max read length is 100
mkdir hg38_idx_genome
STAR --runThreadN 16 \
          --runMode genomeGenerate \
          --genomeDir hg38_idx_genome \
          --genomeFastaFiles hg38.genome.fasta \
          --sjdbGTFfile hg38.genome.gtf \           
          --sjdbOverhang $sjdbOverhang

# Step 2: Per-sample 1st-mapping pass
mkdir aligned_pass1
STAR --genomeDir $genomeDir \
                --readFilesIn trimmed/sample1_G1_R1.fastq.gz trimmed/sample1_G1_R2.fastq.gz \
                --readFilesCommand zcat \
                --outFileNamePrefix  aligned_pass1/sample1_G1. \
                --outSAMtype BAM SortedByCoordinate \
                --runThreadN 16

# Step 3: After running alignment for 4 samples, then
# SJ filtering for 4 samples (Manually filtering junctions  supported by less than 3 UMRs to reduce false positive junctions)
mkdir SJ
cp aligned_pass1/*/*SJ.out.tab SJ
cat $DIR/SJ/*.SJ.out.tab | \
  awk '($7 >= 3 && $5 > 0)'| cut -f1-6| sort| uniq > SJ/merged_UMR_3.SJ.tab

# Step 4: Genome re-indexing with the filtered set of SJ
# Note: Run for all samples across group comparisons
mkdir reindexed_genome
STAR --runThreadN 16 \
            --runMode genomeGenerate \
            --genomeDir reindexed_genome \
            --genomeFastaFiles hg38.genome.fasta \
            --sjdbOverhang $sjdbOverhang \
            --sjdbFileChrStartEnd SJ/merged_UMR_3.SJ.tab

# Step 5: per-sample 2nd-mapping pass
STAR --genomeDir reindexed_genome \
                    --readFilesIn trimmed/sample1_G1_R1.fastq.gz trimmed/sample1_G1_R2.fastq.gz \
                    --readFilesCommand zcat \
                    --outFileNamePrefix aligned_pass2/sample1_G1. \
                    --outSAMtype BAM SortedByCoordinate \
                    --outFilterType BySJout \
                    --outFilterIntronMotifs RemoveNoncanonical \
                    --runThreadN 16

# At the end, we have one BAM file for each sample after 2nd-mapping pass
```

#### 2. Exon-junction read quantification (Rsubread featureCounts)

**Input:**\
`hg38.flat_exon.saf` (Flattened and merged exon annotation), `hg38.SJdatabase.tsv` (Junction database),\
`sample1_G1.Aligned.sortedByCoord.out.bam`, `sample2_G1.Aligned.sortedByCoord.out.bam`, `sample1_G2.Aligned.sortedByCoord.out.bam`, `sample2_G2.Aligned.sortedByCoord.out.bam` in `aligned_pass2` folder

**Output:** Exon-junction count table and annotation is stored in `IE_J_count` and `IE_J_annot` R objects.

```r
# Install Rsubread if not yet installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Rsubread")

library(Rsubread)
BAM_files <- list("aligned_pass2/sample1_G1.Aligned.sortedByCoord.out.bam",
                  "aligned_pass2/sample2_G1.Aligned.sortedByCoord.out.bam",
                  "aligned_pass2/sample1_G2.Aligned.sortedByCoord.out.bam",
                  "aligned_pass2/sample2_G2.Aligned.sortedByCoord.out.bam")

message('Quantify internal exon + exon-exon junction reads ...')
count <- Rsubread::featureCounts(BAM_files, # a list of input BAM files
                                  annot.ext="hg38.flat_exon.saf", # specify relative path to merged and flattened exon annotation file previously generated
                                  useMetaFeatures=FALSE, # summarize exon-level reads
                                  nonSplitOnly=TRUE, splitOnly=FALSE, # quantify internal exon reads
                                  juncCounts=TRUE # quantify exon-exon junction reads
                                  isPairedEnd=FALSE # is paired-end?)

message('Internal exon count and annotation ...')
IE_count <- count$counts
IE_annot <- cbind(count$annotation, Region="Exon", annotated=1)

message('Exon-exon junction count and annotation ...')
J_count <- count$counts_junction[!is.na(count$counts_junction$PrimaryGene),]
J_count$juncID <- paste(J_count$Site1_chr, 
                        J_count$Site1_location,
                        J_count$Site2_location,
                        sep="_")

message('Junction reannotation using junction database ...')
SJ_database <- read.table("hg38.SJdatabase.tsv", header=T) # Specify relative path to SJ database file previously created
uniq_SJ <- SJ_database[SJ_database$freq==1,]
m1 <- match(J_count$juncID, uniq_SJ$juncID)
J_count$PrimaryGene <- ifelse(!is.na(m1), uniq_SJ$geneID[m1], J_count$PrimaryGene)

message("Whether junctions are annotated or not ...")
m2 <- match(J_count$juncID, SJ_database$juncID)
J_count$annotated <- ifelse(!is.na(m2), 1, 0)

message("Processing final junction count table ...")
J_annot <- data.frame(
  GeneID=J_count$PrimaryGene,
  Chr=J_count$Site1_chr,
  Start=J_count$Site1_location,
  End=J_count$Site2_location
)
m <- match(J_annot$GeneID, IE_annot$GeneID)
Strand <- IE_annot$Strand[m]
J_annot <- cbind(J_annot, Strand=Strand, Length=1, Region="Junction", annotated=J_count$annotated)
J_count <- data.frame(J_count[, 9:(ncol(J_count)-2)], row.names = NULL, check.names = FALSE)

message('Combine internal exon counts and junction counts into an exon-junction count table ...')
IE_J_count <- rbind(IE_count, J_count)
IE_J_annot <- rbind(IE_annot, J_annot)
```

An example of R object outputs:

`IE_J_count` contains counts of each feature (exon/junction) in the `IE_J_annot` object with each column representing the corresponding sample (e.g, sample1_G1, sample2_G1, sample1_G1, sample2_G2).\
**Important note:**\
SampleID will do not appear in this object. But, we will know which column belongs to which sample based on the order of BAM files loaded in previous step (The order of samples will be the same as the order of BAM files).\
The row number in `IE_J_count` corresponds to the same row number in `IE_J_annot` object. So, the number of rows in 2 objects are always the same (tells us how many features (exons and junctions) we have.

```tsv
7587  5384  6408  6617
732  567  731  709
329  268  282  294
1  0  0  0
595  479  491  610
5  5  0  1
```

`IE_J_annot` contains 8 columns: GeneID, Chr, Start, End, Strand, Length, Region, Annotated, which provides genomic information for each feature above.
```tsv
GeneID  Chr  Start  End  Strand  Length  Region  Annotated
ENSMUSG00000000001.5  chr3  108014596  108016632  -  2037  Exon  1
ENSMUSG00000000001.5  chr3  108016719  108016928  -  210  Exon  1
ENSMUSG00000000001.5  chr3  108019251  108019404  -  154  Exon  1
ENSMUSG00000000001.5  chr3  108016623  108016719  -  1  Junction  0
ENSMUSG00000000001.5  chr3  108016632  108016719  -  1  Junction  1
ENSMUSG00000000001.5  chr3  108016632  108019251  -  1  Junction  0
```

**Field descriptions:**

`GeneID`, `Chr`, `Strand`: Gene details\
`Region`: if the feature is exon or junction\
`Start`, `End`, `Length`: Start, end position and length of the feature (if feature is junction, length is set to 1)\
`Annotated`: if the feature is annotated or novel (0 means 'novel')

#### 3. Differential exon-junction usage (edgeR diffSpliceDGE)

**Input:** `IE_J_count` (exon-junction counts), `IE_J_annot` (exon-junction annotation), `group` (group names), `contr` (contrast matrix), `design` (design matrix) R objects

**Output:** A list of differentially spliced genes by Simes/F-test and a list of differential exon and junction regions by exon test

First, manually create a tab-separated TSV file `target.tsv` that contains sampleID and groups of samples information like below:

|sampleID|group|
|----|----|
|sample1_G1|G1|
|sample2_G1|G1|
|sample1_G2|G2|
|sample2_G2|G2|

```r
# Install edgeR if not yet installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("edgeR")

library(edgeR)

message("Specify groups of samples and build model/contrast matrix for differential analysis ...")
target <- read.table("target.tsv", sep="\t", header=TRUE, stringsAsFactors=FALSE)
group <- factor(target$group)
samples <- target$sampleID

design <- model.matrix(~ 0 + group)
colnames(design) <- gsub("group", "", colnames(design))

p <- "G2-G1" # specify groups to be compared (in this case, G1 is the reference group)
contr <- makeContrasts(p, levels=design)

message("Constructing DGElist object ...")  
y <- DGEList(counts=IE_J_count, genes=IE_J_annot, group=group)
colnames(y) <- gsub("[.].*$", "", colnames(y))

message("Filtering exons with low mapping reads ...")
keep <- filterByExpr(y, group=group)
y <- y[keep, , keep.lib.sizes=FALSE]

message("Normalizing lib sizes ...")
y <- normLibSizes(y)

message("Estimating dispersion ...")
y <- estimateDisp(y, design, robust=TRUE)

message("Fitting GLM-QL model for design matrix ...")
fit <- glmQLFit(y, design, robust=TRUE)

message("Running diffSpliceDGE ..")
sp <- diffSpliceDGE(fit, contrast=contr, geneid="GeneID", exonid="Start")
```

**Output:**
Results provides ranked genes or exons by evidence for differential splicing, based on sorted adjusted p-values. 
The exon-level tests test for differences between each exon and all the exons for the same gene. 
The gene-level tests test for any differences in exon usage between experimental conditions. 
The Simes method processes the exon-level p-values to give an overall call of differential splicing for each gene. 
It returns the minimum Simes-adjusted p-values for each gene. 
The gene-level tests are likely to be powerful for genes in which several exons are differentially splices. 
The Simes p-values is likely to be more powerful when only a minority of the exons for a gene are differentially spliced. 

- **DEU genes from gene-level Simes test**
```r
topSpliceDGE(sp, test="Simes")
```
|GeneID|Chr|Strand|Symbol|NExons|P.Value|FDR|
|----|----|----|----|----|----|----|
|ENSMUSG00000028337.15|chr4|-|Coro2a|27|5.63409202936316e-18|7.44150875238286e-14|
|ENSMUSG00000052033.14|chr2|+|Pfdn4|11|7.025331391432e-17|4.63952885090169e-13|
|ENSMUSG00000031075.20|chr7|-|Ano1|62|1.87426860695308e-12|8.25177992021212e-09|
|ENSMUSG00000022106.15|chr14|+|Rcbtb2|34|3.01085278939628e-12|9.9418359105865e-09|

**Field descriptions:**

`GeneID`, `Chr`, `Strand`, `Symbol`: Gene details\
`NExons`: The total number of exons and junctions of the gene\
`P.value`: p-value of Simes test\
`FDR`: False discovery rate

- **DEU genes from gene-level F-test**
```r
topSpliceDGE(sp, test="gene")
```
|GeneID|Chr|Strand|Symbol|NExons|gene.F|P.Value|FDR|
|----|----|----|----|----|----|----|----|
|ENSMUSG00000035202.9|chr9|+|Lars2|35|16.1933702582295|1.11593822152576e-34|1.47393120299123e-30|
|ENSMUSG00000028337.15|chr4|-|Coro2a|27|12.5816242961678|7.64969721594055e-23|5.05186004140714e-19|
|ENSMUSG00000022091.7|chr14|-|Sorbs3|38|7.89331056107204|4.92908588386424e-21|2.17011221180263e-17|
|ENSMUSG00000021224.16|chr12|-|Numb|28|10.4648157641059|1.90223935574928e-20|6.28119435268413e-17|

**Field descriptions:**

`GeneID`, `Chr`, `Strand`, `Symbol`: Gene details\
`NExons`: The total number of exons and junctions of the gene\
`gene.F`: F-statistics for gene\
`P.value`: p-value of F-test\
`FDR`: False discovery rate

- **Differentially used exons and splice junctions from exon-level test**
```r
topSpliceDGE(sp, test="exon")
```
|GeneID|Chr|Start|End|Strand|Length|Region|annotated|Symbol|logFC|exon.F|P.Value|FDR|
|----|----|----|----|----|----|----|----|----|----|----|----|----|
|ENSMUSG00000028337.15|chr4|46576626|46581702|-|1|Junction|0|Coro2a|2.23199936023307|118.83727786305|2.08670075161598e-19|4.90410150542534e-14|
|ENSMUSG00000028337.15|chr4|46576626|46583643|-|1|Junction|0|Coro2a|2.19408998911786|114.148756203849|6.73068952852665e-19|7.90913230462874e-14|
|ENSMUSG00000052033.14|chr2|170338348|170338519|+|172|Exon|1|Pfdn4|2.72395320118021|173.077184263658|6.38666490130182e-18|5.00324941703083e-13|
|ENSMUSG00000031075.20|chr7|144292097|144292329|-|233|Exon|1|Ano1|-2.4622662045429|64.9322523521431|3.0230138821824e-14|1.77614913387215e-09|

**Field descriptions:**

`GeneID`, `Chr`, `Strand`, `Symbol`: Gene details\
`Region`: Whether the region is junction or exon\
`Start`, `End`: If region is exon, start/end coordinator of that exonic region; If region is junction, position of two splice sites of that junction\
`annotated`: whether exon/junction is annotated or novel (0 means 'novel')\
`logFC`: log2 fold-change of one exon vs all the exons for the same gene\
`exon.F`: F-statistics for exon/junction\
`P.value`: p-value of exon-level test\
`FDR`: False discovery rate

#### 4. Visualisation of significant gene with differential splicing (DS gene)

Here we provide an example of how DS gene can be visualised.

<p align="center">
  <img src="figures/FGFR1_DS_gene.png" alt="FGFR1 DS gene" width="400"/>
</p>

**To draw the upper panel of DEJU-edgeR**

```bash
# Please download the R script below:

# If run DEJU-edgeR
wget https://raw.githubusercontent.com/TamPham271299/DEJU/refs/heads/main/code/analysis/plotJunc3_diffSpliceDGE.R

source("plotJunc3_diffSpliceDGE.R")
plotJunc(sp, geneid=g, genecol="Symbol", annotation=IE_J$IE_annot)

## If run DEJU-limma
# wget https://raw.githubusercontent.com/TamPham271299/DEJU/refs/heads/main/code/analysis/plotJunc3.R

# source("plotJunc3.R")
# plotJunc(sp, geneid=g, genecol="Symbol", annotation=IE_J$IE_annot)
```

**To draw Sashimi plot (bottom panel)**

- **Step 1: generate annotation file of gene regions and indexed BAM files if not yet**

**Input:** `hg38.genome.gtf`, `sample1_G1.Aligned.sortedByCoord.out.bam`, `sample2_G1.Aligned.sortedByCoord.out.bam`, `sample1_G2.Aligned.sortedByCoord.out.bam`, `sample2_G2.Aligned.sortedByCoord.out.bam` in `aligned_pass2` folder

**Output:** `hg38.genome.genes.bed`, `sample1_G1.Aligned.sortedByCoord.out.bam.bai`, `sample2_G1.Aligned.sortedByCoord.out.bam.bai`, `sample1_G2.Aligned.sortedByCoord.out.bam.bai`, `sample2_G2.Aligned.sortedByCoord.out.bam.bai` in `aligned_pass2` folder

```bash
# Prepare gene coordinates from GTF file (have to be downloaded from Gencode to run this code successfully)
# Otherwise, need to adjust the code accordingly 
(echo -e "chr\tstart\tend\tgeneID\tscore\tstrand\tgeneSymbol"; \
  grep -v "^#" hg38.genome.gtf| \
  awk 'BEGIN {FS=OFS="\t"} {split($9, a, "; ")} {if ($3=="gene") print $1, $4-1, $5, a[1], "0", $7, a[3]}'| \
  sed -e 's/gene_id "//g' -e 's/gene_name "//g' -e 's/"//g') > hg38.genome.genes.bed

# Index bam files (if not yet)
# e.g
# Index BAM files for all 4 samples stored in aligned_pass1 folder
for SAMPLE in sample1_G1 sample2_G1 sample1_G2 sample2_G2; do
  echo $SAMPLE
  BAM="aligned_pass2/${SAMPLE}.Aligned.sortedByCoord.out.bam"
  samtools index $BAM
done
```

- **Step 2: Generate BAM files that contain gene regions of interest**

For more examples, please visit [code/analysis/make_bam_for_visualization.sh](code/analysis/make_bam_for_visualization.sh)

**Input:** `sample1_G1.Aligned.sortedByCoord.out.bam`, `sample2_G1.Aligned.sortedByCoord.out.bam`, `sample1_G2.Aligned.sortedByCoord.out.bam`, `sample2_G2.Aligned.sortedByCoord.out.bam` and their index `.bam.bai` files in `aligned_pass2` folder

**Output:** `sample1_G1.Fgfr1.bam`, `sample2_G1.Fgfr1.bam`, `sample1_G2.Fgfr1.bam`, `sample2_G2.Fgfr1.bam` and their index `.bam.bai` files in `${gene}_${geneID}_${d}` folder (e.g. Fgfr1_ENSMUSG00000031565.19_1000 as an example)

```bash
# Specify DS gene to visualise and set the upstream/downstream distance (in base pairs) to include before the gene’s start coordinate and after its end coordinate.
gene="Fgfr1"
d="1000"

gene_region=$(awk -v g=$gene -v d=$d 'BEGIN {FS=OFS="\t"} {if($7==g) print $1, $2-d, $3+d, $4}' hg38.genome.genes.bed)
chr=$(echo "$gene_region"| cut -f1)
start=$(echo "$gene_region"| awk -v d=$d '{print $2-d}')
end=$(echo "$gene_region"| awk -v d=$d '{print $3+d}')
position=${chr}:${start}-${end}
geneID=$(echo "$gene_region"| cut -f4)

mkdir -p ${gene}_${geneID}_${d}

for SAMPLE in sample1_G1 sample2_G1 sample1_G2 sample2_G2; do
    echo $SAMPLE
    BAM="aligned_pass2/${SAMPLE}.Aligned.sortedByCoord.out.bam"
    samtools view -b -h $BAM $position > ${gene}_${geneID}_${d}/${SAMPLE}.${gene}.bam
    samtools index ${gene}_${geneID}_${d}/${SAMPLE}.${gene}.bam
done
```

- **Step 3: Generate Sashimi plots using Gviz for a neat and nice plot**

Please refer to [Gviz documentation](https://bioconductor.org/packages/devel/bioc/vignettes/Gviz/inst/doc/Gviz.html) for more details.)
Or, we can also visualise Sashimi plot using IGV.

We can also visit [code/analysis/main_figs_codes.R - section Figure 4B](https://github.com/TamPham271299/DEJU/blob/main/code/analysis/main_figs_codes.R#L413-L506) or [code/analysis/supp_figs_codes.R - section Figure S9]([code/analysis/supp_figs_codes.R](https://github.com/TamPham271299/DEJU/blob/main/code/analysis/supp_figs_codes.R#L1155-L1518)) for more Sashimi junction plots (shown in the paper as Figure 4B and the supplementary document [12859_2025_6210_MOESM1_ESM.pdf](https://static-content.springer.com/esm/art%3A10.1186%2Fs12859-025-06210-4/MediaObjects/12859_2025_6210_MOESM1_ESM.pdf) as Figure S9).

```R
# Install Gviz if not yet installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Gviz")

library(Gviz)
options(ucscChromosomeNames=FALSE)

message("Specify gene position and DS region to be zoomed in ...")
# Below is an example to visualise Fgfr1 gene as a DS gene of 2 groups (mouse LP and ML cell types)
geneSymbol <- "Fgfr1"
geneID <- "ENSMUSG00000031565.19"
chr <- "chr8"
g_from <- 26002669 # first coordinate of gene region
g_to <- 26066734 # last coordinate of gene region
z_from <- 26052007 # first coordinate of region to zoom in
z_to <- 26059369 # last coordinate of region to zoom in
d <- 1000 # Upstream/Downstream distance to the first/last coordinate of the gene

message("Load alignmnent track from bam files ...")
samples <- c("sample1_G1", "sample2_G1", "sample1_G2", "sample2_G2")
OUTPUT <- paste0(geneSymbol, "_", geneID, "_", d, "/")
PE_bam_files <- paste0(OUTPUT, samples, ".", geneSymbol, ".bam")

message("Start an alTrack list to store alignment tracks of 4 samples ...")
alTrack <- list()

for (idx in 1:2) {
  s <- samples[idx]
  # set alignment track for sample1_G1 and sample2_G1 (with BLUE color)
  alTrack[[s]] <- AlignmentsTrack(PE_bam_files[idx], 
                                  isPaired = TRUE, 
                                  fill.coverage="blue", 
                                  col.sashimi="blue", 
                                  fill="white", 
                                  fontsize=8, cex.axis=1, col="blue")
}

for (idx in 3:4) {
  s <- samples[idx]
  # set alignment track for sample1_G2 and sample2_G2 (with GREEN color)
  alTrack[[s]] <- AlignmentsTrack(PE_bam_files[idx], 
                                  isPaired = TRUE, 
                                  fill.coverage="green", 
                                  col.sashimi="green", 
                                  fill="white", 
                                  fontsize=8, cex.axis=1, col="green")
}

message("Load transcript annotation track ...")
# Visit https://genome.ucsc.edu/cgi-bin/hgTables for more information of genome of interest
# Here we use the latest version of human genome to visualise
knownGenes <- UcscTrack(genome = "hg38", chromosome = chr, 
                        track = "All GENCODE V48", table="wgEncodeGencodeCompV48", # use UCSC comprehensive genome annotation of Gencode hg38 release 48
                        from = g_from, to = g_to, # Specify start/end position of gene region of interest
                        trackType = "GeneRegionTrack", 
                        rstarts = "exonStarts", rends = "exonEnds", 
                        gene = "name", symbol = "name", 
                        transcript = "name", strand = "strand", 
                        fill = "black", name = "ALL GENCODE V48", col="black")

message("Load genome axis track ...")
idxTrack <- IdeogramTrack(genome="hg38", chromosome=chr)
axTrack <- GenomeAxisTrack()

message("Combine all tracks ...")
plotTracks(c(idxTrack, axTrack, knownGenes, alTrack),
            from = z_from, to = z_to, # Specify start and end region to zoom in
            chromosome = chr, type = c("coverage", "sashimi"),
            sashimiNumbers=TRUE,
            showTitle = FALSE,
            # sashimiScore=9, # optionally, set sashimiScore to reduce some unncessary junctions to be visualised
            lwd.sashimiMax=2,
            sizes = c(0.1,0.2,0.3,rep(0.3, 4)), # To adjust sizes of each panel
            lwd.title=2,
            lwd.border.title=4,
            background.title="white",
            col.axis="black",
            fontcolor="black",
            col.line="black")
```
