# Differential Exon-Junction Usage (DEJU)

This is the repository for the code used to perform analysis and generate figures for the following paper titled "Incorporating exon-exon junction reads enhances differential splicing detection". 

### Introduction

The paper introduced a DEJU analysis workflow implementing a STAR-Rsubread-edgeR-limma framework to analyze DEJU and detect differential splicing between two groups of conditions. Here is a schematic presentation of our proposed DEJU workflow.

![DEJU_workflow](figures/fig_1.png)

The paper also benchmarked the DEJU analysis workflow (DEJU-edgeR, DEJU-limma) against the existing DEU analysis workflow (DEU-edgeR, DEU-limma) and other popular tools (DEXSeq, JunctionSeq) based on simulated RNA-seq datasets. We also performed DEU analysis on RNA-sequencing experiments of [NCBI GEO database (GSE227748)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE227748) using our proposed method and other DEU analysis pipelines benchmarked in the study.

### Citation

If you are using code or pipelines from this repository, please consider citing our associated article:

Pham, M. T., Milevskiy, M. J. G., Visvader, J. E., Chen, Y. Incorporating exon-exon junction reads enhances differential splicing detection. ...

### Repository Structure

```plaintext
DEJU/
├── annotation/
│   ├── gene_info_GRCm39_M32_ensembl109.tsv    # Gene information
│   └── duplicated_sequences.tsv/    # duplicate sequences discarded for the simulation process
├── code/
│   ├── annotation_dl    # Script to prepare reference genome datasets of mm39 (GTF, FASTA, exon-junction SAF...) 
│   ├── alignment      # Script to map reads to reference genome mm39
│   ├── DEU            # Script to quantify exon-junction counts, and detect DEU genes by 6 benchmarked methods
│   ├── simulation    # Main script to generate customized transcriptome, generate and analyze simulation data
│   ├── fastq_dl      # Script to download FASTQ files
│   ├── case_study    # Main script to analyze case-study data
│   └── analysis      # Main script to perform overall analysis, create figures for the paper
├── tools/
│   └── Subread_to_DEXSeq/    # External scripts to generate exon counts for DEXSeq
├── figures/
│   └── fig_S1.png    # DEJU analysis worflow
└── README.md                  # Main documentation for the repository

