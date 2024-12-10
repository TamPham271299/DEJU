# [Differential Exon-Junction Usage (DEJU)]

This is the repository for the code used to perform simulation and case-study DEU analysis and generate figures for the following paper titled "Incorporating exon-exon junction reads enhances differential splicing detection". 

## Introduction

The paper introduce a DEJU analysis workflow implementing a STAR-Rsubread-edgeR-limma framework to analyze DEJU and detect differential splicing between two groups of conditions. Here is a schematic presentation of our proposed DEJU workflow.

[Figure of DEJU workflow]

The paper also benchmarked the DEJU analysis workflow (DEJU-edgeR, DEJU-limma) against the existing DEU analysis workflow (DEU-edgeR, DEU-limma) and other popular methods (DEXSeq, JunctionSeq) based on the simulation RNA-seq datasets.

## Citation

If you are using code or pipelines from this repository, please consider citing our associated article:

Pham, M. T., Milevskiy, M. J. G., Visvader, J. E., Chen, Y. Incorporating exon-exon junction reads enhances differential splicing detection. ...

## Repository Structure

```plaintext
DEJU/
├── annotation/
│   ├── duplicated_sequences.tsv/    # duplicate sequences discarded for the simulation process
│   └── processed/        # Data after cleaning and preprocessing
├── code/
│   ├── annotation_dl    # Script to prepare reference genome datasets of mm39 (GTF, FASTA, SAF...) 
│   ├── alignment      # Script to map reads to reference genome mm39
│   ├── DEU            # Script to quantify exon-junction counts, and detect DEU genes by 6 benchmarked methods
│   ├── simulation    # Main script to generate customized transcriptome and analyze simulation data
│   ├── fastq_dl      # Script to download FASTQ files
│   ├── case_study    # Main script to analyze case-study data
│   └── analysis      # Main script to perform overall analysis, create figures for the paper
├── tools/
│   └── Subread_to_DEXSeq/    # External script written by ... to generate exon counts for DEXSeq
└── README.md                  # Main documentation for the repository

