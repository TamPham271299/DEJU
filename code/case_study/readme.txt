This directory contains scripts used to analyze RNA-seq datasets of case study to detect DEU genes, 
including:
- RNA-seq pre-processing and analysis; and 
- DEU/DEJU analyses.

## Directory structure

case_study/
├── run_00_1.sh                  # Download FASTQ files
├── run_01_1.sh                # FASTQC
├── run_01_2.sh                # Read trimming + FASTQC
├── run_02_0.sh                # Index the genome (run only once) - STAR
├── run_02_1.sh                # 1st-pass mapping - STAR
├── run_02_2_3.sh                # SJ filtering + re-index genome - STAR
├── run_02_4.sh                    # 2nd-pass mapping - STAR
├── run_03_1.sh                    # run edgeR::diffSpliceDGE
├── run_03_2.sh                    # run limma::diffSplice
├── run_03_3.sh                    # run DEXSeq
├── run_03_4.sh                   # run JunctionSeq
├── run_GSE227748.sh               # Main script to run all steps (run_00_1.sh -> run_03_4.sh) for null simulations
└── readme.txt                  # Main documentation for the repository

0 - Please execute bash scripts by running chmod +x *.sh
1 - run the SLURM bash script run_GSE227748.sh with, for example, sbatch run_GSE227748.sh 
    on your HPC (*)

(*) note that, instead of running on-the-go script, we can run each numbered step one-by-one.
