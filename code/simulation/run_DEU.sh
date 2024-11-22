#!/bin/bash 

export noOfSim=20
export sim_mode="DEU"

### Part 1: Simulating RNA-seq reads
job_01_1_id=$(sbatch --parsable run_01_1.sh)

### Part 2: Read alignment using STAR
# 1-pass mapping
job_02_1_id=$(sbatch --parsable --dependency=afterok:$job_01_1_id run_02_1.sh)

### Genome reindexing with filtered SJ
job_02_2_3_id=$(sbatch --parsable --dependency=afterok:$job_02_1_id run_02_2_3.sh)

### 2-pass mapping
job_02_4_id=$(sbatch --parsable --dependency=afterok:$job_02_2_3_id run_02_4.sh)

### Part 3: DEU analysis
# edgeR::diffSpliceDGE
job_03_1_id=$(sbatch --parsable --dependency=afterok:$job_02_2_4_id run_03_1.sh)

# limma::diffSplice
job_03_2_id=$(sbatch --parsable --dependency=afterok:$job_02_2_4_id:$job_03_1_id run_03_2.sh)

# DEXSeq
job_03_3_id=$(sbatch --parsable --dependency=afterok:$job_02_2_4_id run_03_3.sh)

# JunctionSeq
job_03_4_id=$(sbatch --parsable --dependency=afterok:$job_02_2_4_id run_03_4.sh)

### Part 4: Benchmarking analysis
job_04_1_id=$(sbatch --parsable --dependency=afterok:$job_03_1_id:$job_03_2_id:$job_03_3_id:$job_03_4_id run_04_1.sh)
