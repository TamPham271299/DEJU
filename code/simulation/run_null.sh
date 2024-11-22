#!/bin/bash 

export noOfSim=1
export sim_mode="null"

### Simulations
job_01_1_id=$(sbatch --parsable --export=noOfSim="$noOfSim" run_01_1.sh)

### STAR 1-pass mapping
job_02_1_id=$(sbatch --parsable --dependency=afterok:$job_01_1_id run_02_1.sh)

### Genome reindexing with filtered SJ
job_02_2_3_id=$(sbatch --parsable --dependency=afterok:$job_02_1_id run_02_2_3.sh)

### STAR 2-pass mapping
job_02_4_id=$(sbatch --parsable --dependency=afterok:$job_02_2_3_id run_02_4.sh)

### DS using edgeR diffSpliceDGE
job_03_1_id=$(sbatch --parsable --export=noOfSim="$noOfSim" --dependency=afterok:$job_02_2_4_id run_03_1.sh)

### DS using limma diffSplice
job_03_2_id=$(sbatch --parsable --export=noOfSim="$noOfSim" --dependency=afterok:$job_02_2_4_id:$job_03_1_id run_03_2.sh)

### DS using DEXSeq
job_03_3_id=$(sbatch --parsable --export=noOfSim="$noOfSim" --dependency=afterok:$job_02_2_4_id run_03_3.sh)

### DS using JunctionSeq
job_03_4_id=$(sbatch --parsable --export=noOfSim="$noOfSim" --dependency=afterok:$job_02_2_4_id run_03_4.sh)
