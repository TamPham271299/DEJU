This directory contains scripts used to run the complete simulation pipeline, 
including:
- customizing transcriptome with designed splicing patterns;
- simulation of ground truth data with and without DEU genes (null simulation); 
- RNA-seq analysis; and 
- differential exon-junction usage.

0 - Please execute bash scripts by running chmod +x *.sh
1 - run the SLURM bash script run_00.sh with, for example, sbatch run_00.sh on 
    your HPC running SLURM scheduler to customize transcriptome.
2 - In RNA-seq analysis, we use STAR aligner for read mapping. Please index 
    reference genome using STAR aligner, for example, run sbatch run_02_0.sh
3 - run the file parameters_DEU.R and parameters_null.R with, for example, 
    the command Rscript parameters_DEU.R and Rscript parameters_null.R
    to create parameters_DEU.txt and parameters_null.txt, respectively.
4 - run the SLURM bash script run_DEU.sh with, for example, sbatch run_DEU.sh 
    on your HPC to simulate RNA-seq reads with DEU genes followed by read 
    mapping and DEU analysis
5 - run the SLURM bash script run_null.sh with, for example, sbatch run_null.sh
    on your HPC to simulate RNA-seq reads without DEU genes followed by 
    read mapping and DEU analysis
6 - run the SLURM bash script run_test.sh with, for example, sbatch run_test.sh
    on your HPC to test for runtime and computing resource consumption by all 
    compared DEU analysis tools.

Step 1 above will create a customized transcriptome with splicing patterns 
studied in the article. This transcriptome will be used for all simulations.

Step 3 above will create a parameter.txt file with all combinations of
scenarios for the DEU and null simulation. 
- The file parameter_DEU.txt generated in Step 2 will be read when you run 
  run_DEU.sh. It provides the necessary parameters for each DEU simulation. 
- Besides, the file parameter_null.txt generated in Step 2 will be read when 
  you run run_null.sh. It provides the necessary parameters for each null 
  simulation.

Steps 4-5 are run on a per-experiment basis, and each experiment is run
in a SLURM batch job triggered with array batch jobs.
To run the simulation pipeline, do:

Please run step 1-3 prior steps 4-5 which can be run independently. 
Please run step 6 after step 4 is done.

Computing time per step **
1 - approximately 30 seconds
2 - approximately 1 hour
3 - approximately 30 seconds
4 - approximately six days*
5 - approximately 1 day*
6 - approximately 6 hours*

*: note that these steps use parallel computing in several parts of the
    scripts, please adjust accordingly.
**: note that the maximum computing time for 1 simulation.

For 1 simulation (balanced_10_75_3)
run_00 (30s)
run_01_1 (1 hrs)
run_02_0 (...)
run_02_1 (2 hrs)
run_02_2_3 (...)
run_02_4 (...)
run_03_1 (...)
run_03_2 (...)
run_03_3 (...)
run_03_4 (...)

