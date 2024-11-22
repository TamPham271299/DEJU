#!/bin/bash 

Rscript --no-save --no-restore --verbose ../DEU/04_1_performance_analysis.R \
          --DIR=$DIR \
          --DEU_map=$DEU_map \
          --fdr_cutoff=$fdr_cutoff \
          --seed=$seed \
          --workers=$workers \
          --noOfSim=$noOfSim