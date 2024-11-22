#!/bin/bash 

module load R/4.3.3

### Running DS using diffSplice
if [ "$MODE" == "simulation" ]; then
  
  Rscript --no-save --no-restore --verbose ../DEU/03_2_limma_diffSplice.R \
    --mode=$MODE \
    --DIR=$DIR \
    --REF=$REF \
    --target=$target \
    --pair=$pair \
    --isPairedEnd=$isPairedEnd \
    --fdr_cutoff=$fdr_cutoff \
    --seed=$seed \
    --workers=$workers \
    --noOfSim=$noOfSim >> $LOG_03_2 2>&1

elif [ "$MODE" == "case_study" ]; then
  
  Rscript --no-save --no-restore --verbose ../DEU/03_2_limma_diffSplice.R \
    --mode=$MODE \
    --DIR=$DIR \
    --REF=$REF \
    --target=$target \
    --pair=$pair \
    --isPairedEnd=$isPairedEnd \
    --fdr_cutoff=$fdr_cutoff >> $LOG_03_2 2>&1
  
else

  echo "Invalid mode specified. Choose 'simulation' or 'case_study'." >> $LOG_03_2
  
fi
