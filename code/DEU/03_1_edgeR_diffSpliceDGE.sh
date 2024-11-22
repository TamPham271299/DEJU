#!/bin/bash 

module load R/4.3.3

### Running DS using diffSpliceDGE
if [ "$MODE" == "simulation" ]; then

  Rscript --no-save --no-restore --verbose ../DEU/03_1_edgeR_diffSpliceDGE.R \
    --mode=$MODE \
    --DIR=$DIR \
    --REF=$REF \
    --target=$target \
    --pair=$pair \
    --isPairedEnd=$isPairedEnd \
    --fdr_cutoff=$fdr_cutoff \
    --seed=$seed \
    --workers=$workers \
    --noOfSim=$noOfSim >> $LOG_03_1 2>&1

elif [ "$MODE" == "case_study" ]; then

  Rscript --no-save --no-restore --verbose ../DEU/03_1_edgeR_diffSpliceDGE.R \
    --mode=$MODE \
    --DIR=$DIR \
    --REF=$REF \
    --target=$target \
    --pair=$pair \
    --isPairedEnd=$isPairedEnd \
    --fdr_cutoff=$fdr_cutoff >> $LOG_03_1 2>&1
  
else

  echo "Invalid mode specified. Choose 'simulation' or 'case_study'." >> $LOG_03_1

fi

################################################################################
# module load R/4.3.3
# 
# IFS=, read -r -a pair <<< "$pair"
# 
# ### Running DS using diffSpliceDGE
# if [ "$MODE" == "simulation" ]; then
# 
#   echo "====== `date`: Start running diffSpliceDGE ======" >> $LOG_03_1
#   for p in "${pair[@]}"; do
#     echo "====== `date`: $p ======" >> $LOG_03_1
#     Rscript --no-save --no-restore --verbose ../DEU/03_1_edgeR_diffSpliceDGE.R \
#       --mode=$MODE \
#       --DIR=$DIR \
#       --REF=$REF \
#       --target=$target \
#       --pair=$p \
#       --isPairedEnd=$isPairedEnd \
#       --fdr_cutoff=$fdr_cutoff \
#       --seed=$seed \
#       --workers=$workers \
#       --noOfSim=$noOfSim >> $LOG_03_1 2>&1
#   done
#   echo "====== `date`: Running diffSpliceDGE is done! ======" >> $LOG_03_1
# 
# elif [ "$MODE" == "case_study" ]; then
# 
#   echo "====== `date`: Start running diffSpliceDGE ======" >> $LOG_03_1
#   for p in "${pair[@]}"; do
#     echo "====== `date`: $p ======" >> $LOG_03_1
#     Rscript --no-save --no-restore --verbose ../DEU/03_1_edgeR_diffSpliceDGE.R \
#       --mode=$MODE \
#       --DIR=$DIR \
#       --REF=$REF \
#       --target=$target \
#       --pair=$p \
#       --isPairedEnd=$isPairedEnd \
#       --fdr_cutoff=$fdr_cutoff >> $LOG_03_1 2>&1
#   done
#   echo "====== `date`: diffSpliceDGE is done! ======" >> $LOG_03_1
# 
# else
# 
#   echo "Invalid mode specified. Choose 'simulation' or 'case_study'." >> $LOG_03_1
# 
# fi
