#!/bin/bash 
#SBATCH --job-name=00_customized_transcriptome
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --time=00:5:00
#SBATCH --output=slurm-%A.out

module load R/4.3.3

LOG_00="../../data/simulation/customized_transcriptome/00_customized_transcriptome.log"

echo "====== `date`: Customizing transcriptome with designed splicing pattern ======" > $LOG_00

echo "====== `date`: Designing splicing pattern ======" >> $LOG_00
Rscript --no-save --no-restore --verbose customized_transcriptome.R >> $LOG_00 2>&1

echo "====== `date`: From genomic intervals to transcriptome" >> $LOG_00
bash customized_transcriptome.sh >> $LOG_00 2>&1

echo "====== `date`: Customizing transcriptome is done! ======" >> $LOG_00
