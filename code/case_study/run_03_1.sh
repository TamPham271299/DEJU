#!/bin/bash 
#SBATCH --job-name=03_1_edgeR_diffSpliceDGE
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --time=00:45:00
#SBATCH --output=slurm-%A.out

MODE="case_study"
DIR="$DIR"
REF="../../annotation/"
# target="$DIR/target/target.tsv"
target="$DIR/target/target.${pair}.tsv"
pair="$pair"
isPairedEnd="$isPairedEnd"
fdr_cutoff=0.05

### Running DS analysis using diffSpliceDGE
LOG_03_1=$DIR/log/03_1_edgeR_diffSpliceDGE.log

export MODE="$MODE"
export DIR="$DIR" 
export REF="$REF" 
export target="$target" 
export pair="$pair"
export isPairedEnd="$isPairedEnd"
export fdr_cutoff="$fdr_cutoff"
export LOG_03_1="$LOG_03_1"

chmod +x "../DEU/03_1_edgeR_diffSpliceDGE.sh"

echo "====== `date`: edgeR:diffSpliceDGE ======" > $LOG_03_1
bash "../DEU/03_1_edgeR_diffSpliceDGE.sh"
echo "====== `date`: edgeR:diffSplice is done! ======" >> $LOG_03_1
