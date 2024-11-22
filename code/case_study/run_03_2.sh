#!/bin/bash 
#SBATCH --job-name=03_2_limma_diffSplice
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --time=00:20:00
#SBATCH --output=slurm-%A.out

MODE="case_study"
# DIR="../../data/case_study/$projName/"
DIR="$DIR"
REF="../../annotation/"
target="$DIR/target/target.tsv"
pair="$pair"
isPairedEnd="$isPairedEnd"
fdr_cutoff=0.05

### Running DS analysis using diffSplice
LOG_03_2=$DIR/log/03_2_limma_diffSplice.log

export MODE="$MODE"
export DIR="$DIR" 
export REF="$REF" 
export target="$target" 
export pair="$pair"
export isPairedEnd="$isPairedEnd"
export fdr_cutoff="$fdr_cutoff"
export LOG_03_2="$LOG_03_2"

chmod +x "../DEU/03_2_limma_diffSplice.sh"

echo "====== `date`: limma:diffSplice ======" > $LOG_03_2
bash "../DEU/03_2_limma_diffSplice.sh"
echo "====== `date`: limma:diffSplice is done! ======" >> $LOG_03_2
