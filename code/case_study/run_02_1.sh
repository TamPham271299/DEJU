#!/bin/bash 
#SBATCH --job-name=02_1_STAR_aligned_pass1
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --time=2:30:00
#SBATCH --output=slurm-%A.out

MODE="case_study"
DIR="$DIR"
TRIMMED="$DIR/trimmed/"
bioPrjID="$bioPrjID"
dl_list="$DIR/target/sampleInfo_${bioPrjID}.txt"
maxRL=$maxRL
genomeDir="../../data/genomeIndex_${maxRL}bp"
isPairedEnd="$isPairedEnd"

LOG_02_1="$DIR/log/02_1_STAR_aligned_pass1.log"

export DIR="$DIR"
export TRIMMED="$TRIMMED"
export genomeDir="$genomeDir"
export dl_list="$dl_list"
export MODE="$MODE"
export LOG_02_1="$LOG_02_1"

### Start running STAR aligner with 1-pass mapping
echo "====== `date`: 1-pass mapping using STAR ======" > $LOG_02_1

if [ ! -d "$DIR/aligned_pass1" ]; then
  mkdir -p $DIR/aligned_pass1
fi

if [[ $isPairedEnd == TRUE ]]; then
  chmod +x "../alignment/02_1_STAR_aligned_pass1.sh"
  bash "../alignment/02_1_STAR_aligned_pass1.sh"
else
  chmod +x "../alignment/02_1_STAR_aligned_pass1_se.sh"
  bash "../alignment/02_1_STAR_aligned_pass1_se.sh"
fi

echo "====== `date`: 1-pass mapping using STAR is done! ======" >> $LOG_02_1
