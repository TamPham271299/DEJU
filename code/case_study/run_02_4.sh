#!/bin/bash 
#SBATCH --job-name=02_4_STAR_aligned_pass2
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --time=1:30:00
#SBATCH --output=slurm-%A.out

MODE="case_study"
# DIR="../../data/case_study/$projName/"
DIR="$DIR"
TRIMMED="$DIR/trimmed/"
bioPrjID="$bioPrjID"
dl_list="$DIR/target/sampleInfo_${bioPrjID}.txt"
isPairedEnd="$isPairedEnd"

### Running STAR aligner with 2-pass mapping
LOG_02_4=$DIR/log/02_4_STAR_aligned_pass2.log

export DIR="$DIR"
export TRIMMED="$TRIMMED"
export dl_list="$dl_list"
export MODE="$MODE"
export LOG_02_4="$LOG_02_4"

echo "====== `date`: 2-pass mapping using STAR ======" > $LOG_02_4
  
if [ ! -d "$DIR/aligned_pass2" ]; then
  mkdir -p $DIR/aligned_pass2
fi

if [[ $isPairedEnd == TRUE ]]; then
  chmod +x "../alignment/02_4_STAR_aligned_pass2.sh"
  bash "../alignment/02_4_STAR_aligned_pass2.sh"
else
  chmod +x "../alignment/02_4_STAR_aligned_pass2_se.sh"
  bash "../alignment/02_4_STAR_aligned_pass2_se.sh"
fi

echo "====== `date`: 2-pass mapping using STAR is done! ======" >> $LOG_02_4
