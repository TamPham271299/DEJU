#!/bin/bash 
#SBATCH --job-name=01_2_read_trimming
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --time=1:30:00
#SBATCH --output=slurm-%A.out

MODE="case_study"
DIR="$DIR"
RAW="$DIR/raw/"
TRIMMED="$DIR/trimmed/"
bioPrjID="$bioPrjID"
dl_list="$DIR/target/sampleInfo_${bioPrjID}.txt"
isPairedEnd="$isPairedEnd"

LOG_01_2="$DIR/log/01_2_read_trimming.log"

export DIR="$DIR"
export RAW="$RAW"
export TRIMMED="$TRIMMED"
export dl_list="$dl_list"
export MODE="$MODE"
export LOG_01_2="$LOG_01_2"

### Start read trimming
echo "====== `date`: Read trimming ======" > $LOG_01_2

if [[ ! -d "$DIR/trimmed" ]]; then
  mkdir -p $DIR/trimmed
fi

if [[ $isPairedEnd == TRUE ]]; then
  chmod +x "../alignment/01_2_read_trimming.sh"
  bash "../alignment/01_2_read_trimming.sh"
else
  chmod +x "../alignment/01_2_read_trimming_se.sh"
  bash "../alignment/01_2_read_trimming_se.sh"
fi

echo "====== `date`: Read trimming is done! ======" >> $LOG_01_2
