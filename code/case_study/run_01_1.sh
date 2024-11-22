#!/bin/bash 
#SBATCH --job-name=01_1_QC
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --time=00:30:00
#SBATCH --output=slurm-%A.out

MODE="case_study"
# DIR="../../data/case_study/$projName/"
DIR="$DIR"
RAW="$DIR/raw/"
bioPrjID="$bioPrjID"
dl_list="$DIR/target/sampleInfo_${bioPrjID}.txt"
isPairedEnd="$isPairedEnd"

LOG_01_1="$DIR/log/01_1_QC.log"

export DIR="$DIR"
export RAW="$RAW"
export dl_list="$dl_list"
export MODE="$MODE"
export LOG_01_1="$LOG_01_1"

### Start running FASTQC
echo "====== `date`: QC ======" > $LOG_01_1

if [[ $isPairedEnd == TRUE ]]; then
  chmod +x "../alignment/01_1_QC.sh"
  bash "../alignment/01_1_QC.sh"
else
  chmod +x "../alignment/01_1_QC_se.sh"
  bash "../alignment/01_1_QC_se.sh"
fi

echo "====== `date`: QC is done! ======" >> $LOG_01_1
