#!/bin/bash 
#SBATCH --job-name=00_1_fastq_dl_sra
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --time=2:30:00
#SBATCH --output=slurm-%A.out

DIR="$DIR"
RAW="$DIR/raw/"
bioPrjID="$bioPrjID"
dl_list="$DIR/config/sampleInfo_${bioPrjID}.txt"
isPairedEnd="$isPairedEnd"

LOG_00_1="$DIR/log/00_1_fastq_dl_sra.log"

export DIR="$DIR"
export RAW="$RAW"
export dl_list="$dl_list"
export LOG_00_1="$LOG_00_1"

### Start running fastq_dump
echo "====== `date`: FASTQ download using fastq_dump ======" > $LOG_00_1

if [[ $isPairedEnd == TRUE ]]; then
  chmod +x "../fastq_dl/00_1_sra_dl_pe.sh"
  bash "../fastq_dl/00_1_sra_dl_pe.sh"
else
  chmod +x "../fastq_dl/00_1_sra_dl_se.sh"
  bash "../fastq_dl/00_1_sra_dl_se.sh"
fi

echo "====== `date`: Downloading FASTQ is done! ======" >> $LOG_01_1
