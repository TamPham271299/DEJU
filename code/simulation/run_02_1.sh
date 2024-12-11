#!/bin/bash 
#SBATCH --job-name=02_1_STAR_aligned_pass1
#SBATCH --array=1-6
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00
#SBATCH --output=slurm-%A_%a.out

# Reading parameters
if [ "$sim_mode" == "DEU" ]; then
  par="parameters_DEU.txt"
elif [ "$sim_mode" == "null" ]; then
  par="parameters_null.txt"
fi
parVec=$(tail -n +2 $par| sed -n "${SLURM_ARRAY_TASK_ID}"p)
scenario=$(echo $parVec| awk '{print $1}')
libs=$(echo $parVec| awk '{print $2}')
rlen=$(echo $parVec| awk '{print $3}')
fc=$(echo $parVec| awk '{print $4}')

DIR="../../data/simulation/${scenario}_${libs}_${rlen}_${fc}/"
RAW="$DIR/raw/"
maxRL=75
sjdbOverhang=$((maxRL - 1)) # optimal sjdbOverhang
genomeDir="../../data/genomeIndex_${maxRL}bp"
target="../../data/simulation/target/target.${libs}vs${libs}.tsv"
MODE="simulation"

LOG_02_1="$DIR/log/02_1_STAR_aligned_pass1.log"

export DIR="$DIR"
export RAW="$RAW"
export genomeDir="$genomeDir"
export target="$target"
export MODE="$MODE"
export LOG_02_1="$LOG_02_1"

### Start running STAR aligner with 1-pass mapping
echo "====== `date`: 1-pass mapping using STAR ======" > $LOG_02_1

for i in `ls $RAW`; do

  export i="$i"
    
  if [ ! -d "$DIR/aligned_pass1/$i" ]; then
    mkdir -p $DIR/aligned_pass1/$i
  fi
  
  echo "====== `date`: Simulation `echo $i` ======" >> $LOG_02_1
  bash "../alignment/02_1_STAR_aligned_pass1.sh"
    
done

echo "====== `date`: 1-pass mapping using STAR is done! ======" >> $LOG_02_1
