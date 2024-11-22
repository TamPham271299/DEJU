#!/bin/bash 
#SBATCH --job-name=02_4_STAR_aligned_pass2
#SBATCH --array=3,4
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --time=3:00:00
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

# DIR="/vast/projects/MM/tam/Differential_splicing/simulation_3vs3_mouse_genome/RNA-seq/DEU_mix/${scenario}_${libs}_${rlen}_${fc}/"
# RAW="$DIR/raw_fit_RL_75/"
# target="/vast/projects/lab_chen/tam/Differential_splicing/simulation_3vs3_mouse_genome/RNA-seq/config/target.${libs}vs${libs}.tsv"

DIR="../../data/simulation/${scenario}_${libs}_${rlen}_${fc}/"
RAW="$DIR/raw/"
target="../../data/simulation/target/target.${libs}vs${libs}.tsv"
MODE="simulation"

### Running STAR aligner with 2-pass mapping
LOG_02_4=$DIR/log/02_4_STAR_aligned_pass2.log

export DIR="$DIR"
export RAW="$RAW"
export target="$target"
export MODE="$MODE"
export LOG_02_4="$LOG_02_4"

echo "====== `date`: 2-pass mapping using STAR ======" > $LOG_02_4

for i in `ls $RAW`; do

  export i="$i"
  
  if [ ! -d "$DIR/aligned_pass2/$i" ]; then
    mkdir -p $DIR/aligned_pass2/$i
  fi
  
  echo "====== `date`: Simulation $i ======" >> $LOG_02_4
  bash "../alignment/02_4_STAR_aligned_pass2.sh"
    
done

echo "====== `date`: 2-pass mapping using STAR is done! ======" >> $LOG_02_4
