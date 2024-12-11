#!/bin/bash 
#SBATCH --job-name=03_2_limma_diffSplice
#SBATCH --array=1-6
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --time=12:00:00
#SBATCH --output=slurm-%A_%a.out

module load parallel/20240722
module load R/4.3.3 # also load openjdk/21.0.2

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

MODE="simulation"
DIR="../../data/simulation/${scenario}_${libs}_${rlen}_${fc}/"
REF="../../annotation/"
target="../../data/simulation/target/target.${libs}vs${libs}.tsv"
pair="Group_2-Group_1"
isPairedEnd=TRUE
fdr_cutoff=0.05

seed=2024
workers=4

### Running DS analysis using diffSplice
LOG_03_2=$DIR/log/03_2_limma_diffSplice.log

export MODE="$MODE"
export DIR="$DIR" 
export REF="$REF" 
export target="$target" 
export pair="$pair"
export isPairedEnd="$isPairedEnd"
export fdr_cutoff="$fdr_cutoff"
export seed="$seed" 
export workers="$workers" 
export noOfSim="$noOfSim" 
export LOG_03_2="$LOG_03_2"

chmod +x "../DEU/03_2_limma_diffSplice.sh"

echo "====== `date`: limma:diffSplice ======" > $LOG_03_2
bash "../DEU/03_2_limma_diffSplice.sh"
echo "====== `date`: limma:diffSplice is done! ======" >> $LOG_03_2
