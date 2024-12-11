#!/bin/bash 
#SBATCH --job-name=03_3_DEXSeq
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
REF="../../annotation/"
BAM="$DIR/aligned_pass2/"
flat_GTF="$REF/gencode.vM32.annotation.flat.DEXSeq.gtf"
target="../../data/simulation/target/target.${libs}vs${libs}.tsv"
pair="Group_2-Group_1"
MODE="simulation"
fdr_cutoff=0.05

seed=2024
workers=4

### Running DS analysis using DEXSeq
LOG_03_3=$DIR/log/03_3_DEXSeq.log

if [ ! -d "$DIR/log" ]; then
  mkdir -p $DIR/log
fi

export MODE="$MODE"
export DIR="$DIR" 
export REF="$REF"
export BAM="$BAM" 
export flat_GTF="$flat_GTF"
export target="$target"
export pair="$pair"
export fdr_cutoff="$fdr_cutoff"
export seed="$seed" 
export workers="$workers" 
export noOfSim="$noOfSim" 
export LOG_03_3="$LOG_03_3"

chmod +x "../DEU/03_3_DEXSeq.sh"

echo "====== `date`: DEXSeq ======" > $LOG_03_3
bash "../DEU/03_3_DEXSeq.sh"
echo "====== `date`: DEXSeq is done! ======" >> $LOG_03_3
