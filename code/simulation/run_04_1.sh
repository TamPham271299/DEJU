#!/bin/bash 
#SBATCH --job-name=04_1_performance_analysis
#SBATCH --array=1-6
#SBATCH --mem=4G
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
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

DIR="../../data/simulation/${scenario}_${libs}_${rlen}_${fc}/"
DEU_map="../../data/simulation/customized_transcriptome/DEU_genes.info.tsv"
fdr_cutoff=0.05

seed=2024
workers=4

### Performance analysis
LOG_04_1=$DIR/log/04_1_simulation_analysis.log

export DIR="$DIR"
export DEU_map="$DEU_map" 
export fdr_cutoff="$fdr_cutoff"
export seed="$seed" 
export workers="$workers" 
export noOfSim="$noOfSim"

echo "====== `date`: Benchmarking ======" > $LOG_04_1
chmod +x "../DEU/04_1_performance_analysis.sh"
bash "../DEU/04_1_performance_analysis.sh" >> $LOG_04_1 2>&1
echo "====== `date`: Benchmarking is done! ======" >> $LOG_04_1
