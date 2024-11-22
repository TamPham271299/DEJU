#!/bin/bash 
#SBATCH --job-name=03_4_JunctionSeq
#SBATCH --array=4
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --time=2:00:00
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
# REF="/vast/projects/lab_chen/tam/ref_genome/Mus_musculus/Gencode/"
# GTF="$REF/gencode.vM32.annotation.gtf.gz"
# JAR_DIR="/vast/projects/lab_chen/tam/tools/QoRTs"
# target_JunctionSeq="/vast/projects/lab_chen/tam/Differential_splicing/simulation_3vs3_mouse_genome/RNA-seq/config/decoder.bySample.${libs}vs${libs}.txt"

DIR="../../data/simulation/${scenario}_${libs}_${rlen}_${fc}/"
BAM="$DIR/aligned_pass2/"
REF="../../annotation/"
GTF="$REF/gencode.vM32.annotation.gtf.gz"
JAR_DIR="../../tools/QoRTs"
target_JunctionSeq="../../data/simulation/target/decoder.bySample.${libs}vs${libs}.txt"
pair="Group_2-Group_1"
MODE="simulation"
fdr_cutoff=0.05

# noOfSim=20
seed=2024
workers=4
ncores=4

### Running DS analysis using JunctionSeq
LOG_03_4=$DIR/log/03_4_JunctionSeq.log

export MODE="$MODE"
export DIR="$DIR" 
export REF="$REF"
export BAM="$BAM" 
export JAR_DIR="$JAR_DIR" 
export GTF="$GTF" 
export target="$target_JunctionSeq" 
export pair="$pair"
export fdr_cutoff="$fdr_cutoff"
export noOfSim="$noOfSim" 
export ncores="$ncores"
export LOG_03_4="$LOG_03_4"

chmod +x "../DEU/03_4_JunctionSeq.sh"

echo "====== `date`: JunctionSeq ======" >> $LOG_03_4
echo $pair >> $LOG_03_4
bash "../DEU/03_4_JunctionSeq.sh"
echo "====== `date`: JunctionSeq is done! ======" >> $LOG_03_4
