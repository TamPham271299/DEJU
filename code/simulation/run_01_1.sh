#!/bin/bash 
#SBATCH --job-name=01_1_simulation
#SBATCH --array=3,4
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --time=2:00:00
#SBATCH --output=slurm-%A_%a.out

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

### Path to scenario
# DIR="/vast/projects/MM/tam/Differential_splicing/simulation_3vs3_mouse_genome/RNA-seq/DEU_mix/${scenario}_${libs}_${rlen}_${fc}/"
DIR="../../data/simulation/${scenario}_${libs}_${rlen}_${fc}/"
# scripts="/vast/projects/Spatial/tam/Differential_splicing/DEJU/code/simulation"

### Customized transcriptome
# transcriptome_fasta="/vast/projects/lab_chen/tam/ref_genome/Mus_musculus/Gencode/PC_genes_ver2/AS/gencode.vM32.clean.transcriptome.1.fa"
transcriptome_fasta="../../data/simulation/customized_transcriptome/gencode.vM32.custom.transcriptome.fa"
# metadata_of_transcriptome="/vast/projects/lab_chen/tam/ref_genome/Mus_musculus/Gencode/PC_genes_ver2/AS/AS_genes.info.tsv"
metadata_of_transcriptome="../../data/simulation/customized_transcriptome/DEU_genes.info.tsv"

### simReads settings
pe=TRUE
simulate_sequencing_error=FALSE
quality_reference="NA"

### DEU simulation settings
number_of_DEUs_for_ES_pattern=250
number_of_DEUs_for_MXE_pattern=250
number_of_DEUs_for_ASS_pattern=250
number_of_DEUs_for_RI_pattern=250

### general settings
# noOfSim=20
seed=2024
workers=20 # max 96

### Generate .log file
LOG_01_1=$DIR/log/01_1_simulation.log

### Start running simulation
echo "====== `date`: Start simulating RNA-Seq reads ======" > $LOG_01_1

Rscript --no-save --no-restore --verbose main.R \
  --path_to_project=$DIR \
  --transcriptome_fasta=$transcriptome_fasta \
  --metadata_of_transcriptome=$metadata_of_transcriptome \
  --libs=$libs \
  --scenario=$scenario \
  --rlen=$rlen \
  --pe=$pe \
  --simulate_sequencing_error=$simulate_sequencing_error \
  --quality_reference=$quality_reference \
  --fc=$fc \
  --number_of_DEUs_for_ES_pattern=$number_of_DEUs_for_ES_pattern \
  --number_of_DEUs_for_MXE_pattern=$number_of_DEUs_for_MXE_pattern \
  --number_of_DEUs_for_ASS_pattern=$number_of_DEUs_for_ASS_pattern \
  --number_of_DEUs_for_RI_pattern=$number_of_DEUs_for_RI_pattern \
  --number_of_simulation=$noOfSim \
  --seed=$seed \
  --workers=$workers >> $LOG_01_1 2>&1

echo "====== `date`: Simulating RNA-Seq reads is done! ======" >> $LOG_01_1
