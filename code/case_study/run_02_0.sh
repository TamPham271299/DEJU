#!/bin/bash 
#SBATCH --job-name=02_0_genome_index
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --time=2:00:00
#SBATCH --output=slurm-%A.out

module load parallel/20240722
module load STAR/2.7.11b

maxRL="$maxRL"
sjdbOverhang=$((maxRL - 1)) # optimal sjdbOverhang
genomeDir="../../data/genomeIndex_${maxRL}bp"
REF="../../annotation/"
FASTA="$REF/GRCm39.primary_assembly.genome.fa"
GTF="$REF/gencode.vM32.primary_assembly.annotation.gtf"

LOG_02_0="../../data/02_0_STAR_index_genome_${maxRL}bp.log"

export genomeDir="$genomeDir"
export sjdbOverhang="$sjdbOverhang"
export FASTA="$FASTA"
export GTF="$GTF"
export LOG_02_0="$LOG_02_0"

chmod +x "../alignment/02_0_genome_index.sh"
bash "../alignment/02_0_genome_index.sh"