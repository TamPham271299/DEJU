#!/bin/bash 
#SBATCH --job-name=02_2_3_SJfilter_genomeReindex
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --time=2:00:00
#SBATCH --output=slurm-%A.out

MODE="case_study"
DIR="$DIR"
REF="../../annotation/"
FASTA="$REF/GRCm39.primary_assembly.genome.fa"
maxRL=$maxRL
sjdbOverhang=$((maxRL - 1))
            
### Running STAR aligner with SJ filtering and genome reindexing
LOG_02_2=$DIR/log/02_2_SJ_filtering.log
LOG_02_3=$DIR/log/02_3_genome_reindex.log

export DIR="$DIR"
export FASTA="$FASTA"
export sjdbOverhang="$sjdbOverhang"
export MODE="$MODE"
export LOG_02_2="$LOG_02_2"
export LOG_02_3="$LOG_02_3"

echo "====== `date`: SJ filtering ======" > $LOG_02_2
if [ ! -d "$DIR/SJ" ]; then
  mkdir -p $DIR/SJ
fi
chmod +x "../alignment/02_2_SJ_filtering.sh"
bash "../alignment/02_2_SJ_filtering.sh" >> $LOG_02_2 2>&1
echo "====== `date`: SJ filtering is done! ======" >> $LOG_02_2

echo "====== `date`: Genome re-indexing using filtered SJ ======" > $LOG_02_3
if [ ! -d "$DIR/reindexed_genome" ]; then
  mkdir -p $DIR/reindexed_genome
fi
chmod +x "../alignment/02_3_genome_reindex.sh"
bash "../alignment/02_3_genome_reindex.sh" >> $LOG_02_3 2>&1
echo "====== `date`: Genome re-indexing using filtered SJ is done! ======" >> $LOG_02_3
