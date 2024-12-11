#!/bin/bash 
#SBATCH --job-name=02_2_3_SJfilter_genomeReindex
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
REF="../../annotation/"
FASTA="$REF/GRCm39.primary_assembly.genome.fa"
maxRL=75
sjdbOverhang=$((maxRL - 1))
MODE="simulation"

### Running STAR aligner with SJ filtering and genome reindexing
LOG_02_2=$DIR/log/02_2_SJ_filtering.log
LOG_02_3=$DIR/log/02_3_genome_reindex.log

echo "====== `date`: SJ filtering ======" > $LOG_02_2
echo "====== `date`: Genome re-indexing using filtered SJ ======" > $LOG_02_3

export DIR="$DIR"
export FASTA="$FASTA"
export sjdbOverhang="$sjdbOverhang"
export MODE="$MODE"
export LOG_02_2="$LOG_02_2"
export LOG_02_3="$LOG_02_3"

for i in `ls $RAW`; do

  export i="$i"
  
  if [ ! -d "$DIR/SJ/$i" ]; then
    mkdir -p $DIR/SJ/$i
  fi
  
  echo "====== `date`: Simulation `echo $i` ======" >> $LOG_02_2
  bash "../alignment/02_2_SJ_filtering.sh"
  
  if [ ! -d "$DIR/reindexed_genome/$i" ]; then
    mkdir -p $DIR/reindexed_genome/$i
  fi
  
  echo "====== `date`: Simulation `echo $i` ======" >> $LOG_02_3
  bash "../alignment/02_3_genome_reindex.sh"

done

echo "====== `date`: SJ filtering is done! ======" >> $LOG_02_2
echo "====== `date`: Genome re-indexing using filtered SJ is done! ======" >> $LOG_02_3
