#!/bin/bash 
#SBATCH --job-name=03_3_DEXSeq
#SBATCH --array=4
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --time=2:00:00
#SBATCH --output=slurm-%A.out

MODE="case_study"
DIR="$DIR"
BAM="$DIR/aligned_pass2/"
target="$DIR/target/target.tsv"
REF="../../annotation/"
flat_GTF="$REF/gencode.vM32.annotation.flat.DEXSeq.gtf"
pair="$pair"
fdr_cutoff=0.05

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
export LOG_03_3="$LOG_03_3"

chmod +x "../DEU/03_3_DEXSeq.sh"

echo "====== `date`: DEXSeq ======" > $LOG_03_3
bash "../DEU/03_3_DEXSeq.sh"
echo "====== `date`: DEXSeq is done! ======" >> $LOG_03_3
