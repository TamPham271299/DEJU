#!/bin/bash 
#SBATCH --job-name=03_4_JunctionSeq
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --time=2:00:00
#SBATCH --output=slurm-%A.out

MODE="case_study"
# DIR="../../data/case_study/$projName/"
DIR="$DIR"
BAM="$DIR/aligned_pass2/"
maxRL="$maxRL"
JAR_DIR="../../tools/QoRTs"
target_JunctionSeq="$DIR/target/decoder.bySample.txt"
REF="../../annotation/"
GTF="$REF/gencode.vM32.annotation.gtf.gz"
pair="$pair"
fdr_cutoff=0.05
ncores=4

### Running DS analysis using JunctionSeq
LOG_03_4=$DIR/log/03_4_JunctionSeq.log

export MODE="$MODE"
export DIR="$DIR"
export REF="$REF"
export BAM="$BAM"
export maxRL="$maxRL"
export JAR_DIR="$JAR_DIR" 
export GTF="$GTF" 
export target="$target_JunctionSeq" 
export pair="$pair"
export fdr_cutoff="$fdr_cutoff"
export ncores="$ncores"
export LOG_03_4="$LOG_03_4"

chmod +x "../DEU/03_4_JunctionSeq.sh"

echo "====== `date`: JunctionSeq ======" > $LOG_03_4
echo $pair >> $LOG_03_4
bash "../DEU/03_4_JunctionSeq.sh"
echo "====== `date`: JunctionSeq is done! ======" >> $LOG_03_4
