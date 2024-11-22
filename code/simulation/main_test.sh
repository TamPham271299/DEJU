#!/bin/bash 
#SBATCH --job-name=test
#SBATCH --array=1-6
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --time=12:00:00
#SBATCH --output=slurm-%A_%a.out

module load STAR/2.7.11b
module load parallel/20240722
module load subread/2.0.6
module load R/4.3.3 # also load openjdk/21.0.2

# Reading parameters
parVec=$(tail -n +2 parameters_DEU.txt| sed -n "${SLURM_ARRAY_TASK_ID}"p)
scenario=$(echo $parVec| awk '{print $1}')
libs=$(echo $parVec| awk '{print $2}')
rlen=$(echo $parVec| awk '{print $3}')
fc=$(echo $parVec| awk '{print $4}')

### Parameter settings
# DIR="/vast/projects/Spatial/tam/Differential_splicing/simulation_3vs3_mouse_genome/RNA-seq/DEU_mix/${scenario}_${libs}_${rlen}_${fc}/"
# DIR="/vast/projects/MM/tam/Differential_splicing/simulation_3vs3_mouse_genome/RNA-seq/DEU_mix/${scenario}_${libs}_${rlen}_${fc}/"
# 
# RAW="$DIR/raw_fit_RL_75/"
# BAM="$DIR/STAR_aligned_pass2_minUniqSJReads_3/"
# 
# REF="/vast/projects/lab_chen/tam/ref_genome/Mus_musculus/Gencode/"
# FASTA="$REF/GRCm39.primary_assembly.genome.fa"
# GTF="$REF/gencode.vM32.annotation.gtf.gz"
# flat_GTF="$REF/gencode.vM32.annotation.flat.DEXSeq.gtf"
# genomeDir="/vast/projects/lab_chen/tam/Differential_splicing/simulation_3vs3_mouse_genome/RNA-seq/STAR_2_7_9a_genomeIndex/pass1_75bp"
# JAR_DIR="/vast/projects/lab_chen/tam/tools/QoRTs"
# target_JunctionSeq="/vast/projects/lab_chen/tam/Differential_splicing/simulation_3vs3_mouse_genome/RNA-seq/config/decoder.bySample.${libs}vs${libs}.txt"
# target_diffSplice="/vast/projects/lab_chen/tam/Differential_splicing/simulation_3vs3_mouse_genome/RNA-seq/config/target.${libs}vs${libs}.tsv"
# scripts="/vast/projects/lab_chen/tam/Differential_splicing/DEJU/code/test"
# AS_map="/vast/projects/lab_chen/tam/ref_genome/Mus_musculus/Gencode/PC_genes_ver2/AS/AS_genes.info.tsv"

DIR="$currentDir/../../../data/simulation/${scenario}_${libs}_${rlen}_${fc}/"
RAW="$DIR/raw_fit_RL_75/"
BAM="$DIR/STAR_aligned_pass2_minUniqSJReads_3/"

REF="#currentDir/../../../annotations/"
FASTA="$REF/GRCm39.primary_assembly.genome.fa"
GTF="$REF/gencode.vM32.annotation.gtf.gz"
flat_GTF="$REF/gencode.vM32.annotation.flat.DEXSeq.gtf"
genomeDir="$currentDir/../../../data/simulation/STAR_index_75bp/"
JAR_DIR="$currentDir/../../../tools/QoRTs"
target_JunctionSeq="$currentDir/../../../data/simulation/config/decoder.bySample.${libs}vs${libs}.txt"
target="$currentDir/../../../data/simulation/config/target.${libs}vs${libs}.tsv"

mkdir -p $DIR/log
LOG_03_1=$DIR/log/03_1_DS_edgeR_diffSpliceDGE_test.log
LOG_03_2=$DIR/log/03_2_DS_limma_diffSplice_test.log
LOG_03_3=$DIR/log/03_3_DS_DEXSeq_test.log
LOG_03_4=$DIR/log/03_4_DS_JunctionSeq_test.log
cd $DIR

### Serial
ncores=1
workers=1

### 031. edgeR:diffSpliceDGE
echo "====== `date`: edgeR:diffSpliceDGE ======" >> $LOG_03_1
echo "====== `date`: Processing serial DEU ======" >> $LOG_03_1
/usr/bin/time -v Rscript --no-save --no-restore --verbose 03_1_edgeR_diffSpliceDGE.R \
    --DIR=$DIR \
    --REF=$REF \
    --target=$target >> $LOG_03_1 2>&1
echo "====== `date`: edgeR:diffSplice is done! ======" >> $LOG_03_1

### 032. limma:diffSplice
echo "====== `date`: limma:diffSplice ======" >> $LOG_03_2
echo "====== `date`: Processing serial DEU ======" >> $LOG_03_2
/usr/bin/time -v Rscript --no-save --no-restore --verbose 03_2_limma_diffSplice.R \
    --DIR=$DIR \
    --REF=$REF \
    --target=$target >> $LOG_03_2 2>&1
echo "====== `date`: limma:diffSplice is done! ======" >> $LOG_03_2

### 033. DEXSeq
echo "====== `date`: DEXSeq ======" >> $LOG_03_3
echo "====== `date`: Processing serial DEU ======" >> $LOG_03_3
/usr/bin/time -v Rscript --no-save --no-restore --verbose 03_3_DEXSeq.R \
    --DIR=$DIR \
    --flat_GTF=$flat_GTF \
    --libs=$libs \
    --workers=$workers >> $LOG_03_3 2>&1
echo "====== `date`: DEXSeq is done! ======" >> $LOG_03_3

### 034. JunctionSeq
echo "====== `date`: JunctionSeq ======" >> $LOG_03_4
echo "====== `date`: Processing serial DEU ======" >> $LOG_03_4
/usr/bin/time -v Rscript --no-save --no-restore --verbose 03_4_JunctionSeq.R \
    --DIR=$DIR \
    --target=$target_JunctionSeq \
    --ncores=$ncores >> $LOG_03_4 2>&1
echo "====== `date`: JunctionSeq is done! ======" >> $LOG_03_4

### Parallel
ncores=4
workers=4

### 033. DEXSeq
echo "====== `date`: DEXSeq ======" >> $LOG_03_3
echo "====== `date`: Processing parallel DEU ======" >> $LOG_03_3
/usr/bin/time -v Rscript --no-save --no-restore --verbose 03_3_DEXSeq.R \
    --DIR=$DIR \
    --flat_GTF=$flat_GTF \
    --libs=$libs \
    --workers=$workers >> $LOG_03_3 2>&1
echo "====== `date`: DEXSeq is done! ======" >> $LOG_03_3

### 034. JunctionSeq
echo "====== `date`: JunctionSeq ======" >> $LOG_03_4
echo "====== `date`: Processing parallel DEU ======" >> $LOG_03_4
/usr/bin/time -v Rscript --no-save --no-restore --verbose 03_4_JunctionSeq.R \
    --DIR=$DIR \
    --target=$target_JunctionSeq \
    --ncores=$ncores >> $LOG_03_4 2>&1
echo "====== `date`: JunctionSeq is done! ======" >> $LOG_03_4
