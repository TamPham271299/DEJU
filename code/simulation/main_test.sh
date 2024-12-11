#!/bin/bash 
#SBATCH --job-name=test
#SBATCH --array=3,6
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --time=5:00:00
#SBATCH --output=slurm-%A_%a.out

module load R/4.3.3 # also load openjdk/21.0.2

# Reading parameters
parVec=$(tail -n +2 parameters_DEU.txt| sed -n "${SLURM_ARRAY_TASK_ID}"p)
scenario=$(echo $parVec| awk '{print $1}')
libs=$(echo $parVec| awk '{print $2}')
rlen=$(echo $parVec| awk '{print $3}')
fc=$(echo $parVec| awk '{print $4}')

DIR="/vast/projects/MM/tam/Differential_splicing/simulation_3vs3_mouse_genome/RNA-seq/DEU_mix/${scenario}_${libs}_${rlen}_${fc}/"
# DIR="../../data/simulation/${scenario}_${libs}_${rlen}_${fc}/"

REF="../../annotation/"
FASTA="$REF/GRCm39.primary_assembly.genome.fa"
GTF="$REF/gencode.vM32.annotation.gtf.gz"
flat_GTF="$REF/gencode.vM32.annotation.flat.DEXSeq.gtf"
target_JunctionSeq="../../data/simulation/target/decoder.bySample.${libs}vs${libs}.txt"
target="../../data/simulation/target/target.${libs}vs${libs}.tsv"
pair="Group_2-Group_1"

mkdir -p $DIR/log
LOG_03_1=$DIR/log/03_1_DS_edgeR_diffSpliceDGE_test.log
LOG_03_2=$DIR/log/03_2_DS_limma_diffSplice_test.log
LOG_03_3=$DIR/log/03_3_DS_DEXSeq_test.log
LOG_03_4=$DIR/log/03_4_DS_JunctionSeq_test.log

### Serial
ncores=1
workers=1

### 031. edgeR:diffSpliceDGE
echo "====== `date`: edgeR:diffSpliceDGE ======" > $LOG_03_1
echo "====== `date`: Processing serial DEU ======" >> $LOG_03_1
/usr/bin/time -v Rscript --no-save --no-restore --verbose ../DEU/03_1_edgeR_diffSpliceDGE_test.R \
    --DIR=$DIR \
    --REF=$REF \
    --target=$target \
    --pair=$pair >> $LOG_03_1 2>&1
echo "====== `date`: edgeR:diffSplice is done! ======" >> $LOG_03_1

### 032. limma:diffSplice
echo "====== `date`: limma:diffSplice ======" > $LOG_03_2
echo "====== `date`: Processing serial DEU ======" >> $LOG_03_2
/usr/bin/time -v Rscript --no-save --no-restore --verbose ../DEU/03_2_limma_diffSplice_test.R \
    --DIR=$DIR \
    --REF=$REF \
    --target=$target \
    --pair=$pair >> $LOG_03_2 2>&1
echo "====== `date`: limma:diffSplice is done! ======" >> $LOG_03_2

### 033. DEXSeq
echo "====== `date`: DEXSeq ======" > $LOG_03_3
echo "====== `date`: Processing serial DEU ======" >> $LOG_03_3
/usr/bin/time -v Rscript --no-save --no-restore --verbose ../DEU/03_3_DEXSeq_test.R \
    --DIR=$DIR \
    --REF=$REF \
    --flat_GTF=$flat_GTF \
    --target=$target \
    --pair=$pair \
    --workers=$workers >> $LOG_03_3 2>&1
echo "====== `date`: DEXSeq is done! ======" >> $LOG_03_3

### 034. JunctionSeq
echo "====== `date`: JunctionSeq ======" > $LOG_03_4
echo "====== `date`: Processing serial DEU ======" >> $LOG_03_4
/usr/bin/time -v Rscript --no-save --no-restore --verbose ../DEU/03_4_JunctionSeq_test.R \
    --DIR=$DIR \
    --REF=$REF \
    --target=$target_JunctionSeq \
    --ncores=$ncores \
    --pair=$pair >> $LOG_03_4 2>&1
echo "====== `date`: JunctionSeq is done! ======" >> $LOG_03_4

### Parallel
ncores=4
workers=4

### 033. DEXSeq
echo "====== `date`: DEXSeq ======" >> $LOG_03_3
echo "====== `date`: Processing parallel DEU ======" >> $LOG_03_3
/usr/bin/time -v Rscript --no-save --no-restore --verbose ../DEU/03_3_DEXSeq_test.R \
    --DIR=$DIR \
    --REF=$REF \
    --flat_GTF=$flat_GTF \
    --target=$target \
    --pair=$pair \
    --workers=$workers >> $LOG_03_3 2>&1
echo "====== `date`: DEXSeq is done! ======" >> $LOG_03_3

### 034. JunctionSeq
echo "====== `date`: JunctionSeq ======" >> $LOG_03_4
echo "====== `date`: Processing parallel DEU ======" >> $LOG_03_4
/usr/bin/time -v Rscript --no-save --no-restore --verbose ../DEU/03_4_JunctionSeq_test.R \
    --DIR=$DIR \
    --REF=$REF \
    --target=$target_JunctionSeq \
    --ncores=$ncores \
    --pair=$pair >> $LOG_03_4 2>&1
echo "====== `date`: JunctionSeq is done! ======" >> $LOG_03_4
