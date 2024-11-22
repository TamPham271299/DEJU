#!/bin/bash

sampleID_idx=$(head -1 $dl_list| tr '\t' '\n'| grep -n 'sample_title'| cut -d: -f1)

inputs_each=""
for SAMPLE in `sed '1d' $dl_list| cut -f$sampleID_idx`; do
    inputs_each=`echo -e $inputs_each$SAMPLE"\t"`
done

STAR_aligned_pass1 () {
    SAMPLE=`echo $1`
    LANE=`ls $TRIMMED/$SAMPLE`
    R1="$TRIMMED/$SAMPLE/$LANE/${LANE}.fastq.gz"
    
    mkdir -p $ALIGNED/$SAMPLE
    BAM="$ALIGNED/${SAMPLE}/${SAMPLE}."
    
    LOG=$(dirname $BAM)/STAR_aligned_pass1.log
    
    {
        echo "====== `date`: `echo $SAMPLE` ======"
        echo "====== `date`: Start 1-pass mapping ======" 
        CMD="STAR --genomeDir $DIR/reindexed_genome \
                    --readFilesIn $R1 \
                    --readFilesCommand zcat \
                    --outFileNamePrefix $BAM \
                    --outSAMtype BAM SortedByCoordinate \
                    --outFilterType BySJout \
                    --outFilterIntronMotifs RemoveNoncanonical \
                    --runThreadN 16"
        echo $CMD 
        eval $CMD 
        echo "====== `date`: 1-pass mapping is done! ======" 
    } 2>&1 | tee -a "$LOG"

    cat $LOG >> $LOG_02_4
}

export -f STAR_aligned_pass1
parallel -j 2 STAR_aligned_pass1 ::: $inputs_each
