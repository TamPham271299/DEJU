#!/bin/bash

module load fastqc/0.12.1
module load parallel/20240722

sampleID_idx=$(head -1 $dl_list| tr '\t' '\n'| grep -n 'sample_title'| cut -d: -f1)

inputs_each=""
for SAMPLE in `sed '1d' $dl_list| cut -f$sampleID_idx`; do
    inputs_each=`echo -e $inputs_each$SAMPLE"\t"`
done

QC () {
    SAMPLE=`echo $1`
    LANE=`ls $RAW/$SAMPLE`
    R1="$RAW/$SAMPLE/$LANE/${LANE}_1.fastq.gz"
    R2="$RAW/$SAMPLE/$LANE/${LANE}_2.fastq.gz"
    echo $R1 
    echo $R2
    echo $(dirname $R1)
    
    LOG=$(dirname $R1)/QC.log
    
    {
        echo "====== `date`: `echo $SAMPLE` ======" 
        CMD="fastqc --nogroup -t 16 --dir $(dirname $R1) -o $(dirname $R1) $R1 $R2"
        echo $CMD
        eval $CMD
        echo "====== `date`: Quality control is done! ======" 
        
    } 2>&1 | tee -a "$LOG"
    
    cat $LOG >> $LOG_01_1
}

export -f QC
parallel -j 8 QC ::: $inputs_each
