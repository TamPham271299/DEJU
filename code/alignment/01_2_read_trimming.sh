#!/bin/bash

module load fastqc/0.12.1
module load trimgalore/0.6.10
module load parallel/20240722

sampleID_idx=$(head -1 $dl_list| tr '\t' '\n'| grep -n 'sample_title'| cut -d: -f1)

inputs_each=""
for SAMPLE in `sed '1d' $dl_list| cut -f$sampleID_idx| tail -n +5`; do
    inputs_each=`echo -e $inputs_each$SAMPLE"\t"`
done

read_trimming () {
    SAMPLE=`echo $1`
    LANE=`ls $RAW/$SAMPLE`
    R1="$RAW/$SAMPLE/$LANE/${LANE}_1.fastq.gz"
    R2="$RAW/$SAMPLE/$LANE/${LANE}_2.fastq.gz"
    
    OUTPUT="$TRIMMED/$SAMPLE/$LANE"
    mkdir -p $OUTPUT
    
    LOG=$OUTPUT/read_trimming.log
    
    {
        echo "====== `date`: `echo $SAMPLE` ======"
        echo "====== `date`: Read trimming using trim-galore ======"
        CMD1="trim_galore --cores 4 -q 30 --length 20 --paired --fastqc_args \"--nogroup\" --output_dir $OUTPUT $R1 $R2"
        echo $CMD1
        eval $CMD1
        echo "====== `date`: Trimming is done! ======"
        
        echo "====== `date`: `echo $SAMPLE` ======"
        echo "====== `date`: Rename fastq files after trimming ======"
        CMD2='for i in `find $OUTPUT -name "*_val_*.fq.gz"`; do mv "$i" "${i/_val_?.fq.gz/.fastq.gz}"; done'
        echo "$CMD2"
        eval "$CMD2"
        echo "====== `date`: Rename is done! ======"
        
    } 2>&1 | tee -a "$LOG"
    
    cat $LOG >> $LOG_01_2
}


export -f read_trimming
parallel -j 4 read_trimming ::: $inputs_each
