#!/bin/bash

module load parallel/20240722
module load STAR/2.7.11b

### STAR with 2-pass read alignment
if [[ "$MODE" == "simulation" ]]; then
  
  ### Running read mapping for 2 samples simultaneously 
  inputs_each=""
  for SAMPLE in `sed '1d' $target| cut -f1`; do
      inputs_each=`echo -e $inputs_each$SAMPLE"\t"`
  done
  
  STAR_aligned_pass2 () {
      SAMPLE=`echo $1`
      mkdir -p $DIR/aligned_pass2/$i/$SAMPLE
      echo $SAMPLE
      R1="$RAW/$i/${SAMPLE}_R1.fastq.gz"
      R2="$RAW/$i/${SAMPLE}_R2.fastq.gz"
      echo $R1 
      echo $R2
      BAM="$DIR/aligned_pass2/$i/${SAMPLE}/${SAMPLE}."
      echo $BAM
      
      LOG=$(dirname $BAM)/STAR_aligned_pass2.log
      echo "====== `date`: Start 2-pass mapping ======" > $LOG
      echo "====== `date`: `echo $SAMPLE` ======" >> $LOG
      CMD="STAR --genomeDir $DIR/reindexed_genome/$i \
                --readFilesIn $R1 $R2 \
                --readFilesCommand zcat \
                --outFileNamePrefix $BAM \
                --outSAMtype BAM SortedByCoordinate \
                --outFilterType BySJout \
                --runThreadN 16"
      echo $CMD >> $LOG
      eval $CMD >> $LOG
      echo "====== `date`: 2-pass mapping is done! ======" >> $LOG
      cat $LOG >> $LOG_02_4
  }
  
  export -f STAR_aligned_pass2
  parallel -j 2 STAR_aligned_pass2 ::: $inputs_each

elif [[ "$MODE" == "case_study" ]]; then
  
  sampleID_idx=$(head -1 $dl_list| tr '\t' '\n'| grep -n 'sample_title'| cut -d: -f1)
  
  inputs_each=""
  for SAMPLE in `sed '1d' $dl_list| cut -f$sampleID_idx`; do
      inputs_each=`echo -e $inputs_each$SAMPLE"\t"`
  done
  
  STAR_aligned_pass2 () {
      SAMPLE=`echo $1`
      mkdir -p $DIR/aligned_pass2/$SAMPLE      
      LANE=`ls $TRIMMED/$SAMPLE`
      R1="$TRIMMED/$SAMPLE/$LANE/${LANE}_1.fastq.gz"
      R2="$TRIMMED/$SAMPLE/$LANE/${LANE}_2.fastq.gz"
      
      BAM="$DIR/aligned_pass2/${SAMPLE}/${SAMPLE}."
      
      LOG=$(dirname $BAM)/STAR_aligned_pass2.log
      
      {
          echo "====== `date`: `echo $SAMPLE` ======"
          echo "====== `date`: Start 2-pass mapping ======" 
          CMD="STAR --genomeDir $DIR/reindexed_genome \
                    --readFilesIn $R1 $R2 \
                    --readFilesCommand zcat \
                    --outFileNamePrefix $BAM \
                    --outSAMtype BAM SortedByCoordinate \
                    --outFilterType BySJout \
                    --outFilterIntronMotifs RemoveNoncanonical \
                    --runThreadN 16"
          echo $CMD 
          eval $CMD 
          echo "====== `date`: 2-pass mapping is done! ======" 
      } 2>&1 | tee -a "$LOG"
  
      cat $LOG >> $LOG_02_4
  }
  
  export -f STAR_aligned_pass2
  parallel -j 2 STAR_aligned_pass2 ::: $inputs_each

else

  echo "Invalid mode specified. Choose 'simulation' or 'case_study'." >> $LOG_02_4
  
fi