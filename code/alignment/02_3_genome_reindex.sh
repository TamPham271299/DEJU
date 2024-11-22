#!/bin/bash

module load STAR/2.7.11b

if [[ "$MODE" == "simulation" ]]; then
  
  mkdir -p $DIR/reindexed_genome/$i
  
  LOG=$DIR/reindexed_genome/$i/genome_reindexing.log
  echo "====== `date`: Re-indexing genome with filtered SJ ======" > $LOG
  CMD="STAR --runThreadN 16 \
            --runMode genomeGenerate \
            --genomeDir $DIR/reindexed_genome/$i \
            --genomeFastaFiles $FASTA \
            --sjdbOverhang $sjdbOverhang \
            --sjdbFileChrStartEnd $DIR/SJ/$i/merged_UMR_3.SJ.tab"
  echo $CMD >> $LOG
  eval $CMD >> $LOG
  echo "====== `date`: Genome re-indexing is done! ======" >> $LOG
  
  cat $LOG >> $LOG_02_3

elif [[ "$MODE" == "case_study" ]]; then
  
  echo "====== `date`: Re-indexing genome with filtered SJ ======" 
  CMD="STAR --runThreadN 16 \
            --runMode genomeGenerate \
            --genomeDir $DIR/reindexed_genome \
            --genomeFastaFiles $FASTA \
            --sjdbOverhang $sjdbOverhang \
            --sjdbFileChrStartEnd $DIR/SJ/merged_UMR_3.SJ.tab"
  echo $CMD
  eval $CMD
  echo "====== `date`: Genome re-indexing is done! ======"
  
else

  echo "Invalid mode specified. Choose 'simulation' or 'case_study'." >> $LOG_02_3
  
fi

