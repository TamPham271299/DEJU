#!/bin/bash 

module unload openjdk
module load R/4.3.3

if [ "$MODE" == "simulation" ]; then

  inputs_each=""
  for simID in `ls $BAM`; do
      inputs_each=`echo -e $inputs_each$simID"\t"`
  done

  QoRTs_analysis () {
      i=`echo $1`

      ### Generating raw counts via QoRTs
      OUTPUT="$DIR/QoRTs/$i"
      mkdir -p $OUTPUT

      LOG="$OUTPUT/QoRTs.log"

      {
          echo "====== `date`: Running QoRTs ======"
          echo "====== `date`: Processing $i ======"
          for SAMPLE in `ls $BAM/$i`; do
              echo "====== `date`: Start counting reads using QoRTs for $SAMPLE ======"
              mkdir -p $OUTPUT/$SAMPLE
              bam="$BAM/$i/$SAMPLE/${SAMPLE}.Aligned.sortedByCoord.out.bam"
              echo $bam
              CMD1="java -jar $JAR_DIR/QoRTs.jar QC --runFunctions writeKnownSplices,writeNovelSplices,writeSpliceExon $bam $GTF $OUTPUT/$SAMPLE"
              echo $CMD1
              eval $CMD1
          done

          ### Including Novel Splice Junction Loci
          echo "====== `date`: Start counting novel junction reads ======"
          CMD2="java -jar $JAR_DIR/QoRTs.jar mergeNovelSplices --minCount 6 $OUTPUT $target $GTF $OUTPUT"
          echo $CMD2
          eval $CMD2

      } 2>&1 | tee -a "$LOG"

      cat $LOG >> $LOG_03_4
  }

  export -f QoRTs_analysis
  parallel -j 2 QoRTs_analysis ::: $inputs_each
  
  echo $pair >> $LOG_03_4
  echo "====== `date`: Start running JunctionSeq ======" >> $LOG_03_4
  Rscript --no-save --no-restore --verbose ../DEU/03_4_JunctionSeq.R \
            --MODE=$MODE \
            --DIR=$DIR \
            --REF=$REF \
            --target=$target \
            --pair=$pair \
            --fdr_cutoff=$fdr_cutoff \
            --ncores=$ncores \
            --noOfSim=$noOfSim >> $LOG_03_4 2>&1

elif [ "$MODE" == "case_study" ]; then
  
  ### Generating raw counts via QoRTs
  OUTPUT="$DIR/QoRTs/"
  mkdir -p $OUTPUT
  
  for SAMPLE in `ls $BAM`; do
      echo $SAMPLE

      if [[ ! -e "$DIR/QoRTs/$SAMPLE/QC.QORTS_COMPLETED_OK" ]]; then
        mkdir -p $OUTPUT/$SAMPLE
        bam="$BAM/$SAMPLE/${SAMPLE}.Aligned.sortedByCoord.out.bam"
        echo $bam
        LOG="$OUTPUT/$SAMPLE/log.txt"
        echo "====== `date`: Start counting reads using QoRTs for $SAMPLE ======" > $LOG
        CMD="java -jar $JAR_DIR/QoRTs.jar QC --maxReadLength $maxRL --runFunctions writeKnownSplices,writeNovelSplices,writeSpliceExon $bam $GTF $OUTPUT/$SAMPLE"
        echo $CMD >> $LOG
        eval $CMD >> $LOG
        echo "====== `date`: Counting reads is done ======" >> $LOG
        cat $LOG >> $LOG_03_4
      fi

  done

  if [[ ! -e "$DIR/QoRTs/withNovel.forJunctionSeq.gff.gz" ]]; then
  
    echo "====== `date`: Including Novel Splice Junction Loci ======" >> $LOG_03_4
    CMD2="java -jar $JAR_DIR/QoRTs.jar mergeNovelSplices \
                                --minCount 6 \
                                $OUTPUT \
                                $target \
                                $GTF \
                                $OUTPUT"
    echo $CMD2 >> $LOG_03_4
    eval $CMD2 >> $LOG_03_4
    echo "====== `date`: Including novel SJ is done ======" >> $LOG_03_4
  
  fi

  echo "====== `date`: Start running JunctionSeq ======" >> $LOG_03_4
  IFS=',' read -r -a pair <<< "$pair"
    
  for p in "${pair[@]}"; do
    echo "====== $p ======" >> $LOG_03_4
    if [[ ! -e "$DIR/target/decoder.bySample.${p}.txt" ]]; then
      IFS='-' read -r g1 g2 <<< "$p"
      awk -v g1="$g1" -v g2="$g2" 'NR==1 || $2 == g1 || $2 == g2' $target > "$DIR/target/decoder.bySample.${p}.txt"
    fi
    Rscript --no-save --no-restore --verbose ../DEU/03_4_JunctionSeq.R \
              --MODE=$MODE \
              --DIR=$DIR \
              --REF=$REF \
              --target="$DIR/target/decoder.bySample.${p}.txt" \
              --pair=$p \
              --fdr_cutoff=$fdr_cutoff \
              --ncores=$ncores >> $LOG_03_4 2>&1
  done
  
else

  echo "Invalid mode specified. Choose 'simulation' or 'case_study'." >> $LOG_03_4  
  
fi

