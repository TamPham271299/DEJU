#!/bin/bash 

module load subread/2.0.6
module load parallel/20240722
module unload openjdk
module load R/4.3.3 # also load openjdk/21.0.2

if [ "$MODE" == "simulation" ]; then

  echo "====== `date`: Start quantifying exon-level reads using subread ======" >> $LOG_03_3
  inputs_each=""
  for simID in `ls $BAM`; do
      inputs_each=`echo -e $inputs_each$simID"\t"`
  done
  
  featureCounts_subread () {
    i=`echo $1`
    OUTPUT="$DIR/featureCounts_subread/$i"
    mkdir -p $OUTPUT
    
    if [[ ! -e "$OUTPUT/featureCounts_output.txt" ]]; then
    
      LOG="$OUTPUT/featureCounts_subread.log"
  
      bam_files=""
      for sampleID in `sed '1d' $target| cut -f1`; do
          bam_file="$BAM/$i/$sampleID/${sampleID}.Aligned.sortedByCoord.out.bam"
          echo $bam_file
          bam_files="$bam_files $bam_file"
          echo $bam_files
      done
  
      {
        echo "====== `date`: Processing $i ======"
        
        # Run featureCounts
        CMD="featureCounts -p \
                    --countReadPairs \
                    -M \
                    -f \
                    -O \
                    -T 8 \
                    -Q 255 \
                    -F 'GTF' \
                    -a $flat_GTF \
                    -o $OUTPUT/featureCounts_output.txt \
                    $bam_files"
        echo $CMD
        eval $CMD
  
      } 2>&1 | tee -a "$LOG"
  
      cat $LOG >> $LOG_03_3
      
    fi
    
  }

  export -f featureCounts_subread
  parallel -j 5 featureCounts_subread ::: $inputs_each
  
  ### Run DEXSeq
  echo "====== `date`: Start running DEXSeq ======" >> $LOG_03_3
  Rscript --no-save --no-restore --verbose ../DEU/03_3_DEXSeq.R \
            --MODE=$MODE \
            --DIR=$DIR \
            --REF=$REF \
            --flat_GTF=$flat_GTF \
            --target=$target \
            --pair=$pair \
            --fdr_cutoff=$fdr_cutoff \
            --seed=$seed \
            --workers=$workers \
            --noOfSim=$noOfSim >> $LOG_03_3 2>&1
  echo "====== `date`: DEXSeq is done! ======" >> $LOG_03_3  

elif [ "$MODE" == "case_study" ]; then
  
  OUTPUT="$DIR/featureCounts_subread/"
  mkdir -p $OUTPUT

  bam_files=""
  for sampleID in `sed '1d' $target| cut -f1`; do
      bam_file="$BAM/$sampleID/${sampleID}.Aligned.sortedByCoord.out.bam"
      bam_files="$bam_files $bam_file"
  done

  ### Run featureCounts
  echo "====== `date`: Start quantifying exon-level reads using subread ======" >> $LOG_03_3
  CMD="featureCounts -p \
              --countReadPairs \
              -M \
              -f \
              -O \
              -T 8 \
              -Q 255 \
              -F 'GTF' \
              -a $flat_GTF \
              -o $OUTPUT/featureCounts_output.txt \
              $bam_files"
  echo $CMD >> $LOG_03_3
  eval $CMD >> $LOG_03_3 2>&1
  
  echo "====== `date`: Truncate header of featureCounts output ======" >> $LOG_03_3
  samples=$(sed '1d' $target| cut -f1| tr '\n' '\t')
  head -1 $OUTPUT/featureCounts_output.txt > $OUTPUT/tmp.txt
  echo -e "Geneid\tChr\tStart\tEnd\tStrand\tLength\t$samples" >> $OUTPUT/tmp.txt
  tail -n +3 $OUTPUT/featureCounts_output.txt >> $OUTPUT/tmp.txt
  mv $OUTPUT/tmp.txt $OUTPUT/featureCounts_output.txt
  
  echo "====== `date`: Start running DEXSeq ======" >> $LOG_03_3
  IFS=',' read -r -a pair <<< "$pair"

  for p in "${pair[@]}"; do
    echo "====== $p ======" >> $LOG_03_3
    IFS='-' read -r g1 g2 <<< "$p"
    awk -v g1="$g1" -v g2="$g2" 'NR==1 || $2 == g1 || $2 == g2' $target > ${target/.tsv/.${p}.tsv}
    samples=$(sed '1d' ${target/.tsv/.${p}.tsv}| cut -f1|tr '\n' '|'| sed 's/|$//')
    idx=$(sed -n '2p' $OUTPUT/featureCounts_output.txt| tr '\t' '\n'| grep -n -E $samples| cut -d: -f1| tr '\n' ','| sed 's/,$//')
    
    mkdir -p $OUTPUT/$p
    head -1 $OUTPUT/featureCounts_output.txt > $OUTPUT/$p/featureCounts_output.txt
    sed '1d' $OUTPUT/featureCounts_output.txt| cut -f1-6,$idx >> $OUTPUT/$p/featureCounts_output.txt
    
    Rscript --no-save --no-restore --verbose ../DEU/03_3_DEXSeq.R \
              --MODE=$MODE \
              --DIR=$DIR \
              --REF=$REF \
              --flat_GTF=$flat_GTF \
              --target=${target/.tsv/.${p}.tsv} \
              --fdr_cutoff=$fdr_cutoff \
              --pair=$p >> $LOG_03_3 2>&1 
    
  done
    
else

  echo "Invalid mode specified. Choose 'simulation' or 'case_study'." > $LOG_03_3  
  
fi

