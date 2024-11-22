#!/bin/bash

### Junction filtering based on UMR cutoff of 2
if [[ "$MODE" == "simulation" ]]; then

  mkdir -p $DIR/SJ/$i

  ### Manually filtering junctions supported by less than 3 UMRs  
  LOG=$DIR/SJ/$i/SJ_filtering.log
  echo "====== `date`: Start filtering SJ ======" > $LOG
  CMD1="cp $DIR/aligned_pass1/$i/*/*SJ.out.tab $DIR/SJ/$i"
  echo $CMD1 >> $LOG
  eval $CMD1 >> $LOG
  CMD2="cat $DIR/SJ/$i/*.SJ.out.tab | \
        awk '(\$7 >= 3)'| cut -f1-6| sort| uniq > $DIR/SJ/$i/merged_UMR_3.SJ.tab"
  echo $CMD2 >> $LOG
  eval $CMD2 >> $LOG
  echo "====== `date`: SJ filtering is done! ======" >> $LOG
  
  cat $LOG >> $LOG_02_2

elif [[ "$MODE" == "case_study" ]]; then
  
  ### Manually filtering junctions supported by less than 3 UMRs  
  LOG=$DIR/SJ/SJ_filtering.log
  echo "====== `date`: Start filtering SJ ======"
  CMD1="cp $DIR/aligned_pass1/*/*SJ.out.tab $DIR/SJ"
  echo $CMD1
  eval $CMD1
  CMD2="cat $DIR/SJ/*.SJ.out.tab | \
        awk '(\$7 >= 3 && \$5 > 0)'| cut -f1-6| sort| uniq > $DIR/SJ/merged_UMR_3.SJ.tab"
  echo $CMD2
  eval $CMD2
  echo "====== `date`: SJ filtering is done! ======"

else

  echo "Invalid mode specified. Choose 'simulation' or 'case_study'." >> $LOG_02_2
  
fi

