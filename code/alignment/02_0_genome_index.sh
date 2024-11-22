#!/bin/bash

### Indexing reference genome using STAR with optimal sjdbOverhang
echo "====== `date`: Indexing genome using STAR with optimal sjdbOverhang ======" > $LOG_02_0
mkdir -p $genomeDir
CMD="STAR --runThreadN 16 \
          --runMode genomeGenerate \
          --genomeDir $genomeDir \
          --genomeFastaFiles $FASTA \
          --sjdbGTFfile $GTF \           
          --sjdbOverhang $sjdbOverhang"
echo $CMD >> $LOG_02_0
eval $CMD >> $LOG_02_0
echo "====== `date`: Indexing genome is done! ======" >> $LOG_02_0