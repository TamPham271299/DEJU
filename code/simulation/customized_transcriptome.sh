#!/bin/bash 

# Function: From genomic intervals to transcriptome
# Input: Genomic intervals of transcripts designed for alternative splicing pattern (BED12)
# Output: Customized transcriptome (FASTA)

module load bedtools/2.31.1

cd ../../data/simulation/customized_transcriptome

FASTA="../../../annotation/GRCm39.primary_assembly.genome.fa"

### ES
# 1st transcript
bedtools groupby -i <(sed '1d' ES_t1.pre.bed12.tsv| tr '+' 'a') -g 1 -c 2,3,4,5,1,6,7 -o first,min,max,first,count,collapse,collapse| tr 'a' '+'| \
  awk -F'\t' 'BEGIN{OFS="\t";}{print $2, $3, $4, $1, "0", $5, $4, $4, "0", $6, $7, $8}'| \
  sort -k1,1 -k2,2n > ES_t1.bed12.tsv
# 2nd transcript
bedtools groupby -i <(sed '1d' ES_t2.pre.bed12.tsv| tr '+' 'a') -g 1 -c 2,3,4,5,1,6,7 -o first,min,max,first,count,collapse,collapse| tr 'a' '+'| \
  awk -F'\t' 'BEGIN{OFS="\t";}{print $2, $3, $4, $1, "0", $5, $4, $4, "0", $6, $7, $8}'| \
  sort -k1,1 -k2,2n > ES_t2_temp.bed12.tsv
paste <(sort -k4,4 ES_t2_temp.bed12.tsv) <(sed '1d' ES_exon.info.tsv| cut -f1,9| sort -k1,1)| \
  awk '{OFS="\t"}{print $1, $2, $3, $4".skip_exon"$14, $5, $6, $7, $8, $9, $10, $11, $12}'| \
  sort -k1,1 -k2,2n > ES_t2.bed12.tsv

### MXE
# 1st transcript
bedtools groupby -i <(sed '1d' MXE_t1.pre.bed12.tsv) -g 1 -c 2,3,4,5,1,6,7 -o first,min,max,first,count,collapse,collapse| \
  awk -F'\t' 'BEGIN{OFS="\t";}{print $2, $3, $4, $1, "0", $5, $4, $4, "0", $6, $7, $8}'| \
  sort -k1,1 -k2,2n > MXE_t1_temp.bed12.tsv
paste <(sort -k4,4 MXE_t1_temp.bed12.tsv) <(sed '1d' MXE_exon_t1.info.tsv| cut -f1,9| sort -k1,1)| \
  awk '{OFS="\t"}{print $1, $2, $3, $4".skip_exon"$14, $5, $6, $7, $8, $9, $10, $11, $12}'| \
  sort -k1,1 -k2,2n > MXE_t1.bed12.tsv
# 2nd transcript
bedtools groupby -i <(sed '1d' MXE_t2.pre.bed12.tsv) -g 1 -c 2,3,4,5,1,6,7 -o first,min,max,first,count,collapse,collapse| \
  awk -F'\t' 'BEGIN{OFS="\t";}{print $2, $3, $4, $1, "0", $5, $4, $4, "0", $6, $7, $8}'| \
  sort -k1,1 -k2,2n > MXE_t2_temp.bed12.tsv
paste <(sort -k4,4 MXE_t2_temp.bed12.tsv) <(sed '1d' MXE_exon_t2.info.tsv| cut -f1,9| sort -k1,1)| \
  awk '{OFS="\t"}{print $1, $2, $3, $4".skip_exon"$14, $5, $6, $7, $8, $9, $10, $11, $12}'| \
  sort -k1,1 -k2,2n > MXE_t2.bed12.tsv

### ASS
# 1st transcript
bedtools groupby -i <(sed '1d' ASS_t1.pre.bed12.tsv) -g 1 -c 2,3,4,5,1,6,7 -o first,min,max,first,count,collapse,collapse| \
  awk -F'\t' 'BEGIN{OFS="\t";}{print $2, $3, $4, $1, "0", $5, $4, $4, "0", $6, $7, $8}'| \
  sort -k1,1 -k2,2n > ASS_t1.bed12.tsv
# 2nd transcript
bedtools groupby -i <(sed '1d' ASS_t2.pre.bed12.tsv) -g 1 -c 2,3,4,5,1,6,7 -o first,min,max,first,count,collapse,collapse| \
  awk -F'\t' 'BEGIN{OFS="\t";}{print $2, $3, $4, $1, "0", $5, $4, $4, "0", $6, $7, $8}'| \
  sort -k1,1 -k2,2n > ASS_t2_temp.bed12.tsv
paste <(sort -k4,4 ASS_t2_temp.bed12.tsv) <(sed '1d' ASS_exon.info.tsv| cut -f1,9| sort -k1,1)| \
  awk '{OFS="\t"}{print $1, $2, $3, $4".splice_site_at_exon"$14, $5, $6, $7, $8, $9, $10, $11, $12}'| \
  sort -k1,1 -k2,2n > ASS_t2.bed12.tsv

### RI
# 1st transcript
bedtools groupby -i <(sed '1d' RI_t1.pre.bed12.tsv| tr '+' 'a') -g 1 -c 2,3,4,5,1,6,7 -o first,min,max,first,count,collapse,collapse| tr 'a' '+'| \
  awk -F'\t' 'BEGIN{OFS="\t";}{print $2, $3, $4, $1, "0", $5, $4, $4, "0", $6, $7, $8}'| \
  sort -k1,1 -k2,2n > RI_t1.bed12.tsv
# 2nd transcript
bedtools groupby -i <(sed '1d' RI_t2.pre.bed12.tsv| tr '+' 'a') -g 1 -c 2,3,4,5,1,6,7 -o first,min,max,first,count,collapse,collapse| tr 'a' '+'| \
  awk -F'\t' 'BEGIN{OFS="\t";}{print $2, $3, $4, $1, "0", $5, $4, $4, "0", $6, $7, $8}'| \
  sort -k1,1 -k2,2n > RI_t2_temp.bed12.tsv
paste <(sort -k4,4 RI_t2_temp.bed12.tsv) <(sed '1d' RI_exon.info.tsv| cut -f1,9| sort -k1,1)| \
  awk '{OFS="\t"}{print $1, $2, $3, $4".retained_intron_"$14, $5, $6, $7, $8, $9, $10, $11, $12}'| \
  sort -k1,1 -k2,2n > RI_t2.bed12.tsv

### combine ES, MXE, ASS and RI into one BED12 file
cat ES_t1.bed12.tsv > AS.bed12.tsv
cat ES_t2.bed12.tsv >> AS.bed12.tsv
cat MXE_t1.bed12.tsv >> AS.bed12.tsv
cat MXE_t2.bed12.tsv >> AS.bed12.tsv
cat ASS_t1.bed12.tsv >> AS.bed12.tsv
cat ASS_t2.bed12.tsv >> AS.bed12.tsv
cat RI_t1.bed12.tsv >> AS.bed12.tsv
cat RI_t2.bed12.tsv >> AS.bed12.tsv
sort -k4,4 -o AS_sorted.bed12.tsv AS.bed12.tsv

### Extract sequences from reference FASTA based on genomic intervals in BED12
bedtools getfasta -fi $FASTA -bed AS_sorted.bed12.tsv -split -nameOnly -fo gencode.vM32.custom.transcriptome.tmp.fa
fold -w 60 gencode.vM32.custom.transcriptome.tmp.fa > gencode.vM32.custom.transcriptome.fa
