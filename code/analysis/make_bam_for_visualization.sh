#!/bin/bash 

module load samtools

# DIR="/vast/projects/MM/tam/Differential_splicing/milevskiy_2023_GSE227748/"
DIR="../../data/case_study/GSE227748/"
BAM_DIR="$DIR/aligned_pass2"
OUTPUT="$DIR/Gviz/"
REF="../../annotation"
GTF="$REF/gencode.vM32.annotation.gtf"
geneAnno="$REF/gencode.vM32.annotation.genes.bed"

mkdir -p $OUTPUT

### 1. Extract genes
if [[ ! -e ${GTF/.gtf/.genes.bed} ]]; then
  (echo -e "chr\tstart\tend\tgeneID\tscore\tstrand\tgeneSymbol"; \
  grep -v "^#" $GTF| \
  awk 'BEGIN {FS=OFS="\t"} {split($9, a, "; ")} {if ($3=="gene") print $1, $4-1, $5, a[1], "0", $7, a[3]}'| \
  sed -e 's/gene_id "//g' -e 's/gene_name "//g' -e 's/"//g') > ${GTF/.gtf/.genes.bed}
fi

### Index bam files
for SAMPLE in `ls $BAM_DIR`; do
    echo $SAMPLE
    BAM="$DIR/aligned_pass2/${SAMPLE}/${SAMPLE}.Aligned.sortedByCoord.out.bam"
    echo $BAM
    if [[ ! -e ${BAM}.bai ]]; then
      samtools index -@ 16 $BAM
    fi
done

##############################################################################
gene="Numb"
d="1000"

gene_region=$(awk -v g=$gene -v d=$d 'BEGIN {FS=OFS="\t"} {if($7==g) print $1, $2-d, $3+d, $4}' $geneAnno)
chr=$(echo "$gene_region"| cut -f1)
start=$(echo "$gene_region"| awk -v d=$d '{print $2-d}')
end=$(echo "$gene_region"| awk -v d=$d '{print $3+d}')
position=${chr}:${start}-${end}
geneID=$(echo "$gene_region"| cut -f4)

mkdir -p $OUTPUT/${gene}_${geneID}_${d}

for SAMPLE in `ls $BAM_DIR`; do
    echo $SAMPLE
    BAM="$DIR/aligned_pass2/${SAMPLE}/${SAMPLE}.Aligned.sortedByCoord.out.bam"
    echo $BAM
    samtools view -@ 16 -b -h $BAM $position > $OUTPUT/${gene}_${geneID}_${d}/${SAMPLE}.${gene}.bam
    samtools index -@ 16 $OUTPUT/${gene}_${geneID}_${d}/${SAMPLE}.${gene}.bam
done

##############################################################################
gene="Fgfr1"
d="1000"

gene_region=$(awk -v g=$gene -v d=$d 'BEGIN {FS=OFS="\t"} {if($7==g) print $1, $2-d, $3+d, $4}' $geneAnno)
chr=$(echo "$gene_region"| cut -f1)
start=$(echo "$gene_region"| awk -v d=$d '{print $2-d}')
end=$(echo "$gene_region"| awk -v d=$d '{print $3+d}')
position=${chr}:${start}-${end}
geneID=$(echo "$gene_region"| cut -f4)

mkdir -p $OUTPUT/${gene}_${geneID}_${d}

for SAMPLE in `ls $BAM_DIR`; do
    echo $SAMPLE
    BAM="$DIR/aligned_pass2/${SAMPLE}/${SAMPLE}.Aligned.sortedByCoord.out.bam"
    echo $BAM
    samtools view -@ 16 -b -h $BAM $position > $OUTPUT/${gene}_${geneID}_${d}/${SAMPLE}.${gene}.bam
    samtools index -@ 16 $OUTPUT/${gene}_${geneID}_${d}/${SAMPLE}.${gene}.bam
done

##############################################################################
gene="Dusp16"
d="1000"

gene_region=$(awk -v g=$gene -v d=$d 'BEGIN {FS=OFS="\t"} {if($7==g) print $1, $2-d, $3+d, $4}' $geneAnno)
chr=$(echo "$gene_region"| cut -f1)
start=$(echo "$gene_region"| awk -v d=$d '{print $2-d}')
end=$(echo "$gene_region"| awk -v d=$d '{print $3+d}')
position=${chr}:${start}-${end}
geneID=$(echo "$gene_region"| cut -f4)

mkdir -p $OUTPUT/${gene}_${geneID}_${d}

for SAMPLE in `ls $BAM_DIR`; do
    echo $SAMPLE
    BAM="$DIR/aligned_pass2/${SAMPLE}/${SAMPLE}.Aligned.sortedByCoord.out.bam"
    echo $BAM
    samtools view -@ 16 -b -h $BAM $position > $OUTPUT/${gene}_${geneID}_${d}/${SAMPLE}.${gene}.bam
    samtools index -@ 16 $OUTPUT/${gene}_${geneID}_${d}/${SAMPLE}.${gene}.bam
done
