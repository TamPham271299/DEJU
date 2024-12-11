#!/bin/bash

module load R/4.3.3

cd ../../annotation

### Download GTF and FASTA mouse reference genome GRCm39 vM32
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M32/gencode.vM32.annotation.gtf.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M32/gencode.vM32.primary_assembly.annotation.gtf.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M32/GRCm39.primary_assembly.genome.fa.gz
gunzip -k gencode.vM32.annotation.gtf.gz
gunzip -k gencode.vM32.primary_assembly.annotation.gtf.gz
gunzip -k GRCm39.primary_assembly.genome.fa.gz 

### Download gene_info_GRCm39_M32_ensembl109.tsv

### Extracting protein-coding genes to create a customized transcriptome
grep "protein_coding" gencode.vM32.annotation.gtf| gzip > gencode.vM32.PC.annotation.gtf.gz

### Create flattened exon annotation (SAF) from gencode.vM32.annotation.gtf
cd ../code/annotation_dl
Rscript --no-save --no-restore --verbose GTF2SAF.R > GTF2SAF.log 2>&1

### Create SJ database (gencode.vM32.primary_assembly.annotation.SJdatabase.tsv)
GTF="../../annotation/gencode.vM32.primary_assembly.annotation.gtf"
Rscript --no-save --no-restore --verbose GTF2SJdatabase.R \
                                --GTF=$GTF > GTF2SJdatabase.log 2>&1
                                
### Create annotation for DEXseq (gencode.vM32.annotation.flat.DEXSeq.gtf)
# Refer to https://bioconductor.org/packages/release/bioc/vignettes/DEXSeq/inst/doc/DEXSeq.html
# Refer to https://github.com/vivekbhr/Subread_to_DEXSeq
git clone https://github.com/vivekbhr/Subread_to_DEXSeq.git ../../tools/

module load python/3.11.8
# cd /vast/projects/lab_chen/tam/tools
# . DAS/bin/activate
# python --version
Subread_to_DEXSeq="../../tools/Subread_to_DEXSeq"
GTF="../../annotation/gencode.vM32.annotation.gtf"
python $Subread_to_DEXSeq/dexseq_prepare_annotation2.py -f ${GTF/.gtf/.flat.DEXSeq.gtf} $GTF ${GTF/.gtf/.flat.DEXSeq.gff}

### Modify load_SubreadOutput.R to truncate exonID to be suitable for featureCounts
# Modified: From line 43 - 50
# Save as load_SubreadOutput_modified.R




