#!/bin/bash 

# projName=""
bioPrjID="PRJNA946753"
DIR="/vast/projects/MM/tam/Differential_splicing/milevskiy_2023_GSE227748/"
isPairedEnd=TRUE

mkdir -p $DIR/log $DIR/target

# Download sample info from ENA browser
dl_list="$DIR/target/sampleInfo_${bioPrjID}.txt"
wget 'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA946753&result=read_run&fields=sample_accession,run_accession,library_layout,library_strategy,experiment_title,fastq_md5,fastq_ftp,fastq_aspera,sample_title&format=tsv&download=true&limit=0' \
  -O $dl_list
  
sampleID_idx=$(head -1 $dl_list| tr '\t' '\n'| grep -n 'sample_title'| cut -d: -f1)
(head -1 $dl_list; sed '1d' $dl_list| sed 's/ cells, biol rep /_rep/g' | sort -k$sampleID_idx,$sampleID_idx) > ${dl_list/.txt/.1.txt}
mv ${dl_list/.txt/.1.txt} $dl_list

# Make target.tsv (contains sampleID and group)
(echo -e "sampleID\tgroup"; sed '1d' $dl_list| cut -f$sampleID_idx| awk -F"_" '{print $1"_"$2"\t"$1}') > $DIR/target/target.tsv

# Make decoder.bySample.txt (used exclusively for JunctionSeq)
(echo -e "sample.ID\tgroup.ID"; sed '1d' $DIR/target/target.tsv) > $DIR/target/decoder.bySample.txt

export DIR="$DIR"
export RAW="$RAW"
export bioPrjID="$bioPrjID"
export isPairedEnd="$isPairedEnd"

### Part 0: Download fastq files using fastq_dump and SRR... in $dl_list
job_00_1_id=$(sbatch --parsable run_00_1.sh)

### Part 1: QC and read trimming
job_01_1_id=$(sbatch --parsable --dependency=afterok:$job_00_1_id run_01_1.sh)
job_01_2_id=$(sbatch --parsable --dependency=afterok:$job_01_1_id run_01_2.sh)

# Check maxRL from FASTQC report
maxRL=81
use_default=FALSE

if [[ $use_default == FALSE ]]; then
  maxRL="$maxRL"
else 
  maxRL=100
fi

export maxRL="$maxRL"

### Part 2: Read alignment using STAR
# Check if indexed genome exists
genomeDir="../../data/genomeIndex_${maxRL}bp"

if [[ ! -d "$genomeDir" ]]; then

  # Genome index 
  job_02_0_id=$(sbatch --parsable --dependency=afterok:$job_01_2_id run_02_0.sh)
  
  # 1-pass mapping
  job_02_1_id=$(sbatch --parsable --dependency=afterok:$job_02_0_id run_02_1.sh)
  
else

  # 1-pass mapping
  job_02_1_id=$(sbatch --parsable --dependency=afterany:$job_01_2_id run_02_1.sh)

fi

### Genome reindexing with filtered SJ
job_02_2_3_id=$(sbatch --dependency=afterany:$job_02_1_id run_02_2_3.sh)

### 2-pass mapping
job_02_4_id=$(sbatch --dependency=afterany:$job_02_2_3_id run_02_4.sh)

### Part 3: DEU analysis
# pair="Basal-LP,Basal-ML,LP-ML"
pair="Basal-ML"

export pair="$pair"

# edgeR::diffSpliceDGE
job_03_1_id=$(sbatch --parsable --dependency=afterany:$job_02_2_4_id run_03_1.sh)

# limma::diffSplice
job_03_2_id=$(sbatch --parsable --dependency=afterany:$job_02_2_4_id:$job_03_1_id run_03_2.sh)

# DEXSeq
job_03_3_id=$(sbatch --parsable --dependency=afterany:$job_02_2_4_id run_03_3.sh)

# JunctionSeq
job_03_4_id=$(sbatch --parsable --dependency=afterany:$job_02_2_4_id run_03_4.sh)

