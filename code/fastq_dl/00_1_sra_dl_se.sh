#!/bin/bash
module load sra-toolkit/3.1.0
module load parallel/20240722

### Extract sample information
col_samplename=`head -1 $dl_list| tr "\t" "\n"| nl| grep "sample_title"| awk '{print $1}'`
col_lane=`head -1 $dl_list| tr "\t" "\n"| nl| grep "run_accession"| awk '{print $1}'`

inputs=`awk -F'\t' -v col_samplename=$col_samplename -v col_lane=$col_lane 'NR>1{print $col_samplename";"$col_lane}' $dl_list| paste -s -d"\t"`

inputs_each=""
for i in $inputs; do
    sample_name=`echo $i| cut -d";" -f1`
    lane=`echo $i| cut -d";" -f2`
    inputs_each=`echo -e $inputs_each$sample_name";"$lane"\t"`
done

### Download FASTQ files using fastq_dump
download () {
    sample_name=`echo $1| cut -d";" -f1`
    lane=`echo $1| cut -d";" -f2`
    echo $sample_name
    echo $lane
    echo $dl_dir
    mkdir -p $dl_dir/$sample_name/$lane
    cd $dl_dir/$sample_name/$lane
    fastq-dump --gzip $lane
}

export -f download
parallel -j 8 download ::: $inputs_each
