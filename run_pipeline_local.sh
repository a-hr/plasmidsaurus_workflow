#!/usr/bin/env bash

# This script runs the pipeline locally on your machine.

# ---- Input Arguments ----

# alignment options
fastq_dir=fastqs/
fa_path=ref_genome/reference.fa
gtf_path=ref_genome/reference.gtf
bed_path=ref_genome/reference.bed

# sashimi plot options
sashimi_config=inputs/config.tab
palette_path=inputs/palette.txt
bam_tsv=inputs/input_bams.tsv

sashimi_min_cov=3
sashimi_agg=mean_j
sashimi_alpha=0.6

# output options
output_dir=results/

# ---- Step 1: Quality control ----
mkdir -p tmp/
mkdir -p $output_dir/reports/

# Run FastQC on the raw data
fastqc $fastq_dir/*.fastq.gz -o tmp/

# Collect the FastQC reports into a single HTML file
multiqc tmp/ -o $output_dir/reports -n plasmidsaurus_multiqc_report -i "PlasmidSaurus MultiQC Report"

# ---- Step 2: Alignment ----
mkdir -p $output_dir/bams/

for fastq in $fastq_dir/*.fastq.gz; 
do
    outName=$output_dir/bams/$(basename $fastq .fastq.gz)

    minimap2 -t 8 -ax splice --junc-bed ${bed_path} -c --MD ${fa_path} ${fastq} \\
        | samtools view -b - \\
        | samtools sort -o ${outName}.bam

    samtools index ${outName}.bam
done


# ---- Step 3: Plotting ----
mkdir -p $output_dir/plots/

# loop over all the columns of the config file and extract first two columns
while read -r eventID coords
do  
    echo "Plotting $eventID..."
    docker run -w $PWD -v $PWD:$PWD guigolab/ggsashimi \
        -b $bam_tsv \
        -c $coords \
        -g $gtf_path \
        -P $palette \
        -o ${output_dir}/plots/sashimi_$eventID \
        -M $sashimi_min_cov \
        -C 3 \
        -O 3 \
        -A $sashimi_agg \
        --alpha $sashimi_alpha \
        --fix-y-scale \
        --ann-height 5 \
        --width 15
done < $sashimi_config