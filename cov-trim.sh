#! /bin/bash
# tool for trimming all the samples at once
for i in *_R1_001.fastq.gz
do
   base=$(basename ${i} _R1_001.fastq.gz)
   trimmomatic PE ${i} ${base}_R2_001.fastq.gz \
                ${base}_R1_pair.fastq.gz ${base}_R1_unpair.fastq.gz \
                ${base}_R2_pair.fastq.gz ${base}_R2_unpair.fastq.gz \
                ILLUMINACLIP:/home/<user>/mambaforge/envs/covmap-bsl/share/trimmomatic-0.39-2/adapters/adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25
done
