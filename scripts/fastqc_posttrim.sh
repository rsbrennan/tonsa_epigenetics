#!/bin/bash

output=~/tonsa_epigenetics/analysis/fastqc/post_trim/lane1

cd ~/tonsa_epigenetics/data/trimmed/

~/bin/FastQC/fastqc *fq.gz -t 1 -f fastq --noextract -o ${output}

# aggregate fastqc across files:

#multiqc ~/tonsa_epigenetics/analysis/fastqc/post_trim/

output=~/tonsa_epigenetics/analysis/fastqc/post_trim/lane2
cd ~/tonsa_epigenetics/data/trimmed/lane2

~/bin/FastQC/fastqc *fq.gz -t 1 -f fastq --noextract -o ${output}

