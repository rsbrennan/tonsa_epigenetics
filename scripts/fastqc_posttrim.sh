output=~/tonsa_epigenetics/analysis/fastqc/post_trim/

cd ~/tonsa_epigenetics/data/trimmed/

~/bin/FastQC/fastqc *fq.gz -t 4 -f fastq --noextract -o ${output}

# aggregate fastqc across files:

multiqc ~/tonsa_epigenetics/analysis/fastqc/post_trim/
