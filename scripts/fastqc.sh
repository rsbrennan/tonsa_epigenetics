output=~/tonsa_epigenetics/analysis/fastqc/pre_trim/

cd /data/copepods/tonsa_epigenetics

~/bin/FastQC/fastqc *fq.gz -t 4 -f fastq --noextract -o ${output}

# aggregate fastqc across files:

multiqc ~/tonsa_epigenetics/analysis/fastqc/pre_trim/
