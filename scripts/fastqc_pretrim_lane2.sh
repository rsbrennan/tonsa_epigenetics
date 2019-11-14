output=~/tonsa_epigenetics/analysis/fastqc/lane2/pre_trim/

cd /data/copepods/tonsa_epigenetics/lane2

~/bin/FastQC/fastqc *fq.gz -t 4 -f fastq --noextract -o ${output}

# aggregate fastqc across files:

multiqc ~/tonsa_epigenetics/analysis/fastqc/lane2/pre_trim/
