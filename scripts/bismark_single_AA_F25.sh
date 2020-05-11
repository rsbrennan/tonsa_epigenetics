#!/bin/bash

for sample in `ls ~/tonsa_epigenetics/data/trimmed | grep '.fq.gz$' | cut -f 1-3 -d "_"| uniq | grep 'AA_F25'`

do

    cd ~/tonsa_epigenetics/analysis/aligned/single_end_lane1/


    echo "starting lane 1 sample ${sample} read 1"

     ~/bin/Bismark-0.22.3/bismark --bowtie2 --multicore 1 \
        --genome /data/copepods/tonsa_genome/ \
        --output_dir ~/tonsa_epigenetics/analysis/aligned/single_end_lane1/ \
        ~/tonsa_epigenetics/analysis/aligned/lane1/${sample}_1.fq.gz_unmapped_reads_1.fq.gz \
        --rg_tag --rg_id ${sample} --rg_sample ${sample} --gzip --local

    echo "done with sample ${sample} read 1"

    echo "starting lane 1 sample ${sample} read 2"

     ~/bin/Bismark-0.22.3/bismark --bowtie2 --multicore 1 \
        --genome /data/copepods/tonsa_genome/ \
        --output_dir ~/tonsa_epigenetics/analysis/aligned/single_end_lane1/ \
        ~/tonsa_epigenetics/analysis/aligned/lane1/${sample}_2.fq.gz_unmapped_reads_2.fq.gz \
        --rg_tag --rg_id ${sample} --rg_sample ${sample}  --local --pbat

    echo "done with lane 1  sample ${sample} read 2"

### lane 2

    cd ~/tonsa_epigenetics/analysis/aligned/single_end_lane2/
    echo "starting lane 2 sample ${sample} read 1"

     ~/bin/Bismark-0.22.3/bismark --bowtie2 --multicore 1 \
        --genome /data/copepods/tonsa_genome/ \
        --output_dir ~/tonsa_epigenetics/analysis/aligned/single_end_lane2/ \
        ~/tonsa_epigenetics/analysis/aligned/lane2/${sample}_1.fq.gz_unmapped_reads_1.fq.gz \
        --rg_tag --rg_id ${sample} --rg_sample ${sample} --gzip --local

    echo "done with sample ${sample} read 1"

    echo "starting lane 2 sample ${sample} read 2"

     ~/bin/Bismark-0.22.3/bismark --bowtie2 --multicore 1 \
        --genome /data/copepods/tonsa_genome/ \
        --output_dir ~/tonsa_epigenetics/analysis/aligned/single_end_lane2/ \
        ~/tonsa_epigenetics/analysis/aligned/lane2/${sample}_2.fq.gz_unmapped_reads_2.fq.gz \
        --rg_tag --rg_id ${sample} --rg_sample ${sample} --local --pbat

    echo "done with lane 2 sample ${sample} read 2"



done
