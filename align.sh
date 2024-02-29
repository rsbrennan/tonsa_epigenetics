#!/bin/bash

cd ~/tonsa_epigenetics/analysis/aligned_hiRise/lane1/

eval "$(conda shell.bash hook)"
conda activate bismark

for sample in `ls ~/tonsa_epigenetics/data/trimmed | grep '.fq.gz$' | cut -f 1-3 -d "_"| uniq`

do

    echo "starting sample ${sample}"

    ~/bin/Bismark-0.24.0/bismark --bowtie2 --multicore 1 \
        --genome /data/copepods/tonsa_genome/HiRise_denovo \
        --output_dir ~/tonsa_epigenetics/analysis/aligned_hiRise/lane1/ \
        -1 ~/tonsa_epigenetics/data/trimmed/${sample}_1.fq.gz \
        -2 ~/tonsa_epigenetics/data/trimmed/${sample}_2.fq.gz \
        --rg_tag --rg_id ${sample} --rg_sample ${sample} --un --gzip --local --maxins 1000

    echo "done with sample ${sample}"

done


# lane 2

cd ~/tonsa_epigenetics/analysis/aligned_hiRise/lane2/

eval "$(conda shell.bash hook)"
conda activate bismark

for sample in `ls ~/tonsa_epigenetics/data/trimmed/lane2/ | grep '.fq.gz$' | cut -f 1-3 -d "_"| uniq`

do

    echo "starting sample ${sample}"

    ~/bin/Bismark-0.24.0/bismark --bowtie2 --multicore 1 \
        --genome /data/copepods/tonsa_genome/HiRise_denovo \
        --output_dir ~/tonsa_epigenetics/analysis/aligned_hiRise/lane2/ \
        -1 ~/tonsa_epigenetics/data/trimmed/lane2/${sample}_1.fq.gz \
        -2 ~/tonsa_epigenetics/data/trimmed/lane2/${sample}_2.fq.gz \
        --rg_tag --rg_id ${sample} --rg_sample ${sample} --un --gzip --local --maxins 1000

    echo "done with sample ${sample}"

done

# align single end reads to try to improve mapping

for sample in `ls ~/tonsa_epigenetics/data/trimmed | grep '.fq.gz$' | cut -f 1-3 -d "_"| uniq `

do

    cd ~/tonsa_epigenetics/analysis/aligned_hiRise/single_end_lane1/

    echo "starting lane 1 sample ${sample} read 1"

     ~/bin/Bismark-0.24.0/bismark --bowtie2 --multicore 1 \
        --genome /data/copepods/tonsa_genome/HiRise_denovo/ \
        --output_dir ~/tonsa_epigenetics/analysis/aligned_hiRise/single_end_lane1/ \
        ~/tonsa_epigenetics/analysis/aligned_hiRise/lane1/${sample}_1.fq.gz_unmapped_reads_1.fq.gz \
        --rg_tag --rg_id ${sample} --rg_sample ${sample} --gzip --local

    echo "done with sample ${sample} read 1"

    echo "starting lane 1 sample ${sample} read 2"

     ~/bin/Bismark-0.24.0/bismark --bowtie2 --multicore 1 \
                --genome /data/copepods/tonsa_genome/HiRise_denovo/ \
                --output_dir ~/tonsa_epigenetics/analysis/aligned_hiRise/single_end_lane1/ \
        ~/tonsa_epigenetics/analysis/aligned_hiRise/lane1/${sample}_2.fq.gz_unmapped_reads_2.fq.gz \
        --rg_tag --rg_id ${sample} --rg_sample ${sample}  --local --pbat

    echo "done with lane 1  sample ${sample} read 2"

### lane 2

    cd ~/tonsa_epigenetics/analysis/aligned_hiRise/single_end_lane2/
    echo "starting lane 2 sample ${sample} read 1"

     ~/bin/Bismark-0.24.0/bismark --bowtie2 --multicore 1 \
        --genome /data/copepods/tonsa_genome/HiRise_denovo/ \
        --output_dir ~/tonsa_epigenetics/analysis/aligned_hiRise/single_end_lane2/ \
        ~/tonsa_epigenetics/analysis/aligned_hiRise/lane2/${sample}_1.fq.gz_unmapped_reads_1.fq.gz \
        --rg_tag --rg_id ${sample} --rg_sample ${sample} --gzip --local

    echo "done with sample ${sample} read 1"

    echo "starting lane 2 sample ${sample} read 2"

     ~/bin/Bismark-0.24.0/bismark --bowtie2 --multicore 1 \
        --genome /data/copepods/tonsa_genome/HiRise_denovo/ \
        --output_dir ~/tonsa_epigenetics/analysis/aligned_hiRise/single_end_lane2/ \
        ~/tonsa_epigenetics/analysis/aligned_hiRise/lane2/${sample}_2.fq.gz_unmapped_reads_2.fq.gz \
        --rg_tag --rg_id ${sample} --rg_sample ${sample} --local --pbat

    echo "done with lane 2 sample ${sample} read 2"

done
