#!/bin/bash

# merge methylation calls from extracted files

cd ~/tonsa_epigenetics/analysis/methylation_extract/cpg_merged/

for sample in `ls ~/tonsa_epigenetics/analysis/methylation_extract/lane1 | grep 'CpG_context' | cut -f 3-5 -d '_' | cut -f 1 -d '.' | sort | uniq`

do

echo starting ${sample}

zcat ~/tonsa_epigenetics/analysis/methylation_extract/lane1/CpG_context_${sample}_1_bismark_bt2_pe.txt.gz \
    ~/tonsa_epigenetics/analysis/methylation_extract/lane2/CpG_context_${sample}_1_bismark_bt2_pe.txt.gz \
    ~/tonsa_epigenetics/analysis/methylation_extract/single_end_lane1/CpG_context_${sample}_1.fq.gz_unmapped_reads_1_bismark_bt2.txt.gz\
    ~/tonsa_epigenetics/analysis/methylation_extract/single_end_lane1/CpG_context_${sample}_2.fq.gz_unmapped_reads_2_bismark_bt2.txt.gz\
    ~/tonsa_epigenetics/analysis/methylation_extract/single_end_lane2/CpG_context_${sample}_1.fq.gz_unmapped_reads_1_bismark_bt2.txt.gz\
    ~/tonsa_epigenetics/analysis/methylation_extract/single_end_lane2/CpG_context_${sample}_2.fq.gz_unmapped_reads_2_bismark_bt2.txt.gz |\
     gzip \
     > ~/tonsa_epigenetics/analysis/methylation_extract/cpg_merged/CpG_context_merged_${sample}.txt.gz


echo finished ${sample}

done

## check that it worked correctly
 ## total number should add up to sum of indiv files

for sample in `ls ~/tonsa_epigenetics/analysis/methylation_extract/lane1 | grep 'CpG_context' | cut -f 3-5 -d '_' | cut -f 1 -d '.' | sort | uniq`

do

numtot=$(zcat ~/tonsa_epigenetics/analysis/methylation_extract/cpg_merged/CpG_context_merged_${sample}.txt.gz | wc -l)
num1=$(zcat ~/tonsa_epigenetics/analysis/methylation_extract/lane1/CpG_context_${sample}_1_bismark_bt2_pe.txt.gz | wc -l)
num2=$(zcat ~/tonsa_epigenetics/analysis/methylation_extract/lane2/CpG_context_${sample}_1_bismark_bt2_pe.txt.gz| wc -l)
num3=$(zcat ~/tonsa_epigenetics/analysis/methylation_extract/single_end_lane1/CpG_context_${sample}_1.fq.gz_unmapped_reads_1_bismark_bt2.txt.gz| wc -l)
num4=$(zcat ~/tonsa_epigenetics/analysis/methylation_extract/single_end_lane1/CpG_context_${sample}_2.fq.gz_unmapped_reads_2_bismark_bt2.txt.gz| wc -l)
num5=$(zcat ~/tonsa_epigenetics/analysis/methylation_extract/single_end_lane2/CpG_context_${sample}_1.fq.gz_unmapped_reads_1_bismark_bt2.txt.gz| wc -l)
num6=$(zcat ~/tonsa_epigenetics/analysis/methylation_extract/single_end_lane2/CpG_context_${sample}_2.fq.gz_unmapped_reads_2_bismark_bt2.txt.gz| wc -l)

addtot=$((num1 + num2 + num3 + num4 + num5 + num6))

if [ "$addtot" = "$numtot" ]; then
    echo "$sample Values are equal"
else
    echo "$sample Values are not equal"
fi

done

# and use this merged output as input for bismark2bedGraph.
# The only thing that is somewhat inconvenient is that you don't get one nice mapping report / html file because the alignment process was split up

cd ~/tonsa_epigenetics/analysis/methylation_extract/cpg_merged/


for sample in `ls CpG_context_merged*|  cut -f 4-6 -d '_' | cut -f 1 -d '.'`

do

    echo starting ${sample}

    ~/bin/Bismark-0.22.3/bismark2bedGraph --scaffolds --buffer_size 40% -o ${sample}_merged --dir ~/tonsa_epigenetics/analysis/methylation_extract/cpg_merged/ CpG_context_merged_${sample}.txt.gz

    echo finished ${sample}

done

