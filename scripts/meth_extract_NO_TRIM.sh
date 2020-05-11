#!/bin/bash


~/bin/bismark_v0.22.1/bismark_methylation_extractor --bedGraph --scaffolds --gzip \
    --cytosine_report --comprehensive \
    --no_header \
    --parallel 2 \
    --output ~/tonsa_epigenetics/analysis/methylation_extract/bias_check \
    --genome_folder /data/copepods/tonsa_genome/ \
    ~/tonsa_epigenetics/analysis/aligned/lane1/AA_F25_1_1_bismark_bt2_pe.bam

~/bin/bismark_v0.22.1/bismark_methylation_extractor --bedGraph --scaffolds --gzip \
    --cytosine_report --comprehensive \
    --no_header \
    --parallel 2 \
    --output ~/tonsa_epigenetics/analysis/methylation_extract/bias_check \
    --genome_folder /data/copepods/tonsa_genome/ \
    ~/tonsa_epigenetics/analysis/aligned/lane1/AH_F25_1_1_bismark_bt2_pe.bam

~/bin/bismark_v0.22.1/bismark_methylation_extractor --bedGraph --scaffolds --gzip \
    --cytosine_report --comprehensive \
    --no_header \
    --parallel 2 \
    --output ~/tonsa_epigenetics/analysis/methylation_extract/bias_check \
    --genome_folder /data/copepods/tonsa_genome/ \
    ~/tonsa_epigenetics/analysis/aligned/lane1/HA_F25_1_1_bismark_bt2_pe.bam

~/bin/bismark_v0.22.1/bismark_methylation_extractor --bedGraph --scaffolds --gzip \
    --cytosine_report --comprehensive \
    --no_header \
    --parallel 2 \
    --output ~/tonsa_epigenetics/analysis/methylation_extract/bias_check \
    --genome_folder /data/copepods/tonsa_genome/ \
    ~/tonsa_epigenetics/analysis/aligned/lane1/HH_F25_1_1_bismark_bt2_pe.bam
