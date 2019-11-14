~/bin/bismark_v0.22.1/bismark_methylation_extractor --bedGraph --scaffolds --gzip \
    --cytosine_report --comprehensive \
    --ignore 2 --ignore_r2 2 \
    --no_header \
    --parallel 4 \
    --output ~/tonsa_epigenetics/analysis/methylation_extract/lane2/ \
    --genome_folder /data/copepods/tonsa_genome/ \
    ~/tonsa_epigenetics/analysis/aligned/lane2/HA*pe.bam

# unmapped
~/bin/bismark_v0.22.1/bismark_methylation_extractor --bedGraph --scaffolds --gzip \
    --cytosine_report --comprehensive \
    --ignore 2 \
    --no_header \
    --parallel 4 \
    --output ~/tonsa_epigenetics/analysis/methylation_extract/lane2/ \
    --genome_folder /data/copepods/tonsa_genome/ \
    ~/tonsa_epigenetics/analysis/aligned/lane2/HA*bt2.bam
