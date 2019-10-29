cd ~/tonsa_epigenetics/analysis/lambda/aligned



##############
### for script

~/bin/bismark_v0.22.1/bismark_methylation_extractor --gzip --bedGraph --scaffolds \
  --cytosine_report --comprehensive \
  --ignore 2 --ignore_r2 2 \
  --no_header \
  --parallel 4 \
  --output ~/tonsa_epigenetics/analysis/lambda/methylation_extract/ \
  --genome_folder ~/tonsa_epigenetics/analysis/lambda/Bisulfite_Genome/ \
  ~/tonsa_epigenetics/analysis/lambda/aligned/*pe.bam

