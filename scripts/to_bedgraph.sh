cd ~/tonsa_epigenetics/analysis/methylation_extract/cpg_merged/


for sample in `ls CpG_context_merged*|  cut -f 4-6 -d '_' | cut -f 1 -d '.'`

do

    echo starting ${sample}

    ~/bin/Bismark-0.22.3/bismark2bedGraph --scaffolds --buffer_size 40% -o ${sample}_merged --dir ~/tonsa_epigenetics/analysis/methylation_extract/cpg_merged/ CpG_context_merged_${sample}.txt.gz

    echo finished ${sample}

done

