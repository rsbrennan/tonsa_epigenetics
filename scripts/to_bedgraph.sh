cd ~/tonsa_epigenetics/analysis/methylation_extract


for sample in `ls CpG_context_bothlanes*|  cut -f 4-6 -d '_' | cut -f 1 -d '.'`

do

    echo starting ${sample}

    ~/bin/bismark_v0.22.1/bismark2bedGraph --scaffolds -o ${sample}_bothlanes --dir ~/tonsa_epigenetics/analysis/methylation_extract/ CpG_context_bothlanes_${sample}.txt.gz

    echo finished ${sample}

done

