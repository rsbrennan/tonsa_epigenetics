cd ~/tonsa_epigenetics/analysis/non_dir/

sample=AA_F00_1

    echo "starting sample ${sample}"

    ~/bin/bismark_v0.22.1/bismark --bowtie2 --multicore 2 -n 1 \
        --genome /data/copepods/tonsa_genome/ \
        --output_dir ~/tonsa_epigenetics/analysis/non_dir/ \
        -1 ~/tonsa_epigenetics/data/trimmed/${sample}_1.fq.gz \
        -2 ~/tonsa_epigenetics/data/trimmed/${sample}_2.fq.gz \
        --rg_tag --rg_id ${sample} --rg_sample ${sample} --un --score_min L,0,-0.6 -I 0 -X 1000 -non_directional

    echo "done with sampe ${sample}"
