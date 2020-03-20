cd ~/tonsa_epigenetics/analysis/aligned/lane2/

for sample in `ls ~/tonsa_epigenetics/data/trimmed/lane2/ | grep '.fq.gz$' | cut -f 1-3 -d "_"| uniq | grep 'HH_F03'`

do

    echo "starting sample ${sample}"

    ~/bin/bismark_v0.22.1/bismark --bowtie2 --multicore 1 \
        --genome /data/copepods/tonsa_genome/ \
        --output_dir ~/tonsa_epigenetics/analysis/aligned/lane1/ \
        -1 ~/tonsa_epigenetics/data/trimmed/${sample}_1.fq.gz \
        -2 ~/tonsa_epigenetics/data/trimmed/${sample}_2.fq.gz \
        --rg_tag --rg_id ${sample} --rg_sample ${sample} --un --gzip --local

    echo "done with sample ${sample}"

done

