cd ~/tonsa_epigenetics/analysis/aligned/local_align/lane2

for sample in `ls ~/tonsa_epigenetics/data/trimmed/lane2 | grep '.fq.gz$' | cut -f 1-3 -d "_"| uniq | grep 'HA'`

do

    echo "starting sample ${sample}"

    ~/bin/bismark_v0.22.1/bismark --bowtie2 --multicore 2 -n 1 \
        --genome /data/copepods/tonsa_genome/ \
        --output_dir ~/tonsa_epigenetics/analysis/aligned/lane2/ \
        -1 ~/tonsa_epigenetics/data/trimmed/lane2/${sample}_1.fq.gz \
        -2 ~/tonsa_epigenetics/data/trimmed/lane2/${sample}_2.fq.gz \
        --rg_tag --rg_id ${sample}_lane2 --rg_sample ${sample} --un --gzip --local

    echo "done with sampe ${sample}"i

done

