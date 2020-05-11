cd ~/tonsa_epigenetics/analysis/aligned/lane2/

for sample in `ls ~/tonsa_epigenetics/data/trimmed/lane2/ | grep '.fq.gz$' | cut -f 1-3 -d "_"| uniq | grep 'AA_F25'`

do

    echo "starting sample ${sample}"

        ~/bin/Bismark-0.22.3/bismark --bowtie2 --multicore 1 \
        --genome /data/copepods/tonsa_genome/ \
        --output_dir ~/tonsa_epigenetics/analysis/aligned/lane2/ \
        -1 ~/tonsa_epigenetics/data/trimmed/lane2/${sample}_1.fq.gz \
        -2 ~/tonsa_epigenetics/data/trimmed/lane2/${sample}_2.fq.gz \
        --rg_tag --rg_id ${sample} --rg_sample ${sample} --un --gzip --local --maxins 1000

    echo "done with sample ${sample}"

done

