cd ~/tonsa_epigenetics/analysis/aligned/lane2/

for sample in `ls ~/tonsa_epigenetics/data/trimmed/lane2 | grep '.fq.gz$' | cut -f 1-3 -d "_"| uniq | grep 'AA'`

do

    echo "starting sample ${sample} read 1"

    ~/bin/bismark_v0.22.1/bismark --bowtie2 --multicore 2 -n 1 \
        --genome /data/copepods/tonsa_genome/ \
        --output_dir ~/tonsa_epigenetics/analysis/aligned/lane2/ \
        ~/tonsa_epigenetics/analysis/aligned/lane2/${sample}_1.fq.gz_unmapped_reads_1.fq.gz \
        --rg_tag --rg_id ${sample}_lane2 --rg_sample ${sample} --un --gzip --local

    echo "done with sample ${sample} read 1"

    echo "starting sample ${sample} read 2"

    ~/bin/bismark_v0.22.1/bismark --bowtie2 --multicore 2 -n 1 \
        --genome /data/copepods/tonsa_genome/ \
        --output_dir ~/tonsa_epigenetics/analysis/aligned/lane2/ \
        ~/tonsa_epigenetics/analysis/aligned/lane2/${sample}_2.fq.gz_unmapped_reads_2.fq.gz \
        --rg_tag --rg_id ${sample}_lane2 --rg_sample ${sample} --un --local --pbat

    echo "done with sample ${sample} read 2"

done
