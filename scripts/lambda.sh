cd ~/tonsa_epigenetics/analysis/lambda/aligned

for sample in `ls ~/tonsa_epigenetics/data/trimmed | grep '.fq.gz$' | cut -f 1-3 -d "_"| uniq`

do

    echo "starting sample ${sample}"

        ~/bin/Bismark-0.22.3/bismark --bowtie2 --multicore 1 \
        --genome_folder ~/tonsa_epigenetics/analysis/lambda/ \
        --output_dir ~/tonsa_epigenetics/analysis/lambda/aligned \
        -1 ~/tonsa_epigenetics/data/trimmed/${sample}_1.fq.gz \
        -2 ~/tonsa_epigenetics/data/trimmed/${sample}_2.fq.gz \
        --rg_tag --rg_id ${sample} --rg_sample ${sample} --un --gzip --local --maxins 1000

     echo "finished sample ${sample}"

done

