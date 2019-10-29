cd ~/tonsa_epigenetics/analysis/aligned

#for sample in `ls ~/tonsa_epigenetics/analysis/aligned| grep 'unmapped' | cut -f 1-3 -d "_"| uniq`

#do


sample=AA_F00_1

    echo "starting sample ${sample} read1"

    ~/bin/bismark_v0.22.1/bismark --bowtie2 --multicore 2 -n 1 \
        --genome /data/copepods/tonsa_genome/ \
        --output_dir ~/tonsa_epigenetics/analysis/aligned/ \
        ~/tonsa_epigenetics/analysis/aligned/${sample}_1.fq.gz_unmapped_reads_1.fq.gz \
        --rg_tag --rg_id ${sample} --rg_sample ${sample} --gzip --non_directional --upto 1000000
        #--basename ${sample}_unmap_1

    echo "done with sampe ${sample} read1"

#done

# and read 2
    echo "starting sample ${sample} read2"

    #~/bin/bismark_v0.22.1/bismark --bowtie2 --multicore 2 -n 1 \
    #    --genome /data/copepods/tonsa_genome/ \
    #    --output_dir ~/tonsa_epigenetics/analysis/aligned/ \
    #    ~/tonsa_epigenetics/analysis/aligned/${sample}_2.fq.gz_unmapped_reads_2.fq.gz\
    #    --rg_tag --rg_id ${sample} --rg_sample ${sample} --gzip --pbat
        #--basename ${sample}_unmap_2

    echo "done with sampe ${sample} read2"
