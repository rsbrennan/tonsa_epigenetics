for sample in `ls ~/tonsa_epigenetics/data/trimmed | grep '.fq.gz$' | cut -f 1-3 -d "_"| uniq`

do

rg=$(echo \@RG\\tID:${sample}\\tPL:Illumina\\tPU:x\\tLB:${sample}\\tSM:${sample})


bwameth.py --reference /data/copepods/tonsa_genome/GCA_900241095.1_Aton1.0_genomic.fa \
    ~/tonsa_epigenetics/data/trimmed/${sample}_1.fq.gz \
    ~/tonsa_epigenetics/data/trimmed/${sample}_2.fq.gz \
    --read-group $rg -t 6 | \
    ~/bin/samblaster/samblaster | \
    samtools view -h -u - | \
    samtools sort -o ~/tonsa_epigenetics/analysis/aligned/${sample}.bam

done
