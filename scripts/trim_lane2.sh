#!/bin/bash

cd ~/tonsa_epigenetics/data/trimmed/lane2

for i in $(ls /data/copepods/tonsa_epigenetics/lane2 | grep 'fq.gz' | cut -f 1-6 -d '_' | sort | uniq)

do {

 base=$(echo $i | cut -f 1-3 -d "_")

echo $base started

  java -jar /data/programs/Trimmomatic-0.36/trimmomatic-0.36.jar PE \
   -threads 2 \
    /data/copepods/tonsa_epigenetics/lane2/${i}_1.fq.gz \
    /data/copepods/tonsa_epigenetics/lane2/${i}_2.fq.gz \
    ${base}_1.fq.gz ${base}_1_unpaired.fq.gz  \
    ${base}_2.fq.gz ${base}_2_unpaired.fq.gz \
    ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:8:TRUE \
    LEADING:2 TRAILING:2 \
    SLIDINGWINDOW:4:2 \
    MINLEN:31

echo $i done

  }

done


for i in $(ls *.fq.gz | cut -f -3 -d "_" | sort | uniq)

do {

  lines=$(zcat ${i}_1.fq.gz | wc -l )

  echo ${i},$(( lines / 4 ))

} >> ~/tonsa_epigenetics/count_raw_lane2.txt

done

