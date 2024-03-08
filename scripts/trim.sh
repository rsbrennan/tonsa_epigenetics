#!/bin/bash

## trim lane one

cd ~/tonsa_epigenetics/data/trimmed

for i in $(ls /data/copepods/tonsa_epigenetics/ | grep 'fq.gz' | cut -f 1-6 -d '_' | sort | uniq)

do {

 base=$(echo $i | cut -f 1-3 -d "_")

echo $base started

  java -jar /data/programs/Trimmomatic-0.36/trimmomatic-0.36.jar PE \
   -threads 8 \
    /data/copepods/tonsa_epigenetics/${i}_1.fq.gz \
    /data/copepods/tonsa_epigenetics/${i}_2.fq.gz \
    ${base}_1.fq.gz ${base}_1_unpaired.fq.gz  \
    ${base}_2.fq.gz ${base}_2_unpaired.fq.gz \
    ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:8:TRUE \
    SLIDINGWINDOW:5:20 \
    MINLEN:31

echo $i done

  }

done

## trim lane two

cd ~/tonsa_epigenetics/data/trimmed/lane2

for i in $(ls /data/copepods/tonsa_epigenetics/lane2/ | grep 'fq.gz' | cut -f 1-6 -d '_' | sort | uniq)

do {

 base=$(echo $i | cut -f 1-3 -d "_")

echo $base started

  java -jar /data/programs/Trimmomatic-0.36/trimmomatic-0.36.jar PE \
   -threads 8 \
    /data/copepods/tonsa_epigenetics/lane2/${i}_1.fq.gz \
    /data/copepods/tonsa_epigenetics/lane2/${i}_2.fq.gz \
    ${base}_1.fq.gz ${base}_1_unpaired.fq.gz  \
    ${base}_2.fq.gz ${base}_2_unpaired.fq.gz \
    ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:8:TRUE \
    SLIDINGWINDOW:5:20 \
    MINLEN:31

echo $i done

  }

done

