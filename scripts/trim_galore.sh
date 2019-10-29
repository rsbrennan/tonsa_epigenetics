cd ~/tonsa_epigenetics/data/trimmed

indir=/data/copepods/tonsa_epigenetics/
outdir=~/tonsa_epigenetics/data/trimmed/
for i in $(ls /data/copepods/tonsa_epigenetics/ | grep 'fq.gz' | cut -f 1-6 -d '_' | sort | uniq)

do {

 base=$(echo $i | cut -f 1-3 -d "_")

echo $i started

~/bin/TrimGalore-0.6.0/trim_galore --rrbs --paired --non_directional\
    ${indir}${i}_1.fq.gz \
    ${indir}${i}_2.fq.gz \
    --output_dir ${outdir} \
    --basename ${base}

echo $i done

  }

done
