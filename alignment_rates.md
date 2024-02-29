#!/bin/bash



# Because of the hacked together previous steps, need to manually calculate the alignment rates. 

## first count up total number of reads

cd ~/tonsa_epigenetics/data/trimmed

for i in $(ls *.fq.gz | cut -f -3 -d "_" | sort | uniq)

do {

    lines_l1=$(zcat ~/tonsa_epigenetics/data/trimmed/${i}_1.fq.gz | wc -l )
    lines_l2=$(zcat ~/tonsa_epigenetics/data/trimmed/lane2/${i}_1.fq.gz | wc -l )

    echo ${i},$(( lines_l1 / 4 )),$(( lines_l2 / 4 ))

} >> ~/tonsa_epigenetics/count_raw.txt

done

# then caluclate mapped

echo  sample,paired_lane1,single_l1_R1,single_l1_R2,paired_lane2,single_l2_R1,single_l2_R2 > ~/tonsa_epigenetics/count_all_aligned.txt

for sample in `ls ~/tonsa_epigenetics/analysis/aligned/lane1 | grep '.bam' | cut -f 1-3 -d '_' | sort | uniq`

do {

  l1_paired=$(samtools view -c ~/tonsa_epigenetics/analysis/aligned/lane1/${sample}_1_bismark_bt2_pe.bam)
  l2_paired=$(samtools view -c ~/tonsa_epigenetics/analysis/aligned/lane2/${sample}_1_bismark_bt2_pe.bam)
  l1_read1=$(samtools view -c ~/tonsa_epigenetics/analysis/aligned/single_end_lane1/${sample}_1.fq.gz_unmapped_reads_1_bismark_bt2.bam)
  l2_read1=$(samtools view -c ~/tonsa_epigenetics/analysis/aligned/single_end_lane2/${sample}_1.fq.gz_unmapped_reads_1_bismark_bt2.bam)
  l1_read2=$(samtools view -c ~/tonsa_epigenetics/analysis/aligned/single_end_lane1/${sample}_2.fq.gz_unmapped_reads_2_bismark_bt2.bam)
  l2_read2=$(samtools view -c ~/tonsa_epigenetics/analysis/aligned/single_end_lane2/${sample}_2.fq.gz_unmapped_reads_2_bismark_bt2.bam)

    echo ${sample},${l1_paired},${l1_read1},${l1_read2},${l2_paired},${l2_read1},${l2_read1}

} >> ~/tonsa_epigenetics/count_all_aligned.txt

done

# plot results

```r

library(ggplot2)

raw <- read.csv("~/tonsa_epigenetics/count_raw_20200318.txt", header=F)
# remember that this is just the # of reads from the forward read. Need to double
colnames(raw) <- c("sample", "raw_l1", "raw_l2")

# double, bc only one end
raw$rawcount_l1 <- raw$raw_l1*2
raw$rawcount_l2 <- raw$raw_l2*2

aligned <- read.csv("~/tonsa_epigenetics/count_all_aligned.txt", header=T)

# the columns are: aligned paired bam, and other 2 are single end.
# sum them for total
aligned$mapcount_l1 <- aligned$paired_lane1 + aligned$single_l1_R1 + aligned$single_l1_R2
aligned$mapcount_l2 <- aligned$paired_lane2 + aligned$single_l2_R1 + aligned$single_l2_R2

aligned$group <- factor(substr(aligned$sample, 1, 6))

# merge
tot <- merge(aligned, raw, by="sample")

# get total rates

tot$mapcount <- tot$mapcount_l1 + tot$mapcount_l2
tot$rawcount <- tot$rawcount_l1 + tot$rawcount_l2
tot$align_rate <- tot$mapcount/tot$rawcount

p1 <- ggplot(data=tot,aes(x= sample, y=align_rate)) +
   #geom_smooth(method = "lm")+
   geom_bar(stat="identity",color="black", aes(fill=group)) + 
   theme_bw() +
   ylab("alignment rate")+
   #stat_cor(method = "pearson") +
   ggtitle("alignment rate") +
   theme(axis.text.x = element_text(angle = 45, hjust = 1))

p2 <- ggplot(data=tot,aes(x= sample, y=mapcount)) +
   #geom_smooth(method = "lm")+
   geom_bar(stat="identity",color="black", aes(fill=group)) + 
   theme_classic() +
   ylab("mapping count")+
   #stat_cor(method = "pearson") +
   ggtitle("mapping count") +
   theme(axis.text.x = element_text(angle = 45, hjust = 1))

p3 <- ggplot(data=tot,aes(x= sample, y=rawcount)) +
   #geom_smooth(method = "lm")+
   geom_bar(stat="identity",color="black", aes(fill=group)) + 
   theme_classic() +
   ylab("raw count")+
   #stat_cor(method = "pearson") +
   ggtitle("raw count") +
   theme(axis.text.x = element_text(angle = 45, hjust = 1))

p4 <- ggplot(data=tot,aes(x= sample, y=rawcount_l1)) +
   #geom_smooth(method = "lm")+
   geom_bar(stat="identity",color="black", aes(fill=group)) +
   theme_classic() +
   ylab("raw count lane 1")+
   #stat_cor(method = "pearson") +
   ggtitle("raw count lane 1") +
   theme(axis.text.x = element_text(angle = 45, hjust = 1))
p5 <- ggplot(data=tot,aes(x= sample, y=rawcount_l2)) +
   #geom_smooth(method = "lm")+
   geom_bar(stat="identity",color="black", aes(fill=group)) +
   theme_classic() +
   ylab("raw count lane 2") +
   #stat_cor(method = "pearson") +
   ggtitle("raw count lane 2") +
   theme(axis.text.x = element_text(angle = 45, hjust = 1))

library(ggpubr)
ggsave("~/tonsa_epigenetics/figures/counts.png", ggarrange(p3, p2, p1, p4, p5, nrow=2, ncol=3, common.legend=TRUE),
  h=12, w=13)


ggsave("~/tonsa_epigenetics/figures/map_rates.png",p1,h=3, w=5)

write.table(tot, file="~/tonsa_epigenetics/mapping_rates.csv",
              row.names=F, quote=F, sep=",")

```
