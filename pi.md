# pi and selection.
# pi and methylation regions.

### in bash:

calculate pi.

```bash
cd ~/tonsa_epigenetics/analysis/pi

for i in $(ls ~/tonsa_genomics/data/aligned/merged | cut -f 1 -d '.'); do

  echo "Status: starting $i";

  samtools mpileup -Q 20 -B --max-depth 6000 --skip-indels -f /data/copepods/tonsa_genome/GCA_900241095.1_Aton1.0_genomic.fa ~/tonsa_genomics/data/aligned/merged/${i}.bam > ~/tonsa_genomics/analysis/pi/${i}.mpileup;

  echo "Status: $i pileup done; starting popoolation";

  perl ~/bin/popoolation_1.2.2/Variance-sliding.pl --pool-size 50 --min-qual 20 --min-coverage 20 --min-count 2 --max-coverage 3000 --min-covered-fraction 0.5 --input ~/tonsa_genomics/analysis/pi/${i}.mpileup --window-size 100 --step-size 100 --output ~/tonsa_epigenetics/analysis/pi/${i}.l2.pi --measure pi;

  rm ~/tonsa_genomics/analysis/pi/${i}.mpileup;

  echo "Status: $i popoolation done, pileup removed";

done


```

Analyze results


```r
library(tidyverse)
library(stringr)
library(ggplot2)
library(ggpubr)
library(data.table)
library(ggbeeswarm)

#locate the directory containing the files.
#dir <- "~/tonsa_genomics/analysis/pi"
dir <- "~/tonsa_epigenetics/analysis/pi"
files <- file.path(dir, list.files(dir))
files <- files[grep("params", files, invert=TRUE)]
files <- files[grep("F03", files, invert=TRUE)]
files <- files[grep("HH_F00", files, invert=TRUE)]
files <- files[grep("l2.pi", files, invert=FALSE)]

# read in files
d <- lapply(files, fread)

# rename without path, etc.
names(d) <- (str_split_fixed(files, "/", 5)[,5] %>% str_split_fixed( "[.]", 3))[,1]

# assign column names
d <- lapply(d, function(x) {
  colnames(x) <- c("gene", "position", "n_variants", "prop_covered", "pi")
  x
})

# convert pi to numeric
d <- lapply(d, function(x) {
  x$pi <- as.numeric(as.character(x$pi))
  x
})

# make snp name column
d <- lapply(d, function(x) {
  x$snp <- paste(x$gene, x$position, sep="_")
  x
})

# assign group id to each
for(i in 1:length(d)){
 d[[i]]$gp <- names(d)[i]

}

pops <- names(d)

for(i in 1:length(d)){
 print(sum(!is.na(d[[i]]$pi)))

}

# drop na's from the dataset
d2 <- vector(mode = "list", length = 16)

for(i in 1:length(d)){
 d2[[i]] <- d[[i]][!is.na(d[[i]]$pi),]
}



######################
#
# diversity in methylation regions
#
######################

# write output to look for overlap. will use bedtools.

dpi <- bind_rows(d, .id = "column_label")

dout <- dpi[(which(dpi$pi != "na")),]

dat <- read.csv("~/tonsa_epigenetics/analysis/diff_methylation/methylation_summary.txt", sep="\t")

# split snp into chr and snp
write.table(cbind(as.character(dat$CHR), as.numeric(dat$POS)-1, dat$POS, 
    dat$ah_fdr_meth, dat$ha_fdr_meth, dat$hh_fdr_meth,
    dat$ah_large, dat$ha_large, dat$hh_large),
        file="~/tonsa_epigenetics/analysis/all_meth.bed",
        quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

dout$start <- dout$position-50
dout$stop <- dout$position+51


# add pi measure to this. and group
write.table(cbind(dout$gene, dout$start, dout$stop, dout[,4:8]),
        file="~/tonsa_epigenetics/analysis/pi.bed",
        quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

```


look for intersection

```bash
cd ~/tonsa_epigenetics/analysis/
sort -k1,1 -k2,2n pi.bed > pi.sorted.bed
sort -k1,1 -k2,2n all_meth.bed > meth.sorted.bed
#
bedtools intersect -wa -wb \
    -a pi.sorted.bed \
    -b meth.sorted.bed \
    -sorted > ~/tonsa_epigenetics/analysis/pi.overlp.bed


```


```r
###########
## Read in results
###########
dat <- read.csv("~/tonsa_epigenetics/analysis/pi.overlp.bed", header=FALSE, sep="\t")

colnames(dat) <- c("CHR", "start", "stop","nsnp", "cov", "pi",
                    "pi_snp", "gp", "gene", "snp_start", "snp_stop",
                     "ah_fdr", "ha_fdr", "hh_fdr",
                     "ah_large", "ha_large", "hh_large")
dat <- dat[grep("F00",dat$gp, invert=T),]
dat <- dat[grep("AA",dat$gp, invert=T),]

dat$trt <- str_split_fixed(dat$gp, "_", 2)[,1]
dat$trt <- as.factor(dat$trt)

################
###
### delta pi. to look for enrichment
###
################

# the data frame has each replicate of all treatments. so rep1-4 for each.

# want to check if regions that have significant methylation change have
    # lower or high diversity. 

# add in direction of change:
datdir <- read.csv("~/tonsa_epigenetics/analysis/diff_methylation/methylation_summary.txt", sep="\t")

# add mthylation snp id
dat$snp_meth <- paste(dat$CHR, dat$snp_stop, sep=":")

merged <- merge(dat, datdir, by.x="snp_meth", by.y="SNP")
dat <- merged
# add p-val column
dat$pval <- NA
dat$pval[which(dat$trt == "AH")] <- dat$ah_fdr[which(dat$trt == "AH")]
dat$pval[which(dat$trt == "HA")] <- dat$ha_fdr[which(dat$trt == "HA")]
dat$pval[which(dat$trt == "HH")] <- dat$hh_fdr[which(dat$trt == "HH")]

# need to add a significance column that matches each treatment:
dat$ah_sig <- ifelse(dat$ah_fdr < 0.05 & dat$ah_large.x == TRUE, TRUE, FALSE)
dat$ha_sig <- ifelse(dat$ha_fdr < 0.05 & dat$ha_large.x == TRUE, TRUE, FALSE)
dat$hh_sig <- ifelse(dat$hh_fdr < 0.05 & dat$hh_large.x == TRUE, TRUE, FALSE)
dat$sig <- NA
dat$sig[which(dat$trt == "AH")] <- dat$ah_sig[which(dat$trt == "AH")]
dat$sig[which(dat$trt == "HA")] <- dat$ha_sig[which(dat$trt == "HA")]
dat$sig[which(dat$trt == "HH")] <- dat$hh_sig[which(dat$trt == "HH")]

# add direction of change:
dat$direction <- NA
dat$AA_mean_meth <- rowMeans(cbind(dat$aa_F00_mean_meth, dat$aa_F25_mean_meth))
dat$direction[which(dat$trt == "AH")] <- ifelse((dat$aa_F25_mean_meth[which(dat$trt == "AH")] - dat$ah_F25_mean_meth[which(dat$trt == "AH")]) < 0, "increase", "decrease")
dat$direction[which(dat$trt == "HA")] <- ifelse((dat$aa_F25_mean_meth[which(dat$trt == "HA")] - dat$ah_F25_mean_meth[which(dat$trt == "HA")]) < 0, "increase", "decrease")
dat$direction[which(dat$trt == "HH")] <- ifelse((dat$aa_F25_mean_meth[which(dat$trt == "HH")] - dat$ah_F25_mean_meth[which(dat$trt == "HH")]) < 0, "increase", "decrease")



# first, scatter plot of pi with -log10 methylation significance

dat$pval_log <- -log10(dat$pval)

ggplot(dat, aes(x=pi, y=pval_log)) +
        geom_point(alpha=0.3)+
          geom_smooth(method='lm') +
          facet_wrap(~trt)

x <- dat$pi[which(dat$trt == "HH")]
y <- dat$pval_log[which(dat$trt == "HH")]

summary(lm(y ~ x))

x <- dat$pi[which(dat$trt == "AH")]
y <- dat$pval_log[which(dat$trt == "AH")]

summary(lm(y ~ x))

# positive relationship between pi and methylation significance. 

######
# compare the distributions of the different groups.


# dat$ah_sig <- ifelse(dat$ah_fdr < 0.05, TRUE, FALSE)
# dat$ha_sig <- ifelse(dat$ha_fdr < 0.05, TRUE, FALSE)
# dat$hh_sig <- ifelse(dat$hh_fdr < 0.05, TRUE, FALSE)

ah <- dat[dat$trt == "AH",]
ha <- dat[dat$trt == "HA",]
hh <- dat[dat$trt == "HH",]

ks.test((ah$pi[ah$ah_sig == TRUE]), 
        (ah$pi[ah$ah_sig == FALSE]))
# D = 0.23773, p-value = 0.06763
ks.test((ha$pi[ha$ha_sig == TRUE]), 
        (ha$pi[ha$ha_sig == FALSE]))
# D = 0.10828, p-value = 0.2434
ks.test((hh$pi[hh$hh_sig == TRUE]), 
        (hh$pi[hh$hh_sig == FALSE]))
# D = 0.11041, p-value = 3.376e-05

p.adjust(c(0.06763, 0.2434, 3.376e-05), method="bonferroni")
# [1] 0.20289000 0.73020000 0.00010128

stderr <- function(x) sd(x)/sqrt(length(x))
dat %>%
    group_by(trt,sig) %>%
    dplyr::summarize(median = median(pi, na.rm=TRUE),
                     mean = mean(pi, na.rm=TRUE), 
                     stderr = stderr(pi),
                    count = n())

#   trt   sig    median   mean    stderr count
#   <fct> <lgl>   <dbl>  <dbl>     <dbl> <int>
# 1 AH    FALSE 0.0140  0.0165 0.0000880 23931
# 2 AH    TRUE  0.00699 0.0169 0.00312      30
# 3 HA    FALSE 0.0145  0.0169 0.0000756 31892
# 4 HA    TRUE  0.0154  0.0205 0.00186      90
# 5 HH    FALSE 0.0143  0.0171 0.0000697 39032
# 6 HH    TRUE  0.0150  0.0196 0.000726    456


# change labels to match MS:
dat$trt <- as.character(dat$trt)
dat$trt[which(dat$trt =="AH")] <- c("Acidification")
dat$trt[which(dat$trt =="HA")] <- c("Warming")
dat$trt[which(dat$trt =="HH")] <- c("OWA")
dat$trt <- factor(dat$trt, levels=c("Acidification", "Warming", "OWA"))

p <- ggplot(dat, aes(x=gp, y=pi, fill=sig)) +
        geom_boxplot(alpha=0.3, color="black") +
        #geom_quasirandom(alpha=0.3) +
        theme_classic() +
        scale_fill_manual(values=c("grey50", "firebrick3")) +
        stat_summary(fun="mean",geom="point",size=3, shape=18,  position = position_dodge(width = 0.75), color="dodgerblue3") +
        facet_wrap(~trt, nrow = 1, scales = "free_x") +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
        xlab("")

ggsave("~/tonsa_epigenetics/figures/pi_meth_all.pdf",
         p,
        width=8, height=3)




dout <- dat %>%
    group_by(trt,gp,sig) %>%
    dplyr::summarize(median = median(pi, na.rm=TRUE), 
                    count = n())
data.frame(dout)

## !!!!!! n is very very low for sig windows for 
          # AH (< 10) and HA (~20). HH is 100-150 for each group

#   trt          gp   sig      median count
# 1   AH AH_F25_Rep1 FALSE 0.014044980  6287
# 2   AH AH_F25_Rep1  TRUE 0.006437468     9
# 3   AH AH_F25_Rep2 FALSE 0.013641682  6850
# 4   AH AH_F25_Rep2  TRUE 0.006988925     7
# 5   AH AH_F25_Rep3 FALSE 0.015009461  4729
# 6   AH AH_F25_Rep3  TRUE 0.008146615     6
# 7   AH AH_F25_Rep4 FALSE 0.014005072  6065
# 8   AH AH_F25_Rep4  TRUE 0.008218297     8
# 9   HA HA_F25_Rep1 FALSE 0.014555433  6953
# 10  HA HA_F25_Rep1  TRUE 0.016972731    20
# 11  HA HA_F25_Rep2 FALSE 0.014730853  9298
# 12  HA HA_F25_Rep2  TRUE 0.015384985    25
# 13  HA HA_F25_Rep3 FALSE 0.014148585  6801
# 14  HA HA_F25_Rep3  TRUE 0.012854977    17
# 15  HA HA_F25_Rep4 FALSE 0.014294349  8840
# 16  HA HA_F25_Rep4  TRUE 0.014844177    28
# 17  HH HH_F25_Rep1 FALSE 0.014467837  9044
# 18  HH HH_F25_Rep1  TRUE 0.014731648   110
# 19  HH HH_F25_Rep2 FALSE 0.014210314  7727
# 20  HH HH_F25_Rep2  TRUE 0.014589262    99
# 21  HH HH_F25_Rep3 FALSE 0.014437084 13313
# 22  HH HH_F25_Rep3  TRUE 0.015113608   141
# 23  HH HH_F25_Rep4 FALSE 0.014255022  8948
# 24  HH HH_F25_Rep4  TRUE 0.014908401   106

```


