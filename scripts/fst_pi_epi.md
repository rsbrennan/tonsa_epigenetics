# pi and selection.
# pi and methylation regions.


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


By this point, have pi and methylation overlap. Need to get pi and fst overlap. 


```bash

# intersect fst with pi
bedtools intersect \
    -a ~/tonsa_epigenetics/analysis/pi.sorted.bed \
    -b ~/tonsa_epigenetics/analysis/windows_fixed_1500kb_counts_reps_fst.txt -wa -wb   \
    > ~/tonsa_epigenetics/analysis/fst_pi_intersection.txt
        # 26432

# actually, want it with this: ~/tonsa_epigenetics/analysis/windows_fixed_${win_size}kb_reps_fst_epi.txt

bedtools intersect \
    -a ~/tonsa_epigenetics/analysis/pi.sorted.bed \
    -b ~/tonsa_epigenetics/analysis/windows_fixed_1500kb_reps_fst_epi.txt -wa -wb   \
    > ~/tonsa_epigenetics/analysis/fst_methylation_pi_intersection.txt
        # 171668
cat ~/tonsa_epigenetics/analysis/fst_methylation_pi_intersection.txt | wc -l
# 166987
```


now I have a table with Fst, methylation, and pi windows. 

With this, I run the same analysis summarizing fst/methylation. but I also calculate pi. Then I make bins for pi and ask if the patterns are consistent across bins.



```r
###########
## Read in results
###########
dat_in <- read.csv("~/tonsa_epigenetics/analysis/fst_methylation_pi_intersection.txt", header=FALSE, sep="\t")
dim(dat_in)
#[1] 166987     49

colnames(dat_in) <- c(

    #pi labels: should be 8.
    "CHR_pi", "start_pi", "stop_pi",
    "nsnp_pi", "cov", "pi",
    "pi_snp", "gp_pi", 


    # methylation labels
    "chr","start", "stop","SNP","aa_F00_mean_meth",
 "aa_F25_mean_meth","ah_F25_mean_meth","ha_F25_mean_meth","hh_F25_mean_meth",
    "ha_pval_meth","ha_fdr_meth","ah_pval_meth","ah_fdr_meth",
    "hh_pval_meth","hh_fdr_meth","ah_large","ha_large",
    "hh_large",
    "ah_delta_meth_F0","ha_delta_meth_F0","hh_delta_meth_F0",
    "ah_delta_meth_F25","ha_delta_meth_F25","hh_delta_meth_F25",
    # snp labels
    "chr_win", "start_win", "end_win", "window_id", "AH_F25_Rep1_fst", "AH_F25_Rep2_fst", "AH_F25_Rep3_fst", "AH_F25_Rep4_fst", "HA_F25_Rep1_fst", "HA_F25_Rep2_fst", "HA_F25_Rep3_fst", "HA_F25_Rep4_fst", "HH_F25_Rep1_fst", "HH_F25_Rep2_fst", "HH_F25_Rep3_fst", "HH_F25_Rep4_fst", "snp_count")

# add pi window ids:

length(unique(dat_in$pi_snp))
# 2835

length(unique(dat_in$window_id))
#979
# add number of methylation sites per window:

summary_dat_pi <- dat_in %>%
  group_by(pi_snp) %>%
  summarise(
    mean_pi = mean(pi),
    window_id= unique(window_id)
    ) %>%
  group_by(window_id) %>%
    summarise(
    mean_pi_fstwin = mean(mean_pi),
    window_id= unique(window_id),
    pi_win_count = length(unique(pi_snp)), 
    )

summary_dat <- dat_in %>%
  group_by(window_id) %>%
  summarise(
    min_ha_fdr_meth = min(ha_fdr_meth, na.rm = TRUE),
    min_ah_fdr_meth = min(ah_fdr_meth, na.rm = TRUE),
    min_hh_fdr_meth = min(hh_fdr_meth, na.rm = TRUE),
    meth_count = n(), 
  )


# then subset snp data to only relevant cols
snpdat <- dat_in[,33:ncol(dat_in)] %>% distinct(window_id, .keep_all = TRUE)
nrow(snpdat)
# 979

# join the three dfs:
alldat1 <- left_join(snpdat,summary_dat, by = "window_id")
nrow(alldat1)
alldat <- left_join(alldat1, summary_dat_pi, by = "window_id")
nrow(alldat)
# 979
head(alldat)



  new_dat <- alldat
  new_dat$ah_sig <- c("non-sig")
  new_dat$ha_sig <- c("non-sig")
  new_dat$hh_sig <- c("non-sig")
  new_dat$ah_sig[which(new_dat$min_ah_fdr_meth< 0.05)] <- c("sig")
  new_dat$ha_sig[which(new_dat$min_ha_fdr_meth< 0.05)] <- c("sig")
  new_dat$hh_sig[which(new_dat$min_hh_fdr_meth< 0.05)] <- c("sig")

dat_sub <- new_dat[which(new_dat$snp_count > 4),]
dat_sub <- new_dat[which(new_dat$snp_count > 4 & new_dat$meth_count > 4),]
nrow(dat_sub)
# 867

sum(dat_sub$min_ah_fdr_meth< 0.05)
# 53
sum(dat_sub$min_ha_fdr_meth< 0.05)
# 82
sum(dat_sub$min_hh_fdr_meth< 0.05)
# 158




######### ---------------
# make bins for pi

head(dat_sub)

table(dat_sub$pi_win_count)
hist(dat_sub$mean_pi_fstwin, breaks=30)

range(dat_sub$mean_pi_fstwin)

# make bins:
breaks <- c(0, 0.005, 0.01, 0.015, 0.02, 0.03, 0.065)
# add col with bin assignments 
dat_sub$pi_bins <- cut(dat_sub$mean_pi_fstwin, 
                      breaks = breaks,
                      include.lowest = TRUE,
                      labels = paste(breaks[-length(breaks)], 
                                   breaks[-1], 
                                   sep = " - "))
table(dat_sub$pi_bins)
#   0 - 0.005 0.005 - 0.01 0.01 - 0.015 0.015 - 0.02  0.02 - 0.03 0.03 - 0.065
#         103          165          173          173          170           83



# make long format for plotting

repnum <- nrow(dat_sub)
df_p <- data.frame(
        treatment = c(rep(rep("Acidification", repnum),4),rep(rep("Warming", repnum),4),rep(rep("OWA", repnum),4)),

        replicate = c(rep(c(rep("Rep1", repnum),
                                            rep("Rep2", repnum),
                                            rep("Rep3", repnum), 
                                            rep("Rep4", repnum)),
                                          3)),
        
        fst = c(dat_sub$AH_F25_Rep1_fst, dat_sub$AH_F25_Rep2_fst,dat_sub$AH_F25_Rep3_fst,dat_sub$AH_F25_Rep4_fst,
                        dat_sub$HA_F25_Rep1_fst, dat_sub$HA_F25_Rep2_fst,dat_sub$HA_F25_Rep3_fst,dat_sub$HA_F25_Rep4_fst,
                        dat_sub$HH_F25_Rep1_fst, dat_sub$HH_F25_Rep2_fst,dat_sub$HH_F25_Rep3_fst,dat_sub$HH_F25_Rep4_fst),
        significant =c(rep(dat_sub$ah_sig, 4),
                                     rep(dat_sub$ha_sig, 4),
                                     rep(dat_sub$hh_sig, 4)),
        pi_bin = dat_sub$pi_bins,
        pi_mean = dat_sub$mean_pi_fstwin
                            ) 

df_p$treatment <- factor(df_p$treatment, levels = c("Acidification", "Warming", "OWA"))



p <- ggplot(df_p, aes(x = significant, y = fst)) +
  theme_classic()  + 
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_summary(aes(group=treatment), fun = mean, geom = "line", 
              linewidth = 0.8) + 
  stat_summary(aes(group=treatment, shape=treatment, fill=treatment), fun = mean, geom = "point", 
              size = 5) +
  stat_summary(aes(group=replicate, shape=treatment, fill=treatment), fun = mean, geom = "point", 
              size = 2, alpha=0.5) +
  facet_grid(pi_bin ~ treatment) +  # This creates a grid with pi_bin rows and treatment columns
  scale_fill_manual(values=c("#F2AD00", "#00A08A", "#CC3333")) +
  scale_shape_manual(values=c(22,23,24)) +
  theme(legend.position="none") +
  xlab("Significant methylation change") +
  ylab("Fst")


ggsave( p,
          file="~/tonsa_epigenetics/figures/window_reps_sig_FST_PI_1500kb.pdf",
          h=8, w=4.5)

# get summary of values:

mean_fst_df <- df_p %>%
  group_by(treatment, replicate, significant,pi_bin) %>%
  summarise(mean_fst = mean(fst, na.rm = TRUE),
                    median_fst = median(fst, na.rm = TRUE),
                    std_fst = sd(fst, na.rm=TRUE), .groups = 'drop')
as.data.frame(mean_fst_df)

mean_fst_df <- df_p %>%
  group_by(treatment, significant,pi_bin) %>%
  summarise(mean_fst = mean(fst, na.rm = TRUE),
                    median_fst = median(fst, na.rm = TRUE),
                    std_fst = sd(fst, na.rm=TRUE), .groups = 'drop')
as.data.frame(mean_fst_df)

```

