# methytlation window analysis

Looking for relationship of genomic regions with methylation versus genetic change. 

### Do regions that show significant methylation change also show diff allele freq divergence?

Basic approach: 

Make windows around both the epialleles, and the SNPs, then can get means, etc using bedtools groupby

Calculate mean fst in windows, then ask if those windows have sig change in methylation or not.


```r


library(poolfstat)
library(scales)
setwd("~/tonsa_genomics/analysis/")


dat <- read.table("~/tonsa_genomics/analysis/filtered_variants.txt", header=TRUE, stringsAsFactors=FALSE, nrow=2)

pops <- colnames(dat)[11:ncol(dat)]

pdat <- popsync2pooldata(sync.file="variants.sync", poolsizes = rep(100, 28),
    poolnames = pops, min.maf=0, min.rc=0,
     max.cov.per.pool = 1e+07)

fst <- computePairwiseFSTmatrix(pdat, method = "Anova",
  min.cov.per.pool = -1, max.cov.per.pool = 1e+07, min.maf = -1,
  output.snp.values = TRUE)

global.fsts <- computeFST(pdat, method = "Anova", snp.index = NA)


dim(fst$PairwiseSnpFST) 

fst_sub <- as.data.frame(fst$PairwiseSnpFST[,grep("AA_F25",colnames(fst$PairwiseSnpFST))])


all_indiv <- fst_sub[,colnames(fst_sub)[
                      grep(
                        "AA_F25_Rep1_vs_AH_F25_Rep1|AA_F25_Rep2_vs_AH_F25_Rep2|AA_F25_Rep3_vs_AH_F25_Rep3|AA_F25_Rep4_vs_AH_F25_Rep4|AA_F25_Rep1_vs_HA_F25_Rep1|AA_F25_Rep2_vs_HA_F25_Rep2|AA_F25_Rep3_vs_HA_F25_Rep3|AA_F25_Rep4_vs_HA_F25_Rep4|AA_F25_Rep1_vs_HH_F25_Rep1|AA_F25_Rep2_vs_HH_F25_Rep2|AA_F25_Rep3_vs_HH_F25_Rep3|AA_F25_Rep4_vs_HH_F25_Rep4",
                      colnames(fst_sub))]]

desired_order <- c(
  "AA_F25_Rep1_vs_AH_F25_Rep1", "AA_F25_Rep2_vs_AH_F25_Rep2", "AA_F25_Rep3_vs_AH_F25_Rep3", "AA_F25_Rep4_vs_AH_F25_Rep4",
  "AA_F25_Rep1_vs_HA_F25_Rep1", "AA_F25_Rep2_vs_HA_F25_Rep2", "AA_F25_Rep3_vs_HA_F25_Rep3", "AA_F25_Rep4_vs_HA_F25_Rep4",
  "AA_F25_Rep1_vs_HH_F25_Rep1", "AA_F25_Rep2_vs_HH_F25_Rep2", "AA_F25_Rep3_vs_HH_F25_Rep3", "AA_F25_Rep4_vs_HH_F25_Rep4"
)

desired_nms <- substr(desired_order, start = 16, stop = 26)

all_indiv <- all_indiv %>% select(desired_order)

all_indiv[all_indiv < 0] <- 0
all_indiv[is.na(all_indiv)] <- 0 # bc these are invariant
colnames(all_indiv) <- paste(desired_nms, "fst", sep="_")

# Define a regular expression pattern

# then read in variant ids to make bed file.
dat <- read.table("~/tonsa_genomics/analysis/filtered_variants.txt", header=TRUE, stringsAsFactors=FALSE)

# Create a new dataframe with three columns from 'dat'
dfout <- data.frame(
  Chrom = dat$Chrom,
  start = dat$Position - 1,
  end = dat$Position
)

# Add all columns from 'all_indiv' to 'df'
dfout <- cbind(dfout, all_indiv)

write.table(file="~/tonsa_genomics/analysis/fst_reps_all.bed", dfout, sep="\t", 
        quote=F, row.names=F, col.names=F)




```

```bash

cat ~/tonsa_genomics/analysis/fst_reps_all.bed \
 |sort -k1,1 -k2,2n > ~/tonsa_epigenetics/analysis/fst_reps_all.sort.bed

```


```r
# calc methylation differences
dat <- read.table("~/tonsa_epigenetics/analysis/diff_methylation/methylation_summary.txt", header=T)

dat$ah_f0 <- dat$aa_F00_mean_meth - dat$ah_F25_mean_meth
dat$ha_f0 <- dat$aa_F00_mean_meth - dat$ha_F25_mean_meth
dat$hh_f0 <- dat$aa_F00_mean_meth - dat$hh_F25_mean_meth
dat$ah_f25 <- dat$aa_F25_mean_meth - dat$ah_F25_mean_meth
dat$ha_f25 <- dat$aa_F25_mean_meth - dat$ha_F25_mean_meth
dat$hh_f25 <- dat$aa_F25_mean_meth - dat$hh_F25_mean_meth

write.table(dat, file="~/tonsa_epigenetics/analysis/diff_methylation/methylation_summary_withDiffs.txt", row.names=F, col.names=T, quote=F, sep="\t")

```

run windowed analysis. 

```bash
# ~/tonsa_epigenetics/analysis/fst_reps_all.sort.bed
    # "Chrom", "start", "end", "AH_F25_Rep1_fst", "AH_F25_Rep2_fst", "AH_F25_Rep3_fst", "AH_F25_Rep4_fst", "HA_F25_Rep1_fst", "HA_F25_Rep2_fst", "HA_F25_Rep3_fst", "HA_F25_Rep4_fst", "HH_F25_Rep1_fst", "HH_F25_Rep2_fst", "HH_F25_Rep3_fst", "HH_F25_Rep4_fst",

cat ~/tonsa_epigenetics/analysis/diff_methylation/methylation_summary_withDiffs.txt | tail -n +2 | awk '{printf "%s\t%s\t%s\t", $1, $2-1, $2; for (i=3; i<=NF; i++) printf "%s%s", $i, (i==NF ? "\n" : "\t")}' |sort -k1,1 -k2,2n > ~/tonsa_epigenetics/analysis/gene_level_analysis/epi.bed

### they are from ~/tonsa_epigenetics/analysis/diff_methylation/methylation_summary.txt but with the 6 methylation differences added to the end

#######################
# fixed window size:
#######################

# make a loop to cycle through a couple of window sizes
for win_size in 1500 5000; do

    # make the windows
bedtools makewindows -g <(cut -f 1,2 /data/copepods/tonsa_genome/GCA_900241095.1_Aton1.0_genomic.fa.fai) \
    -w $win_size \
    -s $win_size  \
    -i srcwinnum | sort -k1,1 -k2,2n  > ~/tonsa_epigenetics/analysis/windows_${win_size}kb.bed

# Join Fst values and the 'windows.bed' file
    bedtools intersect \
    -a ~/tonsa_epigenetics/analysis/windows_${win_size}kb.bed \
    -b ~/tonsa_epigenetics/analysis/fst_reps_all.sort.bed -wa -wb | sort -k4,4  \
    > ~/tonsa_epigenetics/analysis/windows_reps_fst_${win_size}kb.tab

# Run bedtools groupby command to obtain average values of Fst for each of the windows above
bedtools groupby -i ~/tonsa_epigenetics/analysis/windows_reps_fst_${win_size}kb.tab \
    -g 1,2,3,4 \
    -c 8,9,10,11,12,13,14,15,16,17,18,19 \
    -o mean > ~/tonsa_epigenetics/analysis/windows_mean_reps_fst_${win_size}kb.txt
    # 14399

# get counts for each window. 
    bedtools groupby -i ~/tonsa_epigenetics/analysis/windows_reps_fst_${win_size}kb.tab \
    -g 1,2,3,4 \
    -c 8 \
    -o count  > ~/tonsa_epigenetics/analysis/windows_freq_fst_${win_size}kb.txt

# join the counts with the means
    cut -f 5 ~/tonsa_epigenetics/analysis/windows_freq_fst_${win_size}kb.txt | paste ~/tonsa_epigenetics/analysis/windows_mean_reps_fst_${win_size}kb.txt - | sort -k1,1 -k2,2n  > ~/tonsa_epigenetics/analysis/windows_fixed_${win_size}kb_counts_reps_fst.txt
    # 14399


#   cut -f 5 ~/tonsa_epigenetics/analysis/windows_freq_fst_${win_size}kb.txt | paste ~/tonsa_epigenetics/analysis/windows_mean_fst_${win_size}kb.txt - | sort -k1,1 -k2,2n  > ~/tonsa_epigenetics/analysis/windows_fixed_${win_size}kb_counts_reps_fst.txt
    # 14399



######
# we now have fst in each window. Then get intersect this with the epigenetics.

# intersect fst with epigenetics
    bedtools intersect \
    -a ~/tonsa_epigenetics/analysis/gene_level_analysis/epi.bed \
    -b ~/tonsa_epigenetics/analysis/windows_fixed_${win_size}kb_counts_reps_fst.txt -wa -wb   \
    > ~/tonsa_epigenetics/analysis/windows_fixed_${win_size}kb_reps_fst_epi.txt
        # 26432
#######################

done

```

# plot fst based on window significance. by replicate



```r
library(tidyverse)

dat_in <- read.table("~/tonsa_epigenetics/analysis/windows_fixed_1500kb_reps_fst_epi.txt", header=F) 
ncol(dat_in)
# 41 columns

colnames(dat_in) <- c("chr","start", "stop","SNP","aa_F00_mean_meth",
 "aa_F25_mean_meth","ah_F25_mean_meth","ha_F25_mean_meth","hh_F25_mean_meth",
    "ha_pval_meth","ha_fdr_meth","ah_pval_meth","ah_fdr_meth",
    "hh_pval_meth","hh_fdr_meth","ah_large","ha_large",
    "hh_large",
    "ah_delta_meth_F0","ha_delta_meth_F0","hh_delta_meth_F0",
    "ah_delta_meth_F25","ha_delta_meth_F25","hh_delta_meth_F25",
    # snp labels
    "chr_win", "start_win", "end_win", "window_id", "AH_F25_Rep1_fst", "AH_F25_Rep2_fst", "AH_F25_Rep3_fst", "AH_F25_Rep4_fst", "HA_F25_Rep1_fst", "HA_F25_Rep2_fst", "HA_F25_Rep3_fst", "HA_F25_Rep4_fst", "HH_F25_Rep1_fst", "HH_F25_Rep2_fst", "HH_F25_Rep3_fst", "HH_F25_Rep4_fst", "snp_count")

length(unique(dat_in$window_id))
# add number of methylation sites per window:

summary_dat <- dat_in %>%
  group_by(window_id) %>%
  summarise(
    min_ha_fdr_meth = min(ha_fdr_meth, na.rm = TRUE),
    min_ah_fdr_meth = min(ah_fdr_meth, na.rm = TRUE),
    min_hh_fdr_meth = min(hh_fdr_meth, na.rm = TRUE),
    meth_count = n(), 
  )

# new df removing the methylation sites. keeping only row per snp window
snpdat <- dat_in[,25:ncol(dat_in)] %>% distinct(window_id, .keep_all = TRUE)
nrow(snpdat)
# 1440

# join the two dfs:
alldat <- left_join(snpdat,summary_dat, by = "window_id")
nrow(alldat)
#1440

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
# 910

sum(dat_sub$min_ah_fdr_meth< 0.05)
# 71
sum(dat_sub$min_ha_fdr_meth< 0.05)
# 99
sum(dat_sub$min_hh_fdr_meth< 0.05)
# 187


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
                                     rep(dat_sub$hh_sig, 4))
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
  facet_wrap(vars(treatment), nrow=1) +
  scale_fill_manual(values=c("#F2AD00", "#00A08A", "#CC3333")) +
  scale_shape_manual(values=c(22,23,24)) +
  theme(legend.position="none") +
  xlab("Significant methylation change")


ggsave( p,
          file="~/tonsa_epigenetics/figures/window_reps_sig_FST_1500kb.pdf",
          h=3, w=4.5)

# get summary of values:

mean_fst_df <- df_p %>%
  group_by(treatment, replicate, significant) %>%
  summarise(mean_fst = mean(fst, na.rm = TRUE),
                    median_fst = median(fst, na.rm = TRUE),
                    std_fst = sd(fst, na.rm=TRUE), .groups = 'drop')
as.data.frame(mean_fst_df)

#        treatment replicate significant    mean_fst  median_fst     std_fst
#1  Acidification      Rep1     non-sig 0.007670427 0.002837477 0.010847550
#2  Acidification      Rep1         sig 0.003500428 0.002300155 0.004106526
#3  Acidification      Rep2     non-sig 0.008127625 0.003358448 0.011266150
#4  Acidification      Rep2         sig 0.004587838 0.002507372 0.006960764
#5  Acidification      Rep3     non-sig 0.008358349 0.003220845 0.011491206
#6  Acidification      Rep3         sig 0.003638063 0.001962236 0.004959680
#7  Acidification      Rep4     non-sig 0.007479710 0.002788938 0.010232173
#8  Acidification      Rep4         sig 0.003782533 0.001529349 0.006601008
#9        Warming      Rep1     non-sig 0.010941956 0.003719883 0.017239647
#10       Warming      Rep1         sig 0.003739143 0.001500314 0.005351209
#11       Warming      Rep2     non-sig 0.010146900 0.003908505 0.014028027
#12       Warming      Rep2         sig 0.003669121 0.001745241 0.004610831
#13       Warming      Rep3     non-sig 0.009570715 0.002961785 0.013823709
#14       Warming      Rep3         sig 0.003147567 0.001561752 0.005534449
#15       Warming      Rep4     non-sig 0.009429071 0.003298761 0.013321279
#16       Warming      Rep4         sig 0.003522856 0.001598337 0.004459520
#17           OWA      Rep1     non-sig 0.012542766 0.004289569 0.017202971
#18           OWA      Rep1         sig 0.005006164 0.001586814 0.011458945
#19           OWA      Rep2     non-sig 0.011485992 0.004572018 0.015536193
#20           OWA      Rep2         sig 0.005724609 0.001910180 0.011934641
#21           OWA      Rep3     non-sig 0.011206847 0.004937771 0.014080231
#22           OWA      Rep3         sig 0.005132834 0.002268531 0.010431511
#23           OWA      Rep4     non-sig 0.012428090 0.005479311 0.015465126
#24           OWA      Rep4         sig 0.007350657 0.002670261 0.020024583

mean_fst_df <- df_p %>%
  group_by(treatment, significant) %>%
  summarise(mean_fst = mean(fst, na.rm = TRUE),
                    median_fst = median(fst, na.rm = TRUE),
                    std_fst = sd(fst, na.rm=TRUE), .groups = 'drop')
mean_fst_df
#   treatment     significant mean_fst median_fst std_fst
#   <fct>         <fct>          <dbl>      <dbl>   <dbl>
#1 Acidification non-sig      0.00791    0.00303 0.0110
#2 Acidification sig          0.00388    0.00202 0.00576
#3 Warming       non-sig      0.0100     0.00350 0.0147
#4 Warming       sig          0.00352    0.00170 0.00500
#5 OWA           non-sig      0.0119     0.00488 0.0156
#6 OWA           sig          0.00580    0.00201 0.0140


# set number of sims
nsim <- 10000

# make empty variables to fill in the sims
ah_diff <- rep(NA, nsim)
ha_diff <- rep(NA, nsim)
hh_diff <- rep(NA, nsim)

ah_perm_sig <- rep(NA, nsim)
ha_perm_sig <- rep(NA, nsim)
hh_perm_sig <- rep(NA, nsim)

ah_perm_ns <- rep(NA, nsim)
ha_perm_ns <- rep(NA, nsim)
hh_perm_ns <- rep(NA, nsim)


for(i in 1:nsim){

        # set tmp df that I can just modify directly. 
        tmp_df <- dat_sub
    # randomize significance labels
    tmp_df$perm_ah_sig <- sample(tmp_df$ah_sig)
    tmp_df$perm_ha_sig <- sample(tmp_df$ha_sig)
    tmp_df$perm_hh_sig <- sample(tmp_df$hh_sig)

    # calculate mean for each group:
    ah_sig <- (tmp_df$perm_ah_sig == "sig")
    ha_sig <- (tmp_df$perm_ha_sig == "sig")
    hh_sig <- (tmp_df$perm_hh_sig == "sig")
    ah_perm_sig[i] <- mean(c(tmp_df$AH_F25_Rep1_fst[ah_sig],
                                             tmp_df$AH_F25_Rep2_fst[ah_sig],
                                             tmp_df$AH_F25_Rep3_fst[ah_sig],
                                             tmp_df$AH_F25_Rep4_fst[ah_sig]))
    ah_perm_ns[i] <- mean(c(tmp_df$AH_F25_Rep1_fst[!ah_sig],
                                             tmp_df$AH_F25_Rep2_fst[!ah_sig],
                                             tmp_df$AH_F25_Rep3_fst[!ah_sig],
                                             tmp_df$AH_F25_Rep4_fst[!ah_sig]))
    ha_perm_sig[i] <- mean(c(tmp_df$HA_F25_Rep1_fst[ha_sig],
                                             tmp_df$HA_F25_Rep2_fst[ha_sig],
                                             tmp_df$HA_F25_Rep3_fst[ha_sig],
                                             tmp_df$HA_F25_Rep4_fst[ha_sig]))
    ha_perm_ns[i] <- mean(c( tmp_df$HA_F25_Rep1_fst[!ha_sig],
                                             tmp_df$HA_F25_Rep2_fst[!ha_sig],
                                             tmp_df$HA_F25_Rep3_fst[!ha_sig],
                                             tmp_df$HA_F25_Rep4_fst[!ha_sig]))
    hh_perm_sig[i] <- mean(c(tmp_df$HH_F25_Rep1_fst[hh_sig],
                                             tmp_df$HH_F25_Rep2_fst[hh_sig],
                                             tmp_df$HH_F25_Rep3_fst[hh_sig],
                                             tmp_df$HH_F25_Rep4_fst[hh_sig]))
    hh_perm_ns[i] <-  mean(c(tmp_df$HH_F25_Rep1_fst[!hh_sig],
                                             tmp_df$HH_F25_Rep2_fst[!hh_sig],
                                             tmp_df$HH_F25_Rep3_fst[!hh_sig],
                                             tmp_df$HH_F25_Rep4_fst[!hh_sig]))

    # calculate differences to compate to our observed difference.
        # store this value

    ah_diff[i] <- ah_perm_ns[i] - ah_perm_sig[i]
    ha_diff[i] <- ha_perm_ns[i] - ha_perm_sig[i]
    hh_diff[i] <- hh_perm_ns[i] - hh_perm_sig[i]

    if(i %% 1000 == 0){print(i)}
}


    # calculate p value by comparing to the sampled distribution.
mean_fst_obs <- df_p %>%
  group_by(treatment, replicate, significant) %>%
  summarise(mean_fst = mean(fst, na.rm = TRUE),
                    median_fst = median(fst, na.rm = TRUE),
                    std_fst = sd(fst, na.rm=TRUE)) %>%
  group_by(treatment, significant) %>%
  summarise(mean_fst = mean(mean_fst, na.rm = TRUE),
                    median_fst = median(median_fst, na.rm = TRUE),
                    std_fst = sd(std_fst, na.rm=TRUE))
# # A tibble: 6 x 5
# Groups:   treatment [3]
#  treatment     significant mean_fst median_fst  std_fst
#  <fct>         <fct>          <dbl>      <dbl>    <dbl>
#1 Acidification non-sig      0.00791    0.00303 0.000553
#2 Acidification sig          0.00388    0.00213 0.00135
#3 Warming       non-sig      0.0100     0.00351 0.00178
#4 Warming       sig          0.00352    0.00158 0.000533
#5 OWA           non-sig      0.0119     0.00475 0.00128
#6 OWA           sig          0.00580    0.00209 0.00442


mean_fst_diff <- mean_fst_obs %>%
  select(treatment, significant, mean_fst) %>%
  spread(key = significant, value = mean_fst) %>%
  mutate(diff = `non-sig` - sig)
#  treatment     `non-sig`     sig    diff
# <fct>             <dbl>   <dbl>   <dbl>
#1 Acidification   0.00791 0.00388 0.00403
#2 Warming         0.0100  0.00352 0.00650
#3 OWA             0.0119  0.00580 0.00611

# find the observed difference in fst from the data
ah_diff_obs <- mean_fst_diff$diff[mean_fst_diff$treatment == "Acidification"]
ha_diff_obs <- mean_fst_diff$diff[mean_fst_diff$treatment == "Warming"]
hh_diff_obs <- mean_fst_diff$diff[mean_fst_diff$treatment == "OWA"]

# calc p-value:
pval_ah <- sum(c(ah_diff, ah_diff_obs) >= ah_diff_obs)/nsim  # two-tailed test
pval_ha <- sum(c(ha_diff, ha_diff_obs) >= ha_diff_obs)/nsim  # two-tailed test
pval_hh <- sum(c(hh_diff, hh_diff_obs) >= hh_diff_obs)/nsim  # two-tailed test
# all 1e-04

# make data frames to generate the plot
permPlot <- data.frame(
        treatment = c(rep("Acidification", length(ah_diff)), 
                                rep("Warming", length(ha_diff)),    
                                    rep("OWA", length(hh_diff))), 
        Fst_diff = c(ah_diff, ha_diff, hh_diff)
        )


obsPlot <- data.frame(
        treatment = c("Acidification","Warming","OWA"), 
        Fst_diff = c(ah_diff_obs, ha_diff_obs, hh_diff_obs)
        )

permPlot$treatment <- factor(permPlot$treatment, levels = c("Acidification", "Warming", "OWA"))

p2 <- ggplot(permPlot, aes(x=Fst_diff)) + 
    geom_histogram(color=NA, fill="grey65", bins=50) + 
    theme_classic()+
    facet_wrap(vars(treatment), nrow=1) +
    geom_segment(data=obsPlot, aes(x = Fst_diff, y = 0, xend = Fst_diff, yend = 900, colour = "segment"), 
            lwd=1.1, show.legend = F) +
    # ggtitle(paste("perm pval:", pval)) +
    xlab("Difference in Fst")+ 
    ylab("Frequency")+ 
    theme(plot.title = element_text(size=9),
                axis.text.x = element_text(size = 8)) 

ggsave( p2,
          file="~/tonsa_epigenetics/figures/window_reps_sig_FST_perm_1500kb.pdf",
          h=3, w=6.5)

#### 
# plot the actual distribution of the data. boxplots + points for each

p3 <- ggplot(df_p, aes(x =replicate, y = fst,
                      color=significant, fill=significant)) + 
    gghalves::geom_half_boxplot(side = "r", center = TRUE, errorbar.draw = FALSE,
                                fill="white", outliers=F) +
  gghalves::geom_half_point(
    side = "l", 
    range_scale = .3, 
    alpha = .3, size = 1.5
  ) +
  scale_colour_manual(values = c("non-sig" = "grey23", "sig"="firebrick3")) +
  theme_classic()+
  facet_wrap(vars(treatment), nrow=3)

ggsave( p3,
          file="~/tonsa_epigenetics/figures/window_reps_sig_FST_dotplot_1500kb.pdf",
          h=7, w=5)


#############################################
##
## 5 kb analysis
##
#############################################

# should probably write something so I don't have to repeat all. but good enough for now. 

dat_in <- read.table("~/tonsa_epigenetics/analysis/windows_fixed_5000kb_reps_fst_epi.txt", header=F) 
ncol(dat_in)
# 41 columns

colnames(dat_in) <- c("chr","start", "stop","SNP","aa_F00_mean_meth",
 "aa_F25_mean_meth","ah_F25_mean_meth","ha_F25_mean_meth","hh_F25_mean_meth",
    "ha_pval_meth","ha_fdr_meth","ah_pval_meth","ah_fdr_meth",
    "hh_pval_meth","hh_fdr_meth","ah_large","ha_large",
    "hh_large",
    "ah_delta_meth_F0","ha_delta_meth_F0","hh_delta_meth_F0",
    "ah_delta_meth_F25","ha_delta_meth_F25","hh_delta_meth_F25",
    # snp labels
    "chr_win", "start_win", "end_win", "window_id", "AH_F25_Rep1_fst", "AH_F25_Rep2_fst", "AH_F25_Rep3_fst", "AH_F25_Rep4_fst", "HA_F25_Rep1_fst", "HA_F25_Rep2_fst", "HA_F25_Rep3_fst", "HA_F25_Rep4_fst", "HH_F25_Rep1_fst", "HH_F25_Rep2_fst", "HH_F25_Rep3_fst", "HH_F25_Rep4_fst", "snp_count")

# add number of methylation sites per window:

summary_dat <- dat_in %>%
  group_by(window_id) %>%
  summarise(
    min_ha_fdr_meth = min(ha_fdr_meth, na.rm = TRUE),
    min_ah_fdr_meth = min(ah_fdr_meth, na.rm = TRUE),
    min_hh_fdr_meth = min(hh_fdr_meth, na.rm = TRUE),
    meth_count = n(), 
  )

# new df removing the methylation sites. keeping only row per snp window
snpdat <- dat_in[,25:ncol(dat_in)] %>% distinct(window_id, .keep_all = TRUE)
nrow(snpdat)
# 1733
# there's more... which makes sense bc larger windwos will intersect more methylation sites.

# join the two dfs:
alldat <- left_join(snpdat,summary_dat, by = "window_id")
nrow(alldat)
#1733

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
# 1102

sum(dat_sub$min_ah_fdr_meth< 0.05)
# 79
sum(dat_sub$min_ha_fdr_meth< 0.05)
# 106
sum(dat_sub$min_hh_fdr_meth< 0.05)
# 205

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
                                     rep(dat_sub$hh_sig, 4))
                            ) 

df_p$treatment <- factor(df_p$treatment, levels = c("Acidification", "Warming", "OWA"))

p1 <- ggplot(df_p, aes(x = significant, y = fst)) +
  theme_classic()  + 
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_summary(aes(group=treatment), fun = mean, geom = "line", 
              linewidth = 0.8) + 
  stat_summary(aes(group=treatment, shape=treatment, fill=treatment), fun = mean, geom = "point", 
              size = 5) +
  stat_summary(aes(group=replicate, shape=treatment, fill=treatment), fun = mean, geom = "point", 
              size = 2, alpha=0.5) +
  facet_wrap(vars(treatment), nrow=1) +
  scale_fill_manual(values=c("#F2AD00", "#00A08A", "#CC3333")) +
  scale_shape_manual(values=c(22,23,24)) +
  theme(legend.position="none") +
  xlab("Significant methylation change") +
  ylab("Mean Fst in 5 kb windows")

#ggsave( p,
#          file="~/tonsa_epigenetics/figures/window_reps_sig_FST_5000kb.pdf",
#          h=3, w=4.5)

# get summary of values:

mean_fst_df <- df_p %>%
  group_by(treatment, replicate, significant) %>%
  summarise(mean_fst = mean(fst, na.rm = TRUE),
                    median_fst = median(fst, na.rm = TRUE),
                    std_fst = sd(fst, na.rm=TRUE), .groups = 'drop')
as.data.frame(mean_fst_df)

#       treatment replicate significant    mean_fst  median_fst     std_fst
#1  Acidification      Rep1     non-sig 0.009696744 0.005574732 0.011097436
#2  Acidification      Rep1         sig 0.004229413 0.002300155 0.005408015
#3  Acidification      Rep2     non-sig 0.010278443 0.006243646 0.011746243
#4  Acidification      Rep2         sig 0.004804043 0.002971538 0.006687103
#5  Acidification      Rep3     non-sig 0.010157249 0.005492768 0.011819280
#6  Acidification      Rep3         sig 0.004322696 0.002341025 0.005671827
#7  Acidification      Rep4     non-sig 0.009561650 0.005149200 0.011110530
#8  Acidification      Rep4         sig 0.004255118 0.001697515 0.006664706
#9        Warming      Rep1     non-sig 0.012912676 0.007237421 0.014827305
#10       Warming      Rep1         sig 0.004515489 0.001861700 0.006245486
#11       Warming      Rep2     non-sig 0.012347177 0.006596547 0.014257383
#12       Warming      Rep2         sig 0.004655147 0.002103513 0.006503144
#13       Warming      Rep3     non-sig 0.011537793 0.005775474 0.013817688
#14       Warming      Rep3         sig 0.004168516 0.001903949 0.006511406
#15       Warming      Rep4     non-sig 0.011650611 0.006592637 0.014294673
#16       Warming      Rep4         sig 0.004288830 0.001740912 0.005519168
#17           OWA      Rep1     non-sig 0.016071786 0.010140762 0.017793252
#18           OWA      Rep1         sig 0.006334265 0.001700282 0.012891237
#19           OWA      Rep2     non-sig 0.013939338 0.008121303 0.015458763
#20           OWA      Rep2         sig 0.006814764 0.002226867 0.013042125
#21           OWA      Rep3     non-sig 0.013404692 0.008572654 0.014373027
#22           OWA      Rep3         sig 0.006021581 0.002634144 0.010966272
#23           OWA      Rep4     non-sig 0.015240966 0.010375516 0.017721143
#24           OWA      Rep4         sig 0.008040378 0.002782134 0.019471715

mean_fst_df <- df_p %>%
  group_by(treatment, significant) %>%
  summarise(mean_fst = mean(fst, na.rm = TRUE),
                    median_fst = median(fst, na.rm = TRUE),
                    std_fst = sd(fst, na.rm=TRUE), .groups = 'drop')
mean_fst_df
#  treatment     significant mean_fst median_fst std_fst
#  <fct>         <fct>          <dbl>      <dbl>   <dbl>
#1 Acidification non-sig      0.00992    0.00561 0.0114
#2 Acidification sig          0.00440    0.00236 0.00611
#3 Warming       non-sig      0.0121     0.00661 0.0143
#4 Warming       sig          0.00441    0.00198 0.00619
#5 OWA           non-sig      0.0147     0.00924 0.0164
#6 OWA           sig          0.00680    0.00238 0.0144


# set number of sims
nsim <- 10000

# make empty variables to fill in the sims
ah_diff <- rep(NA, nsim)
ha_diff <- rep(NA, nsim)
hh_diff <- rep(NA, nsim)

ah_perm_sig <- rep(NA, nsim)
ha_perm_sig <- rep(NA, nsim)
hh_perm_sig <- rep(NA, nsim)

ah_perm_ns <- rep(NA, nsim)
ha_perm_ns <- rep(NA, nsim)
hh_perm_ns <- rep(NA, nsim)


for(i in 1:nsim){

        # set tmp df that I can just modify directly. 
        tmp_df <- dat_sub
    # randomize significance labels
    tmp_df$perm_ah_sig <- sample(tmp_df$ah_sig)
    tmp_df$perm_ha_sig <- sample(tmp_df$ha_sig)
    tmp_df$perm_hh_sig <- sample(tmp_df$hh_sig)

    # calculate mean for each group:
    ah_sig <- (tmp_df$perm_ah_sig == "sig")
    ha_sig <- (tmp_df$perm_ha_sig == "sig")
    hh_sig <- (tmp_df$perm_hh_sig == "sig")
    ah_perm_sig[i] <- mean(c(tmp_df$AH_F25_Rep1_fst[ah_sig],
                                             tmp_df$AH_F25_Rep2_fst[ah_sig],
                                             tmp_df$AH_F25_Rep3_fst[ah_sig],
                                             tmp_df$AH_F25_Rep4_fst[ah_sig]))
    ah_perm_ns[i] <- mean(c(tmp_df$AH_F25_Rep1_fst[!ah_sig],
                                             tmp_df$AH_F25_Rep2_fst[!ah_sig],
                                             tmp_df$AH_F25_Rep3_fst[!ah_sig],
                                             tmp_df$AH_F25_Rep4_fst[!ah_sig]))
    ha_perm_sig[i] <- mean(c(tmp_df$HA_F25_Rep1_fst[ha_sig],
                                             tmp_df$HA_F25_Rep2_fst[ha_sig],
                                             tmp_df$HA_F25_Rep3_fst[ha_sig],
                                             tmp_df$HA_F25_Rep4_fst[ha_sig]))
    ha_perm_ns[i] <- mean(c( tmp_df$HA_F25_Rep1_fst[!ha_sig],
                                             tmp_df$HA_F25_Rep2_fst[!ha_sig],
                                             tmp_df$HA_F25_Rep3_fst[!ha_sig],
                                             tmp_df$HA_F25_Rep4_fst[!ha_sig]))
    hh_perm_sig[i] <- mean(c(tmp_df$HH_F25_Rep1_fst[hh_sig],
                                             tmp_df$HH_F25_Rep2_fst[hh_sig],
                                             tmp_df$HH_F25_Rep3_fst[hh_sig],
                                             tmp_df$HH_F25_Rep4_fst[hh_sig]))
    hh_perm_ns[i] <-  mean(c(tmp_df$HH_F25_Rep1_fst[!hh_sig],
                                             tmp_df$HH_F25_Rep2_fst[!hh_sig],
                                             tmp_df$HH_F25_Rep3_fst[!hh_sig],
                                             tmp_df$HH_F25_Rep4_fst[!hh_sig]))

    # calculate differences to compate to our observed difference.
        # store this value

    ah_diff[i] <- ah_perm_ns[i] - ah_perm_sig[i]
    ha_diff[i] <- ha_perm_ns[i] - ha_perm_sig[i]
    hh_diff[i] <- hh_perm_ns[i] - hh_perm_sig[i]

    if(i %% 1000 == 0){print(i)}
}


    # calculate p value by comparing to the sampled distribution.
mean_fst_obs <- df_p %>%
  group_by(treatment, replicate, significant) %>%
  summarise(mean_fst = mean(fst, na.rm = TRUE),
                    median_fst = median(fst, na.rm = TRUE),
                    std_fst = sd(fst, na.rm=TRUE)) %>%
  group_by(treatment, significant) %>%
  summarise(mean_fst = mean(mean_fst, na.rm = TRUE),
                    median_fst = median(median_fst, na.rm = TRUE),
                    std_fst = sd(std_fst, na.rm=TRUE))


mean_fst_diff <- mean_fst_obs %>%
  select(treatment, significant, mean_fst) %>%
  spread(key = significant, value = mean_fst) %>%
  mutate(diff = `non-sig` - sig)
#  treatment     `non-sig`     sig    diff
#  <fct>             <dbl>   <dbl>   <dbl>
#1 Acidification   0.00992 0.00440 0.00552
#2 Warming         0.0121  0.00441 0.00771
#3 OWA             0.0147  0.00680 0.00786

# find the observed difference in fst from the data
ah_diff_obs <- mean_fst_diff$diff[mean_fst_diff$treatment == "Acidification"]
ha_diff_obs <- mean_fst_diff$diff[mean_fst_diff$treatment == "Warming"]
hh_diff_obs <- mean_fst_diff$diff[mean_fst_diff$treatment == "OWA"]

  # calc p-value:
  pval_ah <- sum(c(ah_diff, ah_diff_obs) >= ah_diff_obs)/nsim  # two-tailed test
  pval_ha <- sum(c(ha_diff, ha_diff_obs) >= ha_diff_obs)/nsim  # two-tailed test
  pval_hh <- sum(c(hh_diff, hh_diff_obs) >= hh_diff_obs)/nsim  # two-tailed test
# all 1e-04

# make data frames to generate the plot

permPlot <- data.frame(
        treatment = c(rep("Acidification", length(ah_diff)), 
                                rep("Warming", length(ha_diff)),    
                                    rep("OWA", length(hh_diff))), 
        Fst_diff = c(ah_diff, ha_diff, hh_diff)
        )


obsPlot <- data.frame(
        treatment = c("Acidification","Warming","OWA"), 
        Fst_diff = c(ah_diff_obs, ha_diff_obs, hh_diff_obs)
        )

permPlot$treatment <- factor(permPlot$treatment, levels = c("Acidification", "Warming", "OWA"))

p3 <- ggplot(permPlot, aes(x=Fst_diff)) + 
    geom_histogram(color=NA, fill="grey65", bins=50) + 
    theme_classic()+
    facet_wrap(vars(treatment), nrow=1) +
    geom_segment(data=obsPlot, aes(x = Fst_diff, y = 0, xend = Fst_diff, yend = 900, colour = "segment"), 
            lwd=1.1, show.legend = F) +
    # ggtitle(paste("perm pval:", pval)) +
    xlab("Difference in Fst")+ 
    ylab("Frequency")+ 
    theme(plot.title = element_text(size=9),
                axis.text.x = element_text(size = 8)) 

#ggsave( p2,
#          file="~/tonsa_epigenetics/figures/window_reps_FST_combinedFig_5000kb.pdf",
#          h=3, w=6.5)

#### 
# plot the actual distribution of the data. boxplots + points for each

p2 <- ggplot(df_p, aes(x = treatment, y = fst),
        color=significant, fill=significant) +
   geom_boxplot(aes(color=significant),
    width = .2, fill = "white",
    size = 0.9
  ) +
  #gghalves::geom_half_point(
  #  side = "l", 
  #  range_scale = .3, 
  #  alpha = .3, size = 1.5
  #) +
  scale_colour_manual(values = c("non-sig" = "grey23", "sig"="firebrick3")) +
  theme_classic()+
    xlab("Significant methylation change") +
    ylab(paste0("mean Fst\nin ","5 kb windows"))+ 
    theme(legend.position="none")


#ggsave( p3,
##          file="~/tonsa_epigenetics/figures/window_reps_sig_FST_dotplot_5000kb.pdf",
#          h=7, w=5)
library(ggpubr)
#ggsave(ggarrange(ggarrange(p1, p2, ncol = 2, labels = c("A", "B")), p3, 
#          labels = c("","C"), nrow = 2),
#          file="~/tonsa_epigenetics/figures/window_FST_CombinedFig_5000kb.pdf",
#          h=6.5, w=6.5)

ggsave(ggarrange(p1, p3, ncol = 1, labels = c("A", "B")),
          file="~/tonsa_epigenetics/figures/window_FST_CombinedFig_5000kb.pdf",
          h=5.5, w=5.5)


```


# Flip order

## Significant AF change, whats the methylation change.

Keep replicates in there. 

For each window, calculate mean methylation, then determine if there's a significant af change or not via cmh. To see if this relationship is robust in both directions.


### get p-values from cmh

```bash

sed 's/\,/\t/g' ~/tonsa_genomics/analysis/results_table.csv >~/tonsa_epigenetics/analysis/cmh_results_table.txt

cat ~/tonsa_epigenetics/analysis/cmh_results_table.txt  | tail -n +2 | sort -k1,1 -k2,2n  >~/tonsa_epigenetics/analysis/cmh_results_table_sorted.txt

cut -f 8-11 ~/tonsa_epigenetics/analysis/cmh_results_table_sorted.txt  >~/tonsa_epigenetics/analysis/cmh_pvals_only.txt


cat ~/tonsa_epigenetics/analysis/cmh_results_table_sorted.txt | awk '{print $1 "\t" $2-1 "\t" $2}' | paste - ~/tonsa_epigenetics/analysis/cmh_pvals_only.txt | sort -k1,1 -k2,2n > ~/tonsa_epigenetics/analysis/snps.bed
# cols are chr, start, stop, aa_fdr,ah_fdr, ha_fdr, hh_fdr


```

Already have methylation table from above: `~/tonsa_epigenetics/analysis/gene_level_analysis/epi.bed`


### do window calculations


```bash

#######################
# fixed window size:
#######################

# make a loop to cycle through a couple of window sizes
for win_size in 1500 5000; do

    # make the windows
    bedtools makewindows -g <(cut -f 1,2 /data/copepods/tonsa_genome/GCA_900241095.1_Aton1.0_genomic.fa.fai) \
    -w $win_size \
    -s $win_size  \
    -i srcwinnum | sort -k1,1 -k2,2n  > ~/tonsa_epigenetics/analysis/windows_${win_size}kb.bed

# Join meth values and the 'windows.bed' file
    bedtools intersect \
    -a ~/tonsa_epigenetics/analysis/windows_${win_size}kb.bed \
    -b ~/tonsa_epigenetics/analysis/gene_level_analysis/epi.bed -wa -wb | sort -k4,4  \
    > ~/tonsa_epigenetics/analysis/windows_epi.tab  

# Run bedtools groupby command to obtain average values of methylation for each of the windows above
# columbs of epi.bed are:
    # ### they are from ~/tonsa_epigenetics/analysis/diff_methylation/methylation_summary.txt but with the 6 methylation differences added to the end
  # "chr","start", "stop" "window_id"
    # "chr","start", "stop","SNP","aa_F00_mean_meth",
  # "aa_F25_mean_meth","ah_F25_mean_meth","ha_F25_mean_meth","hh_F25_mean_meth",
    # "ha_pval_meth","ha_fdr_meth","ah_pval_meth","ah_fdr_meth",
    # "hh_pval_meth","hh_fdr_meth","ah_large","ha_large",
    # "hh_large",
    # col 23, 24, 25 "ah_delta_meth_F0","ha_delta_meth_F0","hh_delta_meth_F0",
    ## col 26,27,28 "ah_delta_meth_F25","ha_delta_meth_F25","hh_delta_meth_F25",

    # need to modify to grab the right columns
    bedtools groupby -i ~/tonsa_epigenetics/analysis/windows_epi.tab \
    -g 1,2,3,4 \
    -c  23,24,25 \
    -o mean > ~/tonsa_epigenetics/analysis/windows_mean_epi.txt

# count the number of methylation snps
    bedtools groupby -i ~/tonsa_epigenetics/analysis/windows_epi.tab \
    -g 1,2,3,4 \
    -c 8 \
    -o count  > ~/tonsa_epigenetics/analysis/windows_count_epi.txt

# join everything together:
    cut -f 5 ~/tonsa_epigenetics/analysis/windows_count_epi.txt | paste ~/tonsa_epigenetics/analysis/windows_mean_epi.txt - > ~/tonsa_epigenetics/analysis/epi_windows_${win_size}kb_epi.txt

######
# we now have epi in each window. Then  intersect this with the cmh file from above

# intersect fst with snps
    bedtools intersect \
    -a ~/tonsa_epigenetics/analysis/snps.bed \
    -b ~/tonsa_epigenetics/analysis/epi_windows_${win_size}kb_epi.txt -wa -wb   \
    > ~/tonsa_epigenetics/analysis/methylation_fixed_windows_${win_size}kb.txt

done

```

### plot results


```r

library(tidyverse)

dat_in <- read.table("~/tonsa_epigenetics/analysis/methylation_fixed_windows_1500kb.txt", header=F) 
ncol(dat_in)
# 15 columns

# cols are chr, start, stop, aa_fdr,ah_fdr, ha_fdr, hh_fdr
# # "chr","start", "stop" "window_id"
#  "ah_delta_meth_F0","ha_delta_meth_F0","hh_delta_meth_F0"
# meth_count


colnames(dat_in) <- c("chr","start", "stop",
                                            "aa_fdr","ah_fdr", "ha_fdr", "hh_fdr",
                                            "chr_win", "start_win", "end_win", "window_id",
                                            "ah_delta_meth_F0","ha_delta_meth_F0","hh_delta_meth_F0",
                                            "meth_count")


# add number of methylation sites per window:

summary_dat <- dat_in %>%
  group_by(window_id) %>%
  summarise(
    min_ah_fdr_snp = min(ah_fdr, na.rm = TRUE),
    min_ha_fdr_snp = min(ha_fdr, na.rm = TRUE),
    min_hh_fdr_snp = min(hh_fdr, na.rm = TRUE),
    snp_count = n(), 
  )

# new df removing the snp sites. keeping only row per snp window
methdat <- dat_in[,9:ncol(dat_in)] %>% distinct(window_id, .keep_all = TRUE)
nrow(methdat)
# 1440

# join the two dfs:
alldat <- left_join(methdat,summary_dat, by = "window_id")
nrow(alldat)
#1440

  new_dat <- alldat
  new_dat$ah_sig <- c("non-sig")
  new_dat$ha_sig <- c("non-sig")
  new_dat$hh_sig <- c("non-sig")
  new_dat$ah_sig[which(new_dat$min_ah_fdr< 0.05)] <- c("sig")
  new_dat$ha_sig[which(new_dat$min_ha_fdr< 0.05)] <- c("sig")
  new_dat$hh_sig[which(new_dat$min_hh_fdr< 0.05)] <- c("sig")

sum(new_dat$min_ah_fdr< 0.05)
# 296
sum(new_dat$min_ha_fdr< 0.05)
# 457
sum(new_dat$min_hh_fdr< 0.05)
# 549

dat_sub <- new_dat[which(new_dat$meth_count > 4),]
dat_sub <- new_dat[which(new_dat$snp_count > 4 & new_dat$meth_count > 4),]
nrow(dat_sub)
# 910

# make long format for plotting

repnum <- nrow(dat_sub)
df_p <- data.frame(
        treatment = c(rep("Acidification", repnum),rep("Warming", repnum),rep("OWA", repnum)),
    
        methylation_change = c(dat_sub$ah_delta_meth_F0,
                                                 dat_sub$ha_delta_meth_F0,
                                                 dat_sub$hh_delta_meth_F0),
        significant =c(dat_sub$ah_sig,
                                     dat_sub$ha_sig,
                                     dat_sub$hh_sig)
                            ) 

df_p$treatment <- factor(df_p$treatment, levels = c("Acidification", "Warming", "OWA"))

p1 <- ggplot(df_p, aes(x = significant, y = abs(methylation_change))) +
  theme_classic()  + 
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_summary(aes(group=treatment), fun = mean, geom = "line", 
              linewidth = 0.8) + 
  stat_summary(aes(group=treatment, shape=treatment, fill=treatment), fun = mean, geom = "point", 
              size = 5) +
  facet_wrap(vars(treatment), nrow=1) +
  scale_fill_manual(values=c("#F2AD00", "#00A08A", "#CC3333")) +
  scale_shape_manual(values=c(22,23,24)) +
  theme(legend.position="none") +
  xlab("Significant allele frequency change") +
    ylab(paste0("mean methylation change \nin ", "1.5 kb windows"))


p2 <- ggplot(df_p, aes(x = treatment, y = abs(methylation_change),
                color=significant, fill=significant)) +
     geom_boxplot(aes(color=significant),
    width = .2, fill = "white",
    size = 0.9
  ) +
  #gghalves::geom_half_point(
  #  side = "l", 
  #  range_scale = .3, 
  #  alpha = .3, size = 1.5
  #) +
  scale_colour_manual(values = c("non-sig" = "grey23", "sig"="firebrick3")) +
  theme_classic()+
    xlab("Significant allele frequency change") +
    ylab(paste0("mean methylation change \nin ","1.5 kb windows"))+ 
    theme(legend.position="none")


# get summary of values:

mean_meth_df <- df_p %>%
  group_by(treatment, significant) %>%
  summarise(mean_fst = mean(abs(methylation_change), na.rm = TRUE),
                    median_fst = median(abs(methylation_change), na.rm = TRUE),
                    std_fst = sd(abs(methylation_change), na.rm=TRUE), .groups = 'drop')
as.data.frame(mean_meth_df)

#       treatment significant   mean_fst  median_fst    std_fst
# 1 Acidification     non-sig 0.03173132 0.022956476 0.03043330
# 2 Acidification         sig 0.03084281 0.024814350 0.02912922
# 3       Warming     non-sig 0.03445347 0.023563903 0.03405581
# 4       Warming         sig 0.02262346 0.009145303 0.02819759
# 5           OWA     non-sig 0.04503009 0.038860494 0.03916480
# 6           OWA         sig 0.02870291 0.011097811 0.03625856


# set number of sims
nsim <- 10000

# make empty variables to fill in the sims
ah_diff <- rep(NA, nsim)
ha_diff <- rep(NA, nsim)
hh_diff <- rep(NA, nsim)

ah_perm_sig <- rep(NA, nsim)
ha_perm_sig <- rep(NA, nsim)
hh_perm_sig <- rep(NA, nsim)

ah_perm_ns <- rep(NA, nsim)
ha_perm_ns <- rep(NA, nsim)
hh_perm_ns <- rep(NA, nsim)


for(i in 1:nsim){

        # set tmp df that I can just modify directly. 
        tmp_df <- dat_sub
    # randomize significance labels
    tmp_df$perm_ah_sig <- sample(tmp_df$ah_sig)
    tmp_df$perm_ha_sig <- sample(tmp_df$ha_sig)
    tmp_df$perm_hh_sig <- sample(tmp_df$hh_sig)

    # calculate mean for each group:
    ah_sig <- (tmp_df$perm_ah_sig == "sig")
    ha_sig <- (tmp_df$perm_ha_sig == "sig")
    hh_sig <- (tmp_df$perm_hh_sig == "sig")
    ah_perm_sig[i] <- mean(tmp_df$ah_delta_meth_F0[ah_sig])
    ah_perm_ns[i] <- mean(tmp_df$ah_delta_meth_F0[!ah_sig])
    ha_perm_sig[i] <- mean(tmp_df$ha_delta_meth_F0[ha_sig])
    ha_perm_ns[i] <- mean(tmp_df$ha_delta_meth_F0[!ha_sig])
    hh_perm_sig[i] <- mean(tmp_df$hh_delta_meth_F0[hh_sig])
    hh_perm_ns[i] <- mean(tmp_df$hh_delta_meth_F0[!hh_sig])
    # calculate differences to compate to our observed difference.
        # store this value

    ah_diff[i] <- ah_perm_ns[i] - ah_perm_sig[i]
    ha_diff[i] <- ha_perm_ns[i] - ha_perm_sig[i]
    hh_diff[i] <- hh_perm_ns[i] - hh_perm_sig[i]

    if(i %% 1000 == 0){print(i)}
}


# calculate p value by comparing to the sampled distribution.
mean_meth_obs <- df_p %>%
  group_by(treatment, significant) %>%
  summarise(mean_meth = mean(methylation_change, na.rm = TRUE),
                    median_meth = median(methylation_change, na.rm = TRUE),
                    std_meth = sd(methylation_change, na.rm=TRUE))

#   treatment     significant mean_meth median_meth std_meth
#   <fct>         <fct>           <dbl>       <dbl>    <dbl>
# 1 Acidification non-sig        0.0296     0.0212    0.0325
# 2 Acidification sig            0.0288     0.0235    0.0311
# 3 Warming       non-sig        0.0309     0.0212    0.0373
# 4 Warming       sig            0.0181     0.00629   0.0313
# 5 OWA           non-sig        0.0425     0.0383    0.0419
# 6 OWA           sig            0.0249     0.00853   0.0390


mean_meth_diff <- mean_meth_obs %>%
  select(treatment, significant, mean_meth) %>%
  spread(key = significant, value = mean_meth) %>%
  mutate(diff = `non-sig` - sig)
## Groups:   treatment [3]
#  treatment     `non-sig`    sig     diff
#  <fct>             <dbl>  <dbl>    <dbl>
#1 Acidification    0.0296 0.0288 0.000765
#2 Warming          0.0309 0.0181 0.0129
#3 OWA              0.0425 0.0249 0.0176

# find the observed difference in fst from the data
ah_diff_obs <- mean_meth_diff$diff[mean_meth_diff$treatment == "Acidification"]
ha_diff_obs <- mean_meth_diff$diff[mean_meth_diff$treatment == "Warming"]
hh_diff_obs <- mean_meth_diff$diff[mean_meth_diff$treatment == "OWA"]

  # calc p-value:
  pval_ah <- sum(c(ah_diff, ah_diff_obs) >= ah_diff_obs)/nsim  # two-tailed test
  # 0.3876
  pval_ha <- sum(c(ha_diff, ha_diff_obs) >= ha_diff_obs)/nsim  # two-tailed test
  # 1e-04
  pval_hh <- sum(c(hh_diff, hh_diff_obs) >= hh_diff_obs)/nsim  # two-tailed test
    # 1e-04

# make data frames to generate the plot
permPlot <- data.frame(
        treatment = c(rep("Acidification", length(ah_diff)), 
                                rep("Warming", length(ha_diff)),    
                                    rep("OWA", length(hh_diff))), 
        Meth_diff = c(ah_diff, ha_diff, hh_diff)
        )


obsPlot <- data.frame(
        treatment = c("Acidification","Warming","OWA"), 
        Meth_diff = c(ah_diff_obs, ha_diff_obs, hh_diff_obs)
        )

permPlot$treatment <- factor(permPlot$treatment, levels = c("Acidification", "Warming", "OWA"))

p3 <- ggplot(permPlot, aes(x=Meth_diff)) + 
    geom_histogram(color=NA, fill="grey65", bins=50) + 
    theme_classic()+
    facet_wrap(vars(treatment), nrow=1) +
    geom_segment(data=obsPlot, aes(x = Meth_diff, y = 0, xend = Meth_diff, yend = 900, colour = "segment"), 
            lwd=1.1, show.legend = F) +
    # ggtitle(paste("perm pval:", pval)) +
    xlab("Difference in methylation")+ 
    ylab("Frequency")+ 
    theme(plot.title = element_text(size=9),
                axis.text.x = element_text(size = 8)) 

ggsave(ggarrange(ggarrange(p1, p2, ncol = 2, labels = c("A", "B")), p3, 
          labels = c("","C"), nrow = 2),
          file="~/tonsa_epigenetics/figures/window_Methylation_perm_1500kb.pdf",
          h=6.5, w=6.5)

```





