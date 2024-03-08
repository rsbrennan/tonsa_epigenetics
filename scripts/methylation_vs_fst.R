
library(ggplot2)
library(ggbeeswarm)
library(ggpubr)

plot_windows <- function(dat_in, nm_in, snp_filter, fixed_windows=FALSE, win_size=NA, 
                          pop=NA, summary_stat=NA, meth_direction=NA,
                          fst_ylim_lower=NA, fst_ylim_upper=NA,
                          perm_xlim_lower=NA, perm_xlim_upper=NA){

  print(paste("running", pop, summary_stat))

  new_dat <- dat_in
  new_dat$sig <- c("non-sig")
  new_dat$sig[which(new_dat[paste(pop,"_fdr_af", sep="")]< 0.05)] <- c("sig")

  print(paste("Before snp density filtering: there are", nrow(new_dat), "rows in", nm_in, sep=" "))
  dat_sub <- new_dat[which(new_dat$n_snps > snp_filter),]
  print(paste("After snp density filtering: there are", nrow(dat_sub), "rows in", nm_in, sep=" "))

  # count up number of sig, nonsig

  filtered_ns <- dat_sub[which(dat_sub$sig == "non-sig"),]
  filtered_sig <- dat_sub[which(dat_sub$sig == "sig"),]

  print(paste("After snp density filtering: there are", length(unique(filtered_sig$win_id)), " windows with significant loci and", 
          length(unique(filtered_ns$win_id)), "windows with no significant loci", sep=" "))

    # deduplicate any rows. 
    chromes <- unique(dat_sub$win_id)
    if(length(chromes) == 0){
      stop("there no window id's. Are these windows based on epi positions?")
    }
    # 
    df_dedup <- as.data.frame(matrix(ncol=ncol(dat_sub), nrow=length(chromes)))
    colnames(df_dedup) <- colnames(dat_sub)
    df_dedup$CHR <- as.character(df_dedup$CHR)
    df_dedup$win_id <- as.character(df_dedup$win_id)

    # removing duplicate windows
    # basically, have the mean fst for each window, but lots of methylation values. 
    # so need to get just the one fst value for the window.
    for(i in 1:length(chromes)){
      tmp_df <- dat_sub[which(dat_sub$win_id == chromes[i]),]
      tmp_df$CHR <- as.character(tmp_df$CHR)
      tmp_df$win_id <- as.character(tmp_df$win_id)
      tmp_sig <- which(tmp_df$sig == "sig" )
      if(length(tmp_sig) > 0){
        df_dedup[i,] <- tmp_df[tmp_sig[1],]
      }
      else{
        df_dedup[i,] <- tmp_df[1,]
      }
    }
    # print stats here
    print(paste("After deduplicating: there are", length(unique(df_dedup$win_id)), "windows of which", 
      length(which(df_dedup$sig == "sig")), "have a significant loci and",
      length(which(df_dedup$sig == "non-sig")), "have no significant locus", sep=" "))
    # bc issues if pasting in ggplot y=. 
  df_dedup$fst_plot <- df_dedup[,paste0(pop,"_mean_fst")]
# alternate plot version:
# plot separately 
p1 <-ggplot(df_dedup, aes(x=sig, y=fst_plot,
                      group=sig,color=sig, fill=sig)) + 
  geom_boxplot(
    width = .2, fill = "white",
    size = 0.9, outlier.shape = NA
  ) +
  #ggdist::stat_halfeye(
    #adjust = .33, ## bandwidth
    #width = .67, 
  #  color = NA, ## remove slab interval
  #  position = position_nudge(x = .15)
  #) +
  gghalves::geom_half_point(
    side = "l", 
    range_scale = .3, 
    alpha = .3, size = 1.5
  ) +
  scale_colour_manual(values = c("non-sig" = "grey23", "sig"="firebrick3")) +
  theme_classic()+
  #annotate("text", x = 1, y = -0.007, label = round(mean(df_dedup[which(df_dedup$sig == "non-sig"),paste(pop,summary_stat,"fst", sep="_")]), 3)) +
  #  annotate("text", x = 2, y = -0.007, label = round(mean(df_dedup[which(df_dedup$sig == "sig"),paste(pop,summary_stat,"fst", sep="_")]), 3))+
    xlab("Significant methylation change") +
    ylab(paste0("mean Fst in ", win_size, " kb windows"))+ 
    theme(legend.position="none") +
    #ylim(0,0.12)+
    stat_summary(fun="mean", fill="#00FFFF",color="black", size=3, geom="point", pch=23)

# for manuscript plot: 
  # fst_ylim_lower = 0, fst_ylim_upper = 0.12

  #permutation to check if difference is signfificant

  set.seed(101) ## for reproducibility
  nsim <- 10000
  res <- numeric(nsim) ## set aside space for results

  for (i in 1:nsim) {
    ## standard approach: scramble response value
    tmp_df <- df_dedup
    tmp_df$perm_sig <- sample(tmp_df$sig)
    res[i] <- mean(tmp_df[which(tmp_df$perm_sig == "non-sig"),paste(pop, summary_stat,"fst", sep="_")]) -
          mean(tmp_df[which(tmp_df$perm_sig == "sig"),paste(pop,summary_stat,"fst", sep="_")])

    ## compute & store difference in means; store the value
  }

  obs <- mean(df_dedup[which(df_dedup$sig == "non-sig"),paste(pop,summary_stat,"fst", sep="_")]) -
          mean(df_dedup[which(df_dedup$sig == "sig"),paste(pop,summary_stat,"fst", sep="_")])
  ## append the observed value to the list of results

  # calc p-value:
  pval <- sum(abs(c(res, obs)) >= obs)/nsim  # two-tailed test

  print(paste("pvalue is ", pval))
  resdf <- as.data.frame(res)
  p2 <- ggplot(resdf, aes(res)) + 
    geom_histogram(color=NA, fill="grey65", bins=50) + 
    theme_classic()+
    geom_segment(aes(x = obs, y = 0, xend = obs, yend = 800, colour = "segment"), 
            lwd=1.1, show.legend = F) +
    # ggtitle(paste("perm pval:", pval)) +
    xlab("Difference in Fst")+ 
    ylab("Frequency")+ 
    theme(plot.title = element_text(size=9)) #+
    #xlim(-0.005,0.008)
# for manuscript plot: 
  # perm_xlim_lower = -0.005, perm_xlim_upper = 0.008

  return(ggarrange(p1,p2))

  #ggsave(file=paste0("~/tonsa_epigenetics/figures/sig_epi_fst_",j, "_", win_size,".pdf"), 
   #           p1, w=6.5, h=3.5 )
  #ggsave(file=paste0("~/tonsa_epigenetics/figures/sig_epi_fst_permutation_",j, "_", win_size,".pdf"), 
    #          p2, w=3.5, h=3 )

  }


##### read in files:
### fixed windows

fixed.1 <- read.table("~/tonsa_genomics/analysis/windows_fixed_500kb_fst_epi.txt", header=F)
fixed.2 <- read.table("~/tonsa_genomics/analysis/windows_fixed_1000kb_fst_epi.txt", header=F)
fixed.3 <- read.table("~/tonsa_genomics/analysis/windows_fixed_1500kb_fst_epi.txt", header=F)
fixed.4 <- read.table("~/tonsa_genomics/analysis/windows_fixed_5000kb_fst_epi.txt", header=F)

fixed_colnames <- c("CHR","epi_start", "epi_end",
                    "aa_meth_mean", "ah_meth_mean", "ha_meth_mean", "hh_meth_mean",
                    "ha_fdr_af", "ah_fdr_af", "hh_fdr_af", "win_chr",
                    "win_start", "win_end", "win_id", "ha_mean_fst","ah_mean_fst","hh_mean_fst",
                    "ha_max_fst","ah_max_fst","hh_max_fst", "n_snps")

colnames(fixed.1) <- fixed_colnames
colnames(fixed.2) <- fixed_colnames
colnames(fixed.3) <- fixed_colnames
colnames(fixed.4) <- fixed_colnames

snp_num <- 5
nm_in = "epi windows"
snp_filter=snp_num
fixed_windows=TRUE

# need to maually change these
win = 5 # window size, in kb, for name in file
dat_in <- fixed.4 # specify the correct input file that corresponds to the window size

for(i in c("mean")){

      a <- plot_windows(dat_in, nm_in = nm_in, snp_filter=snp_num,win_size=win,
                      pop="ah", summary_stat=i)
      b <- plot_windows(dat_in, nm_in = nm_in, snp_filter=snp_num, win_size=win,
                      pop="ha", summary_stat=i, meth_direction=m_dir)
      c <- plot_windows(dat_in, nm_in = nm_in, snp_filter=snp_num, win_size=win,
                      pop="hh", summary_stat=i, meth_direction=m_dir)
  
      ggsave(file= paste("~/tonsa_genomics/figures/fst_change",win,i,"fixed_windows.pdf", sep="_"),
            ggarrange(a,b,c, nrow=3, labels="AUTO"), h=7, w=6)
}



################################################################################
#
# check for mean methylation change in windows with sig Allele freq change
#
################################################################################


plot_windows <- function(dat_in, nm_in, snp_filter, fixed_windows=FALSE, win_size=NA, 
                          pop=NA, summary_stat=NA){

  print(paste("running", pop, summary_stat))

  new_dat <- dat_in
  new_dat$sig <- c("non-sig")
  new_dat$sig[which(new_dat[paste(pop,"_cmh", sep="")]< 0.05)] <- c("sig")

  print(paste("Before snp density filtering: there are", nrow(new_dat), "rows in", nm_in, sep=" "))
  dat_sub <- new_dat[which(new_dat$meth_count > snp_filter),]
  print(paste("After snp density filtering: there are", nrow(dat_sub), "rows in", nm_in, sep=" "))

  # count up number of sig, nonsig

  filtered_ns <- dat_sub[which(dat_sub$sig == "non-sig"),]
  filtered_sig <- dat_sub[which(dat_sub$sig == "sig"),]

  print(paste("After snp density filtering: there are", length(unique(filtered_sig$win_id)), " windows with significant loci and", 
          length(unique(filtered_ns$win_id)), "windows with no significant loci", sep=" "))

  # calc diff in methylation:
  dat_sub$meth_diff  <-   (dat_sub[,paste("aa_F00_meth_mean", sep="_")] - 
                              dat_sub[,paste(pop, "F25_meth", summary_stat, sep="_")])

    # deduplicate any rows. 
    chromes <- unique(dat_sub$win_id)
    if(length(chromes) == 0){
      stop("there no window id's. Are these windows based on epi positions?")
    }
    # 
    df_dedup <- as.data.frame(matrix(ncol=ncol(dat_sub), nrow=length(chromes)))
    colnames(df_dedup) <- colnames(dat_sub)
    df_dedup$CHR <- as.character(df_dedup$CHR)
    df_dedup$win_id <- as.character(df_dedup$win_id)

    # removing duplicate windows
    # basically, have the mean fst for each window, but lots of methylation values. 
    # so need to get just the one fst value for the window.
    for(i in 1:length(chromes)){
      tmp_df <- dat_sub[which(dat_sub$win_id == chromes[i]),]
      tmp_df$CHR <- as.character(tmp_df$CHR)
      tmp_df$win_id <- as.character(tmp_df$win_id)
      tmp_sig <- which(tmp_df$sig == "sig" )
      if(length(tmp_sig) > 0){
        df_dedup[i,] <- tmp_df[tmp_sig[1],]
      }
      else{
        df_dedup[i,] <- tmp_df[1,]
      }
    }
    # print stats here
    print(paste("After deduplicating: there are", length(unique(df_dedup$win_id)), "windows of which", 
      length(which(df_dedup$sig == "sig")), "have a significant loci and",
      length(which(df_dedup$sig == "non-sig")), "have no significant locus", sep=" "))

# alternate plot version:
# plot separately 
p1 <-ggplot(df_dedup, aes(x=sig, y=meth_diff,
                      group=sig,color=sig, fill=sig)) + 
  geom_boxplot(
    width = .2, fill = "white",
    size = 0.9, outlier.shape = NA
  ) +
  #ggdist::stat_halfeye(
    #adjust = .33, ## bandwidth
    #width = .67, 
  #  color = NA, ## remove slab interval
  #  position = position_nudge(x = .15)
  #) +
  gghalves::geom_half_point(
    side = "l", 
    range_scale = .3, 
    alpha = .3, size = 1.5
  ) +
  scale_colour_manual(values = c("non-sig" = "grey23", "sig"="firebrick3")) +
  theme_classic()+
  annotate("text", x = 1, y = -0.007, label = round(mean(df_dedup$meth_diff[which(df_dedup$sig == "non-sig")]), 3)) +
    annotate("text", x = 2, y = -0.007, label = round(mean(df_dedup$meth_diff[which(df_dedup$sig == "sig")]), 3))+
    xlab("Significant allele frequency change") +
    ylab(paste0("mean methylation change \nin ", win_size, " kb windows"))+ 
    theme(legend.position="none")

  #permutation to check if difference is signfificant

  set.seed(101) ## for reproducibility
  nsim <- 10000
  res <- numeric(nsim) ## set aside space for results

  for (i in 1:nsim) {
    ## standard approach: scramble response value
    tmp_df <- df_dedup
    tmp_df$perm_sig <- sample(tmp_df$sig)
    res[i] <- mean(tmp_df$meth_diff[which(tmp_df$perm_sig == "non-sig")]) -
          mean(tmp_df$meth_diff[which(tmp_df$perm_sig == "sig")])

    ## compute & store difference in means; store the value
  }

  obs <- mean(df_dedup$meth_diff[which(df_dedup$sig == "non-sig")]) -
          mean(df_dedup$meth_diff[which(df_dedup$sig == "sig")])
  ## append the observed value to the list of results

  # calc p-value:
  pval <- sum(abs(c(res, obs)) >= obs)/nsim  # two-tailed test

  print(paste("pvalue is ", pval))
  resdf <- as.data.frame(res)
  p2 <- ggplot(resdf, aes(res)) + 
    geom_histogram(color=NA, fill="grey65", bins=50) + 
    theme_classic()+
    geom_segment(aes(x = obs, y = 0, xend = obs, yend = 800, colour = "segment"), 
            lwd=1.1, show.legend = F) +
    # ggtitle(paste("perm pval:", pval)) +
    xlab("Difference in methylation")+ 
    ylab("Frequency")+ 
    theme(plot.title = element_text(size=9))

  return(ggarrange(p1,p2))

  }
