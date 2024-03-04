# edgeR analysis of RRBS

library(edgeR)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)

setwd("~/tonsa_epigenetics/analysis/diff_methylation/")

dir <- "/users/r/b/rbrennan/tonsa_epigenetics/analysis/methylation_extract/cpg_merged"
#list.files(dir)
# 
samples <- read.table("~/tonsa_epigenetics/analysis/methylation_extract/sample_id.txt", header=FALSE)
#Sample <- samples$V1[c(1,2,3,4,21,22,23,24)]
Sample <- samples$V1
# drop F3 samples, not informative and not comparable to other analyses.
Sample <- Sample[grep("HH_F03", Sample, invert=T)]
# now point to quant files
files <- file.path(dir, Sample)
names(files) <- gsub("_merged.gz.bismark.cov","",Sample)
all(file.exists(files))

nmlist <- (gsub("_merged.gz.bismark.cov","",Sample))

yall <- readBismark2DGE(files, sample.names=nmlist)
dim(yall)
# [1] 45022327       40
# look at counts
head(yall$counts)

# chromosomes
head(yall$genes)

# and supp info

yall$samples$group<- factor(substr(row.names(yall$samples), 1, 6))
yall$samples

########
# Filtering to remove low counts
########

#sum up the counts of methylated and unmethylatedd reads to get
    # the total read coverage at each CpG site for each sample

Methylation <- gl(2,1,ncol(yall), labels=c("Me","Un"))
Me <- yall$counts[, Methylation=="Me"]
Un <- yall$counts[, Methylation=="Un"]

Coverage <- Me + Un

head(Coverage)

#filter so sites have at least 15 x coverage
HasCoverage <- rowSums(Coverage >= 15) == nrow(yall$samples)/2

# We also filter out CpGs that are never methylated or always methylated as they provide no
    # information about differential methylation:
    # note, there are very few that pass HasCoverage and fail this
HasBoth <- rowSums(Me) > 0 & rowSums(Un) > 0

table(HasCoverage)
#   FALSE     TRUE
# 44923596    98731
table(HasBoth)
#   FALSE     TRUE
#35462655  9559672
table(HasCoverage, HasBoth)
#           HasBoth
#HasCoverage    FALSE     TRUE
#      FALSE 35462598  9460998
#      TRUE        57    98674

y1 <- yall[HasCoverage & HasBoth, keep.lib.sizes=FALSE]
# 98674 remain


## and then look at the excessively high coverage

Methylation <- gl(2,1,ncol(y1), labels=c("Me","Un"))
Me <- y1$counts[, Methylation=="Me"]
Un <- y1$counts[, Methylation=="Un"]

Coverage <- Me + Un

covSums <- (rowSums(Coverage)/(nrow(yall$samples)/2))

hist(covSums, xlim=c(0,2000), breaks=700, col="grey")

mean(covSums)
# 188
median(covSums)
# 79
quantile(covSums, c(0.975))
# 974.1175
HighCov <- (covSums < quantile(covSums, c(0.975)))
# cutoff is 974.1175
hist(covSums[HighCov], breaks=150, col="grey")

y <- y1[HighCov, keep.lib.sizes=FALSE]
# keep.lib.sizes=FALSE causes the library sizes to be recomputed after the filtering.
        # generally recommend this, although the effect on the downstream analysis is usually small.

# re calc mean cov:
Me1 <- y$counts[, grep("Me",colnames(y$counts)) ]
Un1 <- y$counts[, grep("Un",colnames(y$counts))]
Coverage1 <- Me1 + Un1
covSums2 <- (rowSums(Coverage1)/(nrow(y$samples)/2))

mean(covSums2)
#124.8
median(covSums2)
#77.7
min(covSums2)
# 20
sum(covSums2 < 50)
# 16778

16778/length(covSums2)
# 0.1743948

nrow(y$counts)
# [1] 96207

length(unique((y$genes$Chr)))
#[1] 7820

pdf("~/tonsa_epigenetics/figures/coverage.pdf", h=4, w=6)
hist(covSums2, breaks=50, col="grey75",
        main="Coverage of methylation sites",
        xlab="Coverage")
abline(v=mean(covSums2), col="blue", lwd=3)
abline(v=median(covSums2), col="firebrick3", lty=2, lwd=3)

dev.off()

# normalization
### To ensure that the methylated and unmethylated reads for the same sample are treated on the same scale,
    # we need to set the library sizes to be equal for each pair of libraries. We set the library sizes for each
    # sample to be the average of the total read counts for the methylated and unmethylated libraries:
TotalLibSize <- y$samples$lib.size[Methylation=="Me"] +
                y$samples$lib.size[Methylation=="Un"]

y$samples$lib.size <- rep(TotalLibSize, each=2)
y$samples


#####################################################################
#####################################################################
#Design matrix:
#####################################################################
#####################################################################

sampleIN <- read.table("~/tonsa_epigenetics/analysis/sampleID.txt")
colnames(sampleIN) <- c("combined")
sampleIN <- data.frame(combined = sampleIN$combined[grep("HH_F03", sampleIN$combined, invert=T)])

sampleID <- data.frame(
                condition = substr(sampleIN$combined, 1, 2),
                generation = substr(sampleIN$combined, 4, 6)
                )
row.names(sampleID) <- sampleIN$combined

Treatment <- paste(sampleID$condition, sampleID$generation, sep=".")
Treatment <- factor(Treatment)
#Replicate <- factor(sampleID$sample)
#designSL <- model.matrix(~Treatment)

designSL <- model.matrix(~0+Treatment)
design <- modelMatrixMeth(designSL)


#######################
#
# dispersion estimates
#
#######################

# The mean-dispersion relationship of BS-seq data has been studied in the past and no apparent mean-dispersion
    # trend was observed. This is also verified through our own practice.
    # therefore, we would not consider a meandependent dispersion trend as we normally would for RNA-seq data.
    # A common dispersion estimate for all the loci, as well as an empirical Bayes moderated dispersion for
    # each individual locus, can be obtained from the estimateDisp function in edgeR:

y <- estimateDisp(y, design, trend="none")
y$common.dispersion
# [1] 0.02649137

# so BVC is sqrt(dispersion)
sqrt(y$common.dispersion)
# or 0.1627617, which isn't bad.
y$prior.df
# [1] 29.0297

plotBCV(y)

# This returns a DGEList object with additional components (common.dispersion and tagwise.dispersion)
    # added to hold the estimated dispersions. 
        #Here the estimation of trended dispersion has been turned off by setting trend="none".

rm(yall)


#
# save R image to load later
#


save.image(file = "2022-02-21_meth_filtered.RData")



### LOAD R DATA!!
load(file = "~/tonsa_epigenetics/analysis/diff_methylation/2022-02-21_meth_filtered.RData")

library(ggplot2)
library(pheatmap)

##################
##################
#
# data exploration
#
##################
##################

#  M-value, which is defined as M = log2 {(Me + α)/(Un + α)}
    # where Me and Un are the methylated and unmethylated intensities and α is some suitable offset to
    # avoid taking logarithms of zero. The M-value can be interpreted as the base2 logit transformation
    # of the proportion of methylated signal at each locus.

Me <- y$counts[, Methylation=="Me"]
Un <- y$counts[, Methylation=="Un"]

M <- log2(Me + 2) - log2(Un + 2)
colnames(M) <- substr(colnames(Me), 1, 8)

mdsout <- plotMDS(M, top=nrow(M))

data <- mdsout$cmdscale.out

# make a good plot

data <- data.frame(id=row.names(mdsout$cmdscale.out), 
                Line=substr(row.names(mdsout$cmdscale.out), 1,2),
                    gen=substr(row.names(mdsout$cmdscale.out), 4,6),
                    group = substr(row.names(mdsout$cmdscale.out), 1,6),
                    Dimension1 =mdsout$cmdscale.out[,1],  
                    Dimension2= mdsout$cmdscale.out[,2])


d <- ggplot(data, aes(Dimension1, Dimension2, fill=group, shape=Line)) +
        geom_point(size=5) +
        xlab(paste0("Dimension1")) +
        ylab(paste0("Dimension2")) +
        theme_bw() +
        scale_shape_manual(values=c(21,22, 23, 24)) +
        scale_color_manual(values=c('black')) +
        scale_fill_manual(values=c('#D3DDDC','#6699CC', '#F2AD00', '#00A08A', '#CC3333'),
                labels = c( "Founding population", 
                            "Ambient", 
                            "Acidification",
                            "Warming",
                            "OWA"))+
        theme(legend.text=element_text(size=8))+
       # #theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))+
        guides(fill=guide_legend(override.aes=list(
                shape=c(21,21,22, 23, 24), 
                fill=c('#D3DDDC','#6699CC', '#F2AD00', '#00A08A', '#CC3333'))),
                shape= FALSE,
                   size=FALSE)+
               scale_size_manual(values=c(7,5))

d

ggsave("~/tonsa_epigenetics/figures/mds_edgeR.pdf",d, w=5.5, h=3.7)


########################
########################
###
### heatmaps, etc. from methylation percentages
###
########################
########################

# calculate methylation percentages

pm <- as.data.frame(matrix(ncol = (ncol(y$counts)/2), nrow = nrow(y$counts)))
colnames(pm) <- unique(substr(colnames(y$counts),1,8))
row.names(pm) <- row.names(y$counts)

for(i in 1:ncol(pm)){

    mat.tmp <- as.data.frame(y$counts[,grep(colnames(pm)[i], colnames(y$counts))])
    mat.tmp$total <- rowSums(mat.tmp)

    pm[,i] <- mat.tmp[,1]/mat.tmp$total
}


# plot methylation percentages

pm_mean <- colMeans(pm)
data <- data.frame(id=names(pm_mean), 
                    Line=substr(names(pm_mean), 1,2),
                    gen=substr(names(pm_mean), 4,6),
                    group = substr(names(pm_mean), 1,6),
                    mean_methylation = pm_mean)

d2 <- ggplot(data, aes(y=mean_methylation, x= group, fill=group, shape=group)) +
        geom_jitter(size=5, width=0.1) +
        xlab("") +
        ylab("Percent methylation") +
        theme_bw() +
       # ylim(-30, 23) + xlim(-50, 65)+
        scale_shape_manual(values=c(21,21,22, 23, 24)) +
        scale_color_manual(values=c('black')) +
        scale_fill_manual(values=c('#D3DDDC','#6699CC', '#F2AD00', '#00A08A', '#CC3333'),
                labels = c( "Founding population", 
                            "Ambient", 
                            "Acidification",
                            "Warming",
                            "OWA"))+
        theme(legend.text=element_text(size=8))+
       # #theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))+
        guides(fill=guide_legend(override.aes=list(
                shape=c(21,21,22, 23, 24), 
                fill=c('#D3DDDC','#6699CC', '#F2AD00', '#00A08A', '#CC3333'))),
                shape= FALSE,
                   size=FALSE)+
               scale_size_manual(values=c(7,5))

ggsave("~/tonsa_epigenetics/figures/meth_percent.pdf",d2, w=7, h=4.5)

write.table(file = "~/tonsa_epigenetics/analysis/methylation_mean_percent.txt", 
                                data,
                row.names=F, quote=F, sep="\t")

# read in data, if necessary

data <- read.table("~/tonsa_epigenetics/analysis/methylation_mean_percent.txt", header=T)

# run some quick stats on these:
m1 <- aov(data$mean_methylation ~ data$group)
anova(m1)
TukeyHSD(m1)

m1 <- aov(data$mean_methylation ~ data$group)
anova(m1)
TukeyHSD(m1)

mean(data$mean_methylation[data$gen=="F00"])
#0.254
sd(data$mean_methylation[data$gen=="F00"])
# 0.003
mean(data$mean_methylation[data$gen=="F25"])
# 229
sd(data$mean_methylation[data$gen=="F25"])
# 0.009

#######################
#
# Testing for differentially methylated CpG loci, strict method
#
#######################

#fit <- glmFit(y, design)
fit <- glmQLFit(y, design, robust=TRUE)
res <- glmQLFTest(fit)



#######################
#
# Testing for differentially methylated CpG loci, 
#### taking both F0 and F25 into acct

ave.hh <- makeContrasts(TreatmentHH.F25 -
                (TreatmentAA.F00  + TreatmentAA.F25 )/2, levels=design)

de.hh <- glmQLFTest(fit, contrast=ave.hh)
sig_de.hh <- topTags(de.hh, n=sum(summary(decideTests(de.hh))[c(1,3)]))[[1]]

ave.ha <- makeContrasts(TreatmentHA.F25 -
                (TreatmentAA.F00  + TreatmentAA.F25 )/2, levels=design)

de.ha <- glmQLFTest(fit, contrast=ave.ha)
sig_de.ha <- topTags(de.ha, n=sum(summary(decideTests(de.ha))[c(1,3)]))[[1]]

ave.ah <- makeContrasts(TreatmentAH.F25 -
                (TreatmentAA.F00  + TreatmentAA.F25 )/2, levels=design)

de.ah <- glmQLFTest(fit, contrast=ave.ah)
sig_de.ah <- topTags(de.ah, n=sum(summary(decideTests(de.ah))[c(1,3)]))[[1]]

nrow(sig_de.ah)
# [1] 633
nrow(sig_de.ha)
# [1] 531
nrow(sig_de.hh)
# [1] 1548

# parse down to the ones of interest
#hhsig <- pm[( row.names(pm) %in% hhF25_selection),]
hh_pm <- pm[( row.names(pm) %in% row.names(sig_de.hh)),]
ha_pm <- pm[( row.names(pm) %in% row.names(sig_de.ha)),]
ah_pm <- pm[( row.names(pm) %in% row.names(sig_de.ah)),]

# require minimum change

hh_aaf0 <- rowMeans(hh_pm[,1:4])
hh_aaf25 <- rowMeans(hh_pm[,5:8])
hh_hhf25 <- rowMeans(hh_pm[c("HH_F25_1","HH_F25_2","HH_F25_3","HH_F25_4")])

ha_aaf0 <-  rowMeans(ha_pm[,1:4])
ha_aaf25 <- rowMeans(ha_pm[,5:8])
ha_hhf25 <- rowMeans(ha_pm[c("HA_F25_1","HA_F25_2","HA_F25_3","HA_F25_4")])

ah_aaf0 <-  rowMeans(ah_pm[,1:4])
ah_aaf25 <- rowMeans(ah_pm[,5:8])
ah_hhf25 <- rowMeans(ah_pm[c("AH_F25_1","AH_F25_2","AH_F25_3","AH_F25_4")])

hhlarge_avg <- hh_pm[(abs(hh_aaf0 - hh_hhf25) > 0.1 & abs(hh_aaf25 - hh_hhf25) > 0.1 &
                        ((hh_aaf0 - hh_hhf25) > 0 & (hh_aaf25 - hh_hhf25) > 0 |
                         (hh_aaf0 - hh_hhf25) < 0 & (hh_aaf25 - hh_hhf25) < 0 )),]
halarge_avg <- ha_pm[(abs(ha_aaf0 - ha_hhf25) > 0.1 & abs(ha_aaf25 - ha_hhf25) > 0.1 &
                        ((ha_aaf0 - ha_hhf25) > 0 & (ha_aaf25 - ha_hhf25) > 0 |
                         (ha_aaf0 - ha_hhf25) < 0 & (ha_aaf25 - ha_hhf25) < 0 )),]
ahlarge_avg <- ah_pm[(abs(ah_aaf0 - ah_hhf25) > 0.1 & abs(ah_aaf25 - ah_hhf25) > 0.1 &
                        ((ah_aaf0 - ah_hhf25) > 0 & (ah_aaf25 - ah_hhf25) > 0 |
                         (ah_aaf0 - ah_hhf25) < 0 & (ah_aaf25 - ah_hhf25) < 0 )),]

nrow(hhlarge_avg)
# 753
nrow(halarge_avg)
# 161
nrow(ahlarge_avg)
# 128

sum(row.names(hhlarge_avg) %in% row.names(halarge_avg))
# 48
sum(row.names(hhlarge_avg) %in% row.names(ahlarge_avg))
# 50
sum(row.names(halarge_avg) %in% row.names(ahlarge_avg))
# 16

save.image(file = "~/tonsa_epigenetics/analysis/diff_methylation/2023-02-25_edgeR_complete.RData")






### LOAD R DATA
load(file = "~/tonsa_epigenetics/analysis/diff_methylation/2023-02-25_edgeR_complete.RData")

library(edgeR)


# write out table with stats:

#########################
# make a table of everything
#########################

ha_pval <- topTags(de.ha, n=nrow(y$gene))[[1]]$PValue
ha_fdr <- topTags(de.ha, n=nrow(y$gene))[[1]]$FDR
ha_snps <- row.names(topTags(de.ha, n=nrow(y$gene))[[1]])
ha_in <- data.frame(SNP = ha_snps, ha_pval_meth = ha_pval, ha_fdr_meth = ha_fdr)
ha_in$SNP <- gsub("-", ":", ha_in$SNP)

ah_pval <- topTags(de.ah, n=nrow(y$gene))[[1]]$PValue
ah_fdr <- topTags(de.ah, n=nrow(y$gene))[[1]]$FDR
ah_snps <- row.names(topTags(de.ah, n=nrow(y$gene))[[1]])
ah_in <- data.frame(SNP = ah_snps, ah_pval_meth = ah_pval, ah_fdr_meth = ah_fdr)
ah_in$SNP <- gsub("-", ":", ah_in$SNP)

hh_pval <- topTags(de.hh, n=nrow(y$gene))[[1]]$PValue
hh_fdr <- topTags(de.hh, n=nrow(y$gene))[[1]]$FDR
hh_snps <- row.names(topTags(de.hh, n=nrow(y$gene))[[1]])
hh_in <- data.frame(SNP = hh_snps, hh_pval_meth = hh_pval, hh_fdr_meth = hh_fdr)
hh_in$SNP <- gsub("-", ":", hh_in$SNP)

# require minimum change

pm_aaf00 <- rowMeans(pm[,1:4])
pm_aaf25 <- rowMeans(pm[,5:8])
pm_hhf25 <- rowMeans(pm[c("HH_F25_1","HH_F25_2","HH_F25_3","HH_F25_4")])
pm_haf25 <- rowMeans(pm[c("HA_F25_1","HA_F25_2","HA_F25_3","HA_F25_4")])
pm_ahf25 <- rowMeans(pm[c("AH_F25_1","AH_F25_2","AH_F25_3","AH_F25_4")])


outdat <- data.frame(CHR = y$gene$Chr, POS = y$gene$Locus, 
        SNP=paste(y$gene$Chr, y$gene$Locus, sep=":"),
    aa_F00_mean_meth = pm_aaf00, aa_F25_mean_meth = pm_aaf25,
    ah_F25_mean_meth = pm_ahf25, ha_F25_mean_meth = pm_haf25,
    hh_F25_mean_meth = pm_hhf25)

library(dplyr)
# the order of toptags is different!

outdat1 <- outdat %>% left_join(ha_in, by="SNP") %>% left_join(ah_in, by="SNP") %>% left_join(hh_in, by="SNP")

outdat1$ah_large <- outdat1$SNP %in% gsub("-", ":",row.names(ahlarge_avg))
outdat1$ha_large <- outdat1$SNP %in% gsub("-", ":",row.names(halarge_avg))
outdat1$hh_large <- outdat1$SNP %in% gsub("-", ":",row.names(hhlarge_avg))


write.table(file = "~/tonsa_epigenetics/analysis/diff_methylation/methylation_summary.txt", 
                                outdat1,
                row.names=F, quote=F, sep="\t")

outdat <- outdat1
# save pm


# write bed files for other analyses
pm1 <- cbind(outdat[,1:3], pm)
write.table(file = "~/tonsa_epigenetics/analysis/diff_methylation/methylation_percent.txt", pm1,
                row.names=F, quote=F, sep="\t")

## write bed files
sigdat <- outdat[(outdat$hh_large == TRUE),]
write.table(file = "~/tonsa_epigenetics/analysis/diff_methylation/methylation_hh.bed",
                cbind(as.character(sigdat$CHR), sigdat$POS -1, sigdat$POS),
                row.names=F, quote=F, sep="\t", col.names=F)

sigdat <- outdat[(outdat$ha_large == TRUE),]
write.table(file = "~/tonsa_epigenetics/analysis/diff_methylation/methylation_ha.bed",
                cbind(as.character(sigdat$CHR), sigdat$POS -1, sigdat$POS),
                row.names=F, quote=F, sep="\t", col.names=F)

sigdat <- outdat[(outdat$ah_large == TRUE),]
write.table(file = "~/tonsa_epigenetics/analysis/diff_methylation/methylation_ah.bed",
                cbind(as.character(sigdat$CHR), sigdat$POS -1, sigdat$POS),
                row.names=F, quote=F, sep="\t", col.names=F)


####################################################################
####################################################################
# UPset plot
####################################################################
####################################################################


# find unique significant
ah <- sum(outdat$ah_large == TRUE & outdat$ha_large == FALSE & outdat$hh_large == FALSE)
ha <- sum(outdat$ah_large == FALSE & outdat$ha_large == TRUE & outdat$hh_large == FALSE)
hh <- sum(outdat$ah_large == FALSE & outdat$ha_large == FALSE & outdat$hh_large == TRUE)

# pairwise overlaps

ah_ha <- sum(outdat$ah_large == TRUE & outdat$ha_large == TRUE & outdat$hh_large == FALSE)
ah_hh <- sum(outdat$ah_large == TRUE & outdat$ha_large == FALSE & outdat$hh_large == TRUE)
ha_hh <- sum(outdat$ah_large == FALSE & outdat$ha_large == TRUE & outdat$hh_large == TRUE)

# 3-way overlaps
ah_ha_hh <- sum(outdat$ah_large == TRUE & outdat$ha_large == TRUE & outdat$hh_large == TRUE)

all_input <- all_overlaps <- c(
                        "Acidic" = ah,
                        "Warm" = ha,
                        "OWA" = hh,
                "Acidic&Warm" = ah_ha,
                "Acidic&OWA" = ah_hh,
                "Warm&OWA" = ha_hh,
                "Acidic&Warm&OWA" = ah_ha_hh)




pdf("~/tonsa_epigenetics/figures/all_upSet.pdf",
    height = 4, width = 5)

upset(fromExpression(all_input),
        #order.by = "degree", 
        #group.by = "sets",
        keep.order = TRUE, empty.intersections = "on",
        sets = c("OWA","Warm", "Acidic"),
        mainbar.y.label = "Loci intersecting",
        sets.x.label = "Number of significant loci",
        point.size = 3.4, line.size = 1.2, 
          sets.bar.color=rev(c("#F2AD00", "#00A08A", "#CC3333")),
           text.scale = c(1.1, 1.1, 1, 1, 1.5, 1)
)

dev.off()

