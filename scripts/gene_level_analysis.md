# gene level analysis

annotation gff: `~/tonsa_annotation/gawn/05_results/GCA_900241095.1_Aton1.0_genomic.fa.gff3`
methylation file: `~/tonsa_epigenetics/analysis/diff_methylation/methylation_summary.txt`
snp table: `~/tonsa_genomics/analysis/results_table.txt`


```bash
## make bed files

# remove genes and mRNA and cds from gff
cat ~/tonsa_annotation/gawn/05_results/GCA_900241095.1_Aton1.0_genomic.SORTED.gff3 |\
    awk '$3 != "gene"' |\
    awk '$3 != "mRNA"' |\
    awk '$3 != "CDS"' > ~/tonsa_annotation/gawn/05_results/GCA_900241095.1_Aton1.0_genomic.SORTED.exon.gff3

# use only genes gff
cat ~/tonsa_annotation/gawn/05_results/GCA_900241095.1_Aton1.0_genomic.SORTED.gff3 |\
    awk '$3 != "exon"' |\
    awk '$3 != "mRNA"' |\
    awk '$3 != "CDS"' > ~/tonsa_annotation/gawn/05_results/GCA_900241095.1_Aton1.0_genomic.SORTED.gene.gff3


##############################
##### repeat for epigenetics
##############################

## make bed files
#cat ~/tonsa_epigenetics/analysis/diff_methylation/methylation_summary.txt | tail -n +2 | awk '{print $1 "\t" $2-1 "\t" $2}' |sort -k1,1 -k2,2n > ~/tonsa_epigenetics/analysis/gene_level_analysis/snps.epi.bed

cat ~/tonsa_epigenetics/analysis/diff_methylation/methylation_summary.txt | tail -n +2 | awk '{print $1 "\t" $2-1 "\t" $2}' | sort -k1,1 -k2,2n > ~/tonsa_epigenetics/analysis/gene_level_analysis/meth.epi.bed

# first get closest of exon only

bedtools closest -D b -t all -b ~/tonsa_annotation/gawn/05_results/GCA_900241095.1_Aton1.0_genomic.SORTED.exon.gff3 \
                    -a ~/tonsa_epigenetics/analysis/gene_level_analysis/meth.epi.bed \
                    > ~/tonsa_epigenetics/analysis/gene_level_analysis/snp.exons.epi.bed
# this can return multiple matches for a single snp.

# Then closest for gene
bedtools closest -D b -t all -b ~/tonsa_annotation/gawn/05_results/GCA_900241095.1_Aton1.0_genomic.SORTED.gene.gff3 \
                    -a ~/tonsa_epigenetics/analysis/gene_level_analysis/meth.epi.bed \
                    > ~/tonsa_epigenetics/analysis/gene_level_analysis/snp.genes.epi.bed

```

now need to figure out which loci intersect a gene but not a exon. Build a big annotation table from this where I categorize each snp.

categories:
- exon
- intron
- upstream 
- downstream

Already done for SNPs, from 2022 paper.


```python
####################################################
####################################################
## assign genes for epigenetics
####################################################
####################################################

import pandas as pd
import numpy as np
import os

# files:
snp_file = os.path.expanduser("~/tonsa_epigenetics/analysis/diff_methylation/methylation_summary.txt")
annotation_file = os.path.expanduser("~/tonsa_annotation/gawn/annotation_merged_2020_05_27.txt")
gene_file = os.path.expanduser("~/tonsa_epigenetics/analysis/gene_level_analysis/snp.genes.epi.bed")
exon_file = os.path.expanduser("~/tonsa_epigenetics/analysis/gene_level_analysis/snp.exons.epi.bed")

############
# link each snp to its category
############

# read in snp file
snp_in = pd.read_csv(snp_file, delimiter="\t",header='infer')
# subset to only the first 3 cols
snp_ids = snp_in.iloc[:, range(3)]

# read in gene file
gene_in = pd.read_csv(gene_file, delimiter="\t",header=None)
# read in snp file
exon_in = pd.read_csv(exon_file, delimiter="\t",header=None)

gene_in['SNP'] = gene_in[0].astype(str) + ":" + gene_in[2].astype(str)
exon_in['SNP'] = exon_in[0].astype(str) + ":" + exon_in[2].astype(str)

# split out name:
gene_in['annotation'] = gene_in[11].str.split("Name=",n = 1, expand = True)[1]
exon_in['annotation'] = exon_in[11].str.split("Name=",n = 1, expand = True)[1]
gene_in.annotation.fillna(value="-", inplace=True)
exon_in.annotation.fillna(value="-", inplace=True)

gene_in[5] = gene_in[5].replace(".", "-")
exon_in[5] = exon_in[5].replace(".", "-")

snp_in["class"] = np.nan
snp_in["annotation"] = np.nan
snp_in["distance"] = np.nan

for idx, row in snp_in.iterrows():
    gene_match = gene_in[gene_in['SNP'] == (row['SNP'])].reset_index(drop=True)
    exon_match = exon_in[exon_in['SNP'] == (row['SNP'])].reset_index(drop=True)
    #### check for not matches.
    if (exon_match.empty and gene_match.empty):
        #print("no match!")
        snp_in.at[idx,'class'] = "-"
        snp_in.at[idx,'annotation'] = "-"
        snp_in.at[idx,'distance'] = "-"
    #### check if there's more than 1 match. multiple genes
    if len(gene_match.index) > 1 and len(exon_match.index) > 1:
    # if both match and not missing see which is closer. or if value is 0, take exon
        ### if both equal to zero, take exon. if exon in one, intron in another. call it a exon.
        if any(x == 0 for x in exon_match[12]):
            snp_in.at[idx,'class'] = "exon"
            # find which are 0
            idx_zero = exon_match.index[exon_match[12] == 0].tolist()
            for ex_row in idx_zero:
                snp_in.at[idx,'annotation'] = snp_in.at[idx,'annotation'] + ";" + exon_match['annotation'][ex_row].split(";")[0]
            snp_in.at[idx,'distance'] = 0
            # remove the nan
            snp_in.at[idx,'annotation'] = snp_in.at[idx,'annotation'].replace('nan;','')
        #### check for 0 in gene, and no 0 in exon these are introns. assing to both genes
        if all(x != 0 for x in exon_match[12]) and any(x == 0 for x in gene_match[12]):
            snp_in.at[idx,'class'] = "intron"
            # find the annotation that is 0
            idx_zero = gene_match.index[gene_match[12] == 0].tolist()
            for ex_row in idx_zero:
                snp_in.at[idx,'annotation'] = snp_in.at[idx,'annotation'] + ";" + gene_match['annotation'][ex_row].split(";")[0]
            snp_in.at[idx,'distance'] = 0
            # remove the nan
            snp_in.at[idx,'annotation'] = snp_in.at[idx,'annotation'].replace('nan;','')
        # check if no zeros, then find lowest.
        # In this case, it can't be an exon. must be either up or downstream
        if all(x  != 0 for x in exon_match[12]) and all(x != 0 for x in gene_match[12]):
            #print("none are zero")
            snp_in.at[idx,'annotation'] = gene_match['annotation'].iloc[gene_match[12].abs().idxmin()].split(";")[0]
            snp_in.at[idx,'distance'] = gene_match[12].iloc[gene_match[12].abs().idxmin()]
            if snp_in.at[idx,'distance'] < 0:
                snp_in.at[idx,'class'] = "promoter"
            if snp_in.at[idx,'distance'] > 0:
                snp_in.at[idx,'class'] = "downstream"
    #### check if one exon, but two genes:
    if len(gene_match.index) > 1 and len(exon_match.index) == 1:
        # if exon is 0, choose this one.
        if any(x == 0 for x in exon_match[12]):
            snp_in.at[idx,'class'] = "exon"
            idx_zero = exon_match.index[exon_match[12] == 0].tolist()
            snp_in.at[idx,'annotation'] = exon_match['annotation'].iloc[0].split(";")[0]
            snp_in.at[idx,'distance'] = 0
        #  check for 0 in gene, and no 0 in exon these are introns. assing to both genes
        if all(x != 0 for x in exon_match[12]) and any(x == 0 for x in gene_match[12]):
            #print("exon not 0, gene is 0")
            snp_in.at[idx,'class'] = "intron"
            # find the annotation that is 0
            idx_zero = gene_match.index[gene_match[12] == 0].tolist()
            for ex_row in idx_zero:
                snp_in.at[idx,'annotation'] = snp_in.at[idx,'annotation'] + ";" + gene_match['annotation'][ex_row].split(";")[0]
            snp_in.at[idx,'distance'] = 0
            # remove the nan
            snp_in.at[idx,'annotation'] = snp_in.at[idx,'annotation'].replace('nan;','')
        # check if no zeros, then find lowest.
        # In this case, it can't be an exon. must be either up or downstream
        if all(x  != 0 for x in exon_match[12]) and all(x != 0 for x in gene_match[12]):
            snp_in.at[idx,'annotation'] = gene_match['annotation'].iloc[gene_match[12].abs().idxmin()].split(";")[0]
            snp_in.at[idx,'distance'] = gene_match[12].iloc[gene_match[12].abs().idxmin()]
            if snp_in.at[idx,'distance'] < 0:
                snp_in.at[idx,'class'] = "promoter"
            if snp_in.at[idx,'distance'] > 0:
                snp_in.at[idx,'class'] = "downstream"
    # check if two exons, but one genes. only a few hundred of these:
    if len(gene_match.index) == 1 and len(exon_match.index) > 1:
        # if both equal to zero, take exon. if exon in one, intron in another. call it a exon.
        if any(x == 0 for x in exon_match[12]):
            snp_in.at[idx,'class'] = "exon"
            idx_zero = exon_match.index[exon_match[12] == 0].tolist()
            for ex_row in idx_zero:
                snp_in.at[idx,'annotation'] = snp_in.at[idx,'annotation'] + ";" + exon_match['annotation'][ex_row].split(";")[0]
            snp_in.at[idx,'annotation'] = ";".join(pd.unique(snp_in.at[idx,'annotation'].split(";")))
            snp_in.at[idx,'annotation'] = snp_in.at[idx,'annotation'].replace('nan;','')
            snp_in.at[idx,'distance'] = 0
        # first check for 0 in gene, and no 0 in exon these are introns. assing to both genes
        if all(x != 0 for x in exon_match[12]) and any(x == 0 for x in gene_match[12]):
            #print("exon not 0, gene is 0")
            snp_in.at[idx,'class'] = "intron"
            # find the annotation that is 0
            idx_zero = gene_match.index[gene_match[12] == 0].tolist()
            snp_in.at[idx,'annotation'] = gene_match['annotation'].iloc[0].split(";")[0]
            snp_in.at[idx,'distance'] = 0
        # check if no zeros, then find lowest.
        # In this case, it can't be an exon. must be either up or downstream
        if all(x  != 0 for x in exon_match[12]) and all(x != 0 for x in gene_match[12]):
            snp_in.at[idx,'annotation'] = gene_match['annotation'].iloc[0]
            snp_in.at[idx,'distance'] = gene_match[12].iloc[0]
            if snp_in.at[idx,'distance'] < 0:
                snp_in.at[idx,'class'] = "promoter"
            if snp_in.at[idx,'distance'] > 0:
                snp_in.at[idx,'class'] = "downstream"
    #if both match, do they have any actual annotation?
    if len(gene_match.index) == 1 and len(exon_match.index) == 1:
        if (exon_match['annotation'].iloc[0] == "-" and gene_match['annotation'].iloc[0] == "-"):
            snp_in.at[idx,'class'] = "-"
            snp_in.at[idx,'annotation'] = "-"
            snp_in.at[idx,'distance'] = "-"
    # if both match and not missing take exon if its 0.
    # otherwise, take intron, up, or down stream
        if (exon_match['annotation'].iloc[0] != "-" and gene_match['annotation'].iloc[0] != "-"):
            # if both equal to zero, take exon
            if (exon_match[12].iloc[0] == 0 and gene_match[12].iloc[0] == 0 ):
                snp_in.at[idx,'class'] = "exon"
                snp_in.at[idx,'annotation'] = exon_match['annotation'].iloc[0].split(";")[0]
                snp_in.at[idx,'distance'] = 0
            # if one is zero, take that one.
            else:
                # first check for 0 in one
                if (exon_match[12].iloc[0] == 0):
                    snp_in.at[idx,'class'] = "exon"
                    snp_in.at[idx,'annotation'] = exon_match['annotation'].iloc[0].split(";")[0]
                    snp_in.at[idx,'distance'] = 0
                if (exon_match[12].iloc[0] != 0 and gene_match[12].iloc[0] == 0 ):
                    snp_in.at[idx,'class'] = "intron"
                    snp_in.at[idx,'annotation'] = gene_match['annotation'].iloc[0].split(";")[0]
                    snp_in.at[idx,'distance'] = 0
                # check if upstream
                # if negative,
                if (gene_match[12].iloc[0] < 0 ):
                    snp_in.at[idx,'class'] = "promoter"
                    snp_in.at[idx,'annotation'] = gene_match['annotation'].iloc[0].split(";")[0]
                    snp_in.at[idx,'distance'] = gene_match[12].iloc[0]
                # check if downstream
                if (gene_match[12].iloc[0] > 0 ):
                    snp_in.at[idx,'class'] = "downstream"
                    snp_in.at[idx,'annotation'] = gene_match['annotation'].iloc[0].split(";")[0]
                    snp_in.at[idx,'distance'] = gene_match[12].iloc[0]
    if idx % 5000 == 0:
        print(idx)


snp_in.to_csv('/users/r/b/rbrennan/tonsa_genomics/analysis/gene_level_analysis/annotation_table_methylation.txt', index=False, sep='\t')


snp_in['class'].value_counts()

#-             43442
#exon          34594
#promoter       6558
#downstream     6251
#intron         5362


# snp_in.loc[snp_in['annotation'] == 'nan']

```

output file is: `/tonsa_genomics/analysis/gene_level_analysis/annotation_table_methylation.txt'`


## look at gene level results

```r

library(dplyr)

# read in snp table, from pnas paper
snp <- read.csv("~/tonsa_genomics/analysis/results_table.txt", header=T, sep="\t")
nrow(snp)
table(snp$class)

# read in methylation annotation table
meth <- read.csv("~/tonsa_genomics/analysis/gene_level_analysis/annotation_table_methylation.txt", header=T, sep="\t")
colnames(meth) <- gsub("_af", "",colnames(meth))
table(meth$class)
#          - downstream       exon     intron   promoter
#     43442       6251      34594       5362       6558

# read in percent methylation
pm <- read.csv("~/tonsa_epigenetics/analysis/diff_methylation/methylation_percent.txt", header=T, sep="\t")

meth_merged <- merge(pm, meth, by="SNP")
meth <- meth_merged

drops <- c("POS.y","CHR.y")
meth <- meth[ , !(names(meth) %in% drops)]
colnames(meth) <- c("SNP", "CHR", "POS", colnames(meth)[4:ncol(meth)])

# calculate change in methylation so these are carried over.
meth$aa_delta_F0_meth <- meth$aa_F00_mean_meth - meth$aa_F25_mean_meth
meth$ah_delta_F0_meth <- meth$aa_F00_mean_meth - meth$ah_F25_mean_meth
meth$ha_delta_F0_meth <- meth$aa_F00_mean_meth - meth$ha_F25_mean_meth
meth$hh_delta_F0_meth <- meth$aa_F00_mean_meth - meth$hh_F25_mean_meth

meth$ah_delta_F25_meth <- meth$aa_F25_mean_meth - meth$ah_F25_mean_meth
meth$ha_delta_F25_meth <- meth$aa_F25_mean_meth - meth$ha_F25_mean_meth
meth$hh_delta_F25_meth <- meth$aa_F25_mean_meth - meth$hh_F25_mean_meth


meth_relocated <- meth %>% relocate(class, .after = last_col()) %>% relocate(annotation, .after = last_col()) %>% relocate(distance, .after = last_col())

meth <- meth_relocated

write.table(meth , "~/tonsa_genomics/analysis/gene_level_analysis/methylation_results_table.txt", row.names=F, quote=F, sep="\t")

# drop some columns from meth, these are unnecessary
drops <- c("ah_large", "ha_large", "hh_large", "POS.y","CHR.y")
meth <- meth[ , !(names(meth) %in% drops)]
meth <- meth[ , !(names(meth) %in% drops)]

colnames(meth)[2:3] <- c("CHR", "POS")
colnames(meth)[4:23] <- paste(colnames(meth)[4:23], "_meth", sep="")
colnames(meth)[34:39] <- paste(colnames(meth)[34:39], "_meth", sep="")

######
# need to summarize data on a gene by gene basis to compare to gene expression.
# I think best here to do this separately for SNPs, methylation
# first do anything that hits a gene.
# then do genic (exon and intron)
# then exon only
# then promoter only

genes_meth <- (unique(unlist(strsplit(as.character(meth$annotation),split = ";"))))
genes_snp <- (unique(unlist(strsplit(as.character(snp$annotation),split = ";"))))
genes <- unique(c(genes_meth, genes_snp))

#####
# Methylation
#####

class_all <- unique(snp$class)
class_ids <- list(all = class_all, genic = c("exon", "intron"), exon = c("exon"), promoter = c("promoter"))

# Create empty list to hold dataframes
meth_out <- vector("list", length(class_ids))
names(meth_out) <- names(class_ids)

# Iterate over each class
for (class_name in names(class_ids)) {
  # Filter loci by class
  tmp_filt <- meth %>% filter(class %in% class_ids[[class_name]])
  
  # Get unique genes
  genes_filt <- unique(unlist(strsplit(as.character(tmp_filt$annotation), split = ";")))
  
  # Create the dataframe
    tmp_df <- as.data.frame(matrix(nrow=length(genes_filt), ncol=ncol(meth)+2))
    colnames(tmp_df) <- c("gene", "n_meth",
                        colnames(meth)[4:(ncol(meth)-3)],
                        "ah_fdr_meth_min","ha_fdr_meth_min","hh_fdr_meth_min",
                        "ah_delta_meth_min","ha_delta_meth_min","hh_delta_meth_min"
                        )
    tmp_df$gene <- genes_filt

  # Iterate over each gene
  for (i in seq_along(genes_filt)) {
    tmp_gene <- tmp_filt[grep(genes_filt[i], tmp_filt$annotation), ]
    
    # Populate the dataframe
    tmp_df$gene[i] <- genes_filt[i]
    tmp_df$n_meth[i] <- nrow(tmp_gene)
    tmp_df[i, 3:(ncol(tmp_df) - 6)] <- colMeans(tmp_gene[, 4:(ncol(tmp_gene) - 3)], na.rm = TRUE)
    tmp_df$ah_fdr_meth_min[i] <- min(tmp_gene$ah_fdr_meth, na.rm = TRUE)
    tmp_df$ha_fdr_meth_min[i] <- min(tmp_gene$ha_fdr_meth, na.rm = TRUE)
    tmp_df$hh_fdr_meth_min[i] <- min(tmp_gene$hh_fdr_meth, na.rm = TRUE)
    # Compute the index of the maximum absolute value
    max_abs_ah <- which.max(abs(tmp_gene$ah_delta_F25_meth))
    max_abs_ha <- which.max(abs(tmp_gene$ha_delta_F25_meth))
    max_abs_hh <- which.max(abs(tmp_gene$hh_delta_F25_meth))
    tmp_df$ah_delta_meth_min[i] <- tmp_gene$ah_delta_F25_meth[max_abs_ah]
    tmp_df$ha_delta_meth_min[i] <- tmp_gene$ha_delta_F25_meth[max_abs_ha]
    tmp_df$hh_delta_meth_min[i] <- tmp_gene$hh_delta_F25_meth[max_abs_hh]
  
    if (i %% 500 == 0) {
      print(i)
    }
  }
  
  # Add to the list of dataframes
  meth_out[[class_name]] <- tmp_df
  print(paste("done with", class_name))
}


#############
# repeat for snps:
#############

drops <- c("aa_sig", "ah_sig", "ha_sig", "hh_sig", "ah_sig_nolab", "ha_sig_nolab", "hh_sig_nolab", "lab_sln_sig")
snp <- snp[ , !(names(snp) %in% drops)]

# make empty list to hold dataframes that we'll fill in
snp_out <- vector(mode = "list", length = 4)
names(snp_out) <- c("all", "genic", "exon", "promoter")

for(class_index in 1:length(class_ids)){
    # filter to get only the loci in the class of interest.
    tmp_filt <- snp %>% filter(class %in% as.vector(class_ids[[class_index]]) )
    # get the unique genes from the filtered dataset we just made
    genes_filt <- (unique(unlist(strsplit(as.character(tmp_filt$annotation),split = ";"))))
    # make the dataframe that we need to fill in
    tmp_df <- as.data.frame(matrix(nrow=length(genes_filt), ncol=ncol(snp)-1))
    colnames(tmp_df) <- c("gene", "n_snp",
                        colnames(snp)[4:(ncol(snp)-4)],
                        "aa_pval_AF_min","ah_pval_AF_min","ha_pval_AF_min","hh_pval_AF_min"
                        )
    # loop over all the genes in this df, save mean, number snps, etc.
    for(i in 1:length(genes_filt)){
        # pull down loci that match each gene
        tmp_gene <- tmp_filt[grep(genes_filt[i], tmp_filt$annotation),]
        # add gene to output df
        tmp_df$gene[i] <- genes_filt[i]

        # save number of loci
        tmp_df$n_snp[i] <- nrow(tmp_gene)

        # take colmeans 
        tmp_df[i, 3:(ncol(tmp_df)-4)] <- colMeans(tmp_gene[,4:(ncol(tmp_gene)-4)], na.rm=TRUE)

        # let's also take the min p value for each set
        tmp_df$aa_pval_AF_min[i] <- min(tmp_gene$aa_fdr, na.rm=TRUE)
        tmp_df$ah_pval_AF_min[i] <- min(tmp_gene$ah_fdr, na.rm=TRUE)
        tmp_df$ha_pval_AF_min[i] <- min(tmp_gene$ha_fdr, na.rm=TRUE)
        tmp_df$hh_pval_AF_min[i] <- min(tmp_gene$hh_fdr, na.rm=TRUE)

        if(i %% 500 == 0){print(i)}
    }

    # add to the list of dataframes
    snp_out[[names(class_ids)[class_index]]] <- tmp_df
    print("done with")
    print(names(class_ids)[class_index])
}

save.image(file = "~/tonsa_genomics/analysis/gene_level_analysis/2023-02_24_gene_level.RData")


#################################################################################


load(file = "~/tonsa_genomics/analysis/gene_level_analysis/2023-02_24_gene_level.RData")

library(ggplot2)
library(dplyr)
library(ggExtra)
library(ggpubr)
library(ggbeeswarm)
library(plyr)
library(ggridges)

# make merged list to compare

snp_meth_list <- list()
snp_meth_list[[1]] <- merge(snp_out[['all']], meth_out[['all']], by="gene")

names(snp_meth_list) <- "all"
#names(snp_meth_list) <- names(snp_out)


##################################################################################
####################
# plot by pvalues of methylation vs af changes
####################
##################################################################################

# below will also plot by category for each gene.

# make df to hold Rsq, slope, p-value

##### Figure S2


lm_df <- as.data.frame(matrix(ncol=6, nrow=0))
colnames(lm_df) <- c("treatment","gene_set", "measure", "n_snps", "rsq", "pval")


for(i in 4){ # this sets the number of snps. greater than or equal to the values here

    for(j in 1:1){ # cycle over the snp classes.

    fout <- snp_meth_list[[j]] %>% filter(n_snp >= i & n_meth >= i )
    fout[,2:ncol(fout)] <- sapply(fout[,(2:ncol(fout))],as.numeric)

    # convert all p-values to -log10(p)
    fout[,grep("fdr",colnames(fout))] <-    lapply(fout[,grep("fdr",colnames(fout))], function(x)  {
                                                            x <- -log10(x)
                                                            x
                                                        }
                                                    )

    # add to df to save:

    # ah
        tmp_lm <-summary(lm(fout$ah_fdr_meth ~ fout$ah_fdr))
        lm_trt <- c("AH")
        lm_gen <- c("F0-F25")
        lm_df[nrow(lm_df)+1,] <- c(lm_trt,names(snp_meth_list)[j],lm_gen, i, round(tmp_lm$r.squared, 4), round(tmp_lm$coefficients[2,4], 3))

        p2 <- ggplot(fout, aes(x=ah_fdr, y=ah_fdr_meth)) +
                geom_point(pch=21, alpha=0.3, fill="black", size=2) +   theme_classic() +
                geom_smooth(method = lm, se = TRUE) + 
#               coord_cartesian(xlim=c(0, 3), ylim= c(0,2.5) ) +
                ggtitle(paste("Acidification: rsq:",
                                round(tmp_lm$r.squared, 3), ";",
                                "pval:", round(tmp_lm$coefficients[2,4], 3) ))+
                theme(plot.title = element_text(size=10))



        tmp_lm <-summary(lm(fout$ha_fdr_meth ~ fout$ha_fdr))
        lm_trt <- c("HA")
        lm_gen <- c("F0-F25")
        lm_df[nrow(lm_df)+1,] <- c(lm_trt,names(snp_meth_list)[j],lm_gen, i, round(tmp_lm$r.squared, 4), round(tmp_lm$coefficients[2,4], 3))

        p3 <- ggplot(fout, aes(x=ha_fdr, y=(ha_fdr_meth))) +
                geom_point(pch=21, alpha=0.3, fill="black", size=2) +   theme_classic() +
                geom_smooth(method = lm, se = TRUE) + 
#               coord_cartesian(xlim=c(0, 3) ) +
                ggtitle(paste("Warming: rsq:",
                                round(tmp_lm$r.squared, 3), ";",
                                "pval:", round(tmp_lm$coefficients[2,4], 3) ))+
                theme(plot.title = element_text(size=10))

        tmp_lm <-summary(lm(fout$hh_fdr_meth ~ fout$hh_fdr))
        lm_trt <- c("HH")
        lm_gen <- c("F0-F25")
        lm_df[nrow(lm_df)+1,] <- c(lm_trt,names(snp_meth_list)[j],lm_gen, i, round(tmp_lm$r.squared, 4), round(tmp_lm$coefficients[2,4], 3))
        p4 <- ggplot(fout, aes(x=(hh_fdr), y=(hh_fdr_meth))) + 
                geom_point(pch=21, alpha=0.3, fill="black", size=2) +   theme_classic() +
                geom_smooth(method = lm, se = TRUE) + 
    #           coord_cartesian(xlim=c(0, 3)) +
                ggtitle(paste("Greenhouse: rsq:",
                                round(tmp_lm$r.squared, 3), ";",
                                "pval:", round(tmp_lm$coefficients[2,4], 3) ))+
                theme(plot.title = element_text(size=10))

        ggsave(paste("~/tonsa_epigenetics/figures/meth_vs_cmh_F0_", names(snp_meth_list)[j], "_",  i, "_loci.pdf", sep=""), ggarrange(p2, p3, p4, nrow=1), w=10, h=3)
}
}

########################################################################
########################
# plot by methylation vs DGE
########################
########################################################################

# read in dge:
dge <- read.csv("~/tonsa_epigenetics/analysis/DGE_HHHH_vs_AAAAA.txt", header=T, sep="\t")
colnames(dge) <- c("baseMean", colnames(dge)[2:ncol(dge)] )
dge$gene <- row.names(dge)


dge_meth_list <- list()
dge_meth_list[[1]] <- inner_join(dge, meth_out[['all']], by="gene")

names(dge_meth_list) <- "all"
sapply(dge_meth_list, nrow)
#     all    genic     exon promoter
#    2824     2340     1762      413


# can categorize each gene as having a significant snp/methylation or not. this is arbitrary for snps, can just use quantile or something?

dge_meth_list[[1]]$category <- "all"

dge_meth_all <- bind_rows(dge_meth_list)

### DGE sig or not, Meth sig or not
dge_meth <- dge_meth_list[[1]]
# drop na's
dge_meth <- dge_meth[!is.na(dge_meth$padj),]

# assign sig or not:
dge_meth$Meth_hh_sig <- ifelse(dge_meth$hh_fdr_meth_min < 0.05, TRUE, FALSE)
dge_meth$DGE_sig <- ifelse(dge_meth$padj < 0.05, TRUE, FALSE)

dge_meth <- dge_meth[which(dge_meth$n_meth > 4),]
nrow(dge_meth)
#1772

chisq.test(table(dge_meth$DGE_sig, dge_meth$Meth_hh_sig))$expected
chisq.test(table(dge_meth$DGE_sig, dge_meth$Meth_hh_sig),correct = TRUE)
# p-value = 0.02876, but lower than expected. 

################################
# make scatter plot
################################

dge_meth_sig <- dge_meth[which(dge_meth$Meth_hh_sig == TRUE),]
nrow(dge_meth_sig)
# 197

head(dge_meth_sig[order(abs(dge_meth_sig$hh_delta_F25_meth)), ])


fit1 <- lm(log2FoldChange ~ hh_delta_F25_meth, data = dge_meth_sig)
summary(fit1)
#Multiple R-squared:  0.03003,  Adjusted R-squared:  0.02505
# F-statistic: 6.036 on 1 and 195 DF,  p-value: 0.01489

# change in methylation is 
      #  meth$hh_delta_F0_meth <- meth$aa_F00_mean_meth - meth$hh_F25_mean_meth
      # so positive is decrease in HH, negative is increase in HH.

p1_scatter <- ggplot(dge_meth_sig, aes(y=log2FoldChange, x = (hh_delta_F25_meth*-1))) +
        geom_hline(yintercept=0, linetype="dashed", color="grey60") +
        geom_vline(xintercept=0, linetype="dashed", color="grey60") +
        geom_point(alpha=0.8, shape=21, fill="black", size=1.5)+
        stat_smooth(method = "lm") +
        # geom_beeswarm() +
        theme_classic() +
      xlab("Mean change in OWA methylation") +
      labs(y = expression(log[2]~(OWA/Ambient)))  # Y-axis label with subscript

ggsave(p1_scatter, file= ("~/tonsa_epigenetics/figures/epi_dge_scatter_sig.pdf"), h=3, w=4)



#########----------------------------------------------------------------------------
# max rather than mean


dge_meth_sig2 <- dge_meth[which(dge_meth$Meth_hh_sig == TRUE & abs(dge_meth$hh_delta_meth_min) > 0.1),]
nrow(dge_meth_sig2)
#154

fit1 <- lm(log2FoldChange ~ hh_delta_meth_min, data = dge_meth_sig2)
summary(fit1)
# Multiple R-squared:  0.02743, Adjusted R-squared:  0.02103
# F-statistic: 4.287 on 1 and 152 DF,  p-value: 0.0401

p2_scatter <- ggplot(dge_meth_sig2, aes(y=log2FoldChange, x = (hh_delta_meth_min*-1))) +
        geom_hline(yintercept=0, linetype="dashed", color="grey60") +
        geom_vline(xintercept=0, linetype="dashed", color="grey60") +
        geom_point(alpha=0.8, shape=21, fill="black", size=1.5)+
        stat_smooth(method = "lm") +
        # geom_beeswarm() +
        theme_classic() +
      xlab("Max change in OWA methylation") +
      labs(y = expression(log[2]~(OWA/Ambient))) # Add inherit.aes = FALSE

ggsave(p2_scatter, file= ("~/tonsa_epigenetics/figures/epi_dge_scatter_sig_MAX.pdf"), h=3, w=4)




########################################################################################################
# GE plasticity specific genes, any relationship to methylation?
########################################################################################################

# ~/tonsa_epigenetics/analysis/DGE_HHHH_vs_AAHH.txt
#~/tonsa_epigenetics/analysis/DGE_AAAA_vs_HHAA.txt
dge <- read.csv("~/tonsa_epigenetics/analysis/DGE_HHHH_vs_AAHH.txt", header=T, sep="\t")
colnames(dge) <- c("baseMean", colnames(dge)[2:ncol(dge)] )
dge$gene <- row.names(dge)

dge_meth_list <- list()
dge_meth_list[[1]] <- inner_join(dge, meth_out[['all']],      by="gene")

# calc some summary info:
dge_meth <- dge_meth_list[[1]]
dge_meth <- dge_meth[!is.na(dge_meth$padj),]

# assign sig or not:
dge_meth$DGE_sig <- ifelse(dge_meth$padj < 0.05, TRUE, FALSE)

dge_meth <- dge_meth[which(dge_meth$n_meth > 4),]
nrow(dge_meth)
#1721
mu <- ddply(dge_meth, "DGE_sig", summarise, grp.mean=mean((hh_F25_mean_meth)))
#  DGE_sig   grp.mean
#1   FALSE 0.11467556
#2    TRUE 0.02728438
ster <- function(x) sd(x)/sqrt(length(x))
mu2 <- ddply(dge_meth, "DGE_sig", summarise, grp.se=ster((hh_F25_mean_meth)))
#  DGE_sig      grp.se
#1   FALSE 0.004657732
#2    TRUE 0.004277920
ks.test((dge_meth$hh_F25_mean_meth[dge_meth$DGE_sig == TRUE]), 
        (dge_meth$hh_F25_mean_meth[dge_meth$DGE_sig == FALSE]))
# D = 0.25659, p-value = 0.01929

p1 <- ggplot(dge_meth, aes(x=hh_F25_mean_meth,y=DGE_sig,fill=DGE_sig, group= DGE_sig)) +
  geom_density_ridges(scale = 4.5, alpha=0.3,
    jittered_points = TRUE,
    position = position_points_jitter(width = 0.0, height = 0),
    point_shape = '', point_size = 1.5, point_alpha = 0.5
  ) + 
  scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
  scale_x_continuous(expand = c(0, 0)) + 
  coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
  theme_ridges() +
        scale_fill_manual(values=c("gray50", "#9B1E03")) +
        scale_color_manual(values=c("gray50", "#9B1E03")) +
    #geom_vline(data=mu, aes(xintercept=grp.mean, color=DGE_sig),
    #         linetype="dashed", size=1.5) +
    theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold")) +
    stat_summary(size=2, shape=21)
    #geom_point(aes(x = grp.mean, y = DGE_sig, fill=DGE_sig),size=6, shape=21, data = mu, inherit.aes = F)


ggsave(p1, file= ("~/tonsa_epigenetics/figures/epi_plasticity_GH.pdf"), h=3.5, w=5.5)


# make violin plot:

dge_meth$plasticity <- ifelse(dge_meth$DGE_sig == TRUE, "Plastic\nGenes", "Non-plastic\nGenes")

pviolin <- ggplot(dge_meth, aes(y=hh_F25_mean_meth, x=plasticity, fill=plasticity, group= plasticity)) +
  geom_violin(alpha=0.5) + 
  #scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
  #scale_x_continuous(expand = c(0, 0)) + 
  #coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
  theme_classic() +
   #geom_quasirandom() +
        scale_fill_manual(values=c("gray50", "#9B1E03")) +
        scale_color_manual(values=c("gray50", "#9B1E03")) +
    #geom_vline(data=mu, aes(xintercept=grp.mean, color=DGE_sig),
    #         linetype="dashed", size=1.5) +
    theme(legend.position="none") +
  #stat_summary(fun = mean, geom = "crossbar", 
  #             fun.min = mean, fun.max = mean, 
  #             width = 0.3, size = 0.75, color = "black") +
    stat_summary(size=5, fun = mean, geom="point", shape=21, color="black")+
    xlab(NULL) +
    ylab("OWA % Methylation")

ggsave(pviolin, file= ("~/tonsa_epigenetics/figures/epi_plasticity_violin_GH.pdf"), h=3.5, w=3.5)

ggsave(ggarrange(p1_scatter, pviolin, labels="AUTO"), 
            file="~/tonsa_epigenetics/figures/epi_plasticity_combined_GH.pdf",
            h=3, w=6)


#### ambient

dge <- read.csv("~/tonsa_epigenetics/analysis/DGE_AAAA_vs_HHAA.txt", header=T, sep="\t")
colnames(dge) <- c("baseMean", colnames(dge)[2:ncol(dge)] )
dge$gene <- row.names(dge)

dge_meth_list <- list()
dge_meth_list[[1]] <- inner_join(dge, meth_out[['all']],      by="gene")

# calc some summary info:
dge_meth <- dge_meth_list[[1]]
dge_meth <- dge_meth[!is.na(dge_meth$padj),]

# assign sig or not:
dge_meth$DGE_sig <- ifelse(dge_meth$padj < 0.05, TRUE, FALSE)

dge_meth <- dge_meth[which(dge_meth$n_meth > 4),]
nrow(dge_meth)
#1789

mu <- ddply(dge_meth, "DGE_sig", summarise, grp.mean=mean((aa_F25_mean_meth)))
#   DGE_sig   grp.mean
# 1   FALSE 0.13348975
# 2    TRUE 0.06781858
mu2 <- ddply(dge_meth, "DGE_sig", summarise, grp.se=ster((aa_F25_mean_meth)))
#  DGE_sig      grp.se
#1   FALSE 0.005732615
#2    TRUE 0.006781181

ks.test((dge_meth$aa_F25_mean_meth[dge_meth$DGE_sig == TRUE]), 
        (dge_meth$aa_F25_mean_meth[dge_meth$DGE_sig == FALSE]))
# D = 0.16136, p-value = 8.66e-08

dge_meth$plasticity <- ifelse(dge_meth$DGE_sig == TRUE, "Plastic\nGenes", "Non-plastic\nGenes")


pviolin2 <- ggplot(dge_meth, aes(y=aa_F25_mean_meth, x=plasticity, fill=plasticity, group= plasticity)) +
  geom_violin(alpha=0.5) + 
  #scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
  #scale_x_continuous(expand = c(0, 0)) + 
  #coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
  theme_classic() +
   #geom_quasirandom() +
        scale_fill_manual(values=c("gray50", "#9B1E03")) +
        scale_color_manual(values=c("gray50", "#9B1E03")) +
    #geom_vline(data=mu, aes(xintercept=grp.mean, color=DGE_sig),
    #         linetype="dashed", size=1.5) +
    theme(legend.position="none") +
  #stat_summary(fun = mean, geom = "crossbar", 
  #             fun.min = mean, fun.max = mean, 
  #             width = 0.3, size = 0.75, color = "black") +
    stat_summary(size=5, fun = mean, geom="point", shape=21, color="black")+
    xlab(NULL) +
    ylab("Ambient % Methylation")

ggsave(pviolin2, file= ("~/tonsa_epigenetics/figures/epi_plasticity_AM.pdf"), h=3, w=4)



#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
#
# link methylaiton to expression level overall

dge <- read.csv("~/tonsa_epigenetics/analysis/DGE_norm_counts.txt", header=T, sep="\t")
dge$gene <- row.names(dge)

dge_meth_list <- list()
dge_meth_list[[1]] <- inner_join(dge, meth_out[['all']],      by="gene")
dge_meth_list[[2]] <- inner_join(dge, meth_out[['genic']],    by="gene")
dge_meth_list[[3]] <- inner_join(dge, meth_out[['exon']],     by="gene")
dge_meth_list[[4]] <- inner_join(dge, meth_out[['promoter']], by="gene")

names(dge_meth_list) <- names(snp_out)
sapply(dge_meth_list, nrow)
#     all    genic     exon promoter
#    2824     2340     1762  


dge_meth <- dge_meth_list[[3]]

# assign sig or not:

dge_meth <- dge_meth[which(dge_meth$n_meth > 5),]
nrow(dge_meth)
# 1193
#plot(dge_meth$HH_F25_1_meth, log2(dge_meth$HHHH_F1_REP1 ))

# I should calculate the mean expression
# then use these for the plots.

dge_meth$AA_mean_expression <- rowMeans(dge_meth[, c("AAAA_F1_REP1", 
                                                     "AAAA_F1_REP2", 
                                                     "AAAA_F1_REP3",
                                                     "AAAA_F1_REP4")])


summary(lm(log2(dge_meth$AA_mean_expression ) ~dge_meth$aa_F25_mean_meth))
#                          Estimate Std. Error t value Pr(>|t|)
#(Intercept)                8.86156    0.07603  116.56   <2e-16 ***
#dge_meth$aa_F25_mean_meth -3.99833    0.29128  -13.73   <2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Residual standard error: 2.245 on 1191 degrees of freedom
#Multiple R-squared:  0.1366,	Adjusted R-squared:  0.1359
#F-statistic: 188.4 on 1 and 1191 DF,  p-value: < 2.2e-16

dge_meth$log2 <- log2(dge_meth$AA_mean_expression)
p1 <- ggplot(dge_meth, aes(y=log2, x = aa_F25_mean_meth)) +
		geom_hline(yintercept=0, linetype="dashed", color="grey60") +
		geom_vline(xintercept=0, linetype="dashed", color="grey60") +
		geom_point(alpha=0.8, shape=21, fill="black")+
	    stat_smooth(method = "lm") +
	    # geom_beeswarm() +
		theme_classic() +
		ylab("log2 expression") +
		xlab("percent methylation")

#ggsave(p1, file= ("~/tonsa_epigenetics/figures/epi_dge_scatter_expression.pdf"), h=4, w=4)


########

dge <- read.csv("~/tonsa_genomics/analysis/DGE_norm_counts.txt", sep="\t")
dge$gene <- row.names(dge)
dfp <- merge(meth_out[['exon']], dge, by="gene")

fout <- dfp %>% filter(n_meth > 4)
nrow(fout)
# 1308
x=log2(fout$AAAA_F1_REP2)
y=fout$AA_F25_2_meth
df_plot <- data.frame(log2_expression_count = x, percent_methylation = y)


fout$methylation_category <- "low"
fout$methylation_category[which(fout$aa_F25_mean_meth > 0.2)] <- "high"

fout$AA_mean_expression <- log2(rowMeans(fout[,grep("AAAA_F1_REP", colnames(fout))]))

p2 <- ggplot(fout, aes(AA_mean_expression, fill=methylation_category)) +
	    geom_histogram(position="dodge")+
		theme_bw() + 
		xlab("log2 mean expression")
#p2

ggpubr::ggarrange(p1, p2, widths = c(0.4, 0.6), labels="AUTO")

ggsave("~/tonsa_epigenetics/figures/perc_meth_vs_expression_S12.png",ggpubr::ggarrange(p1, p2, widths = c(0.4, 0.6), labels="AUTO"), h=4, w=7)
ggsave("~/tonsa_epigenetics/figures/perc_meth_vs_expression_S12.pdf",ggpubr::ggarrange(p1, p2, widths = c(0.4, 0.6), labels="AUTO"), h=4, w=7)

```
