# GO enrichment:

```bash
  ~/bin/Trinotate-Trinotate-v4.0.2/util/extract_GO_assignments_from_Trinotate_xls.pl \
                         --Trinotate_xls /data/copepods/tonsa_transcriptome/tonsa_trinotate_report.txt \
                         -G \
                         > /data/copepods/tonsa_transcriptome/atonsa_trinotate_go_annotations.txt
                         # --include_ancestral_terms \


```

```r

library(dplyr)
library(tidyr)

goterms <- read.csv("/data/copepods/tonsa_transcriptome/atonsa_trinotate_go_annotations.txt", sep="\t", header=F)
colnames(goterms) <- c("gene", "go_terms")
dat <- read.table("~/tonsa_genomics/analysis/gene_level_analysis/methylation_results_table.txt", header=T)

#Separate the annotation column
dat_separated <- dat %>%
  separate_rows(annotation, sep = ";")

#Join with goterms
dat_with_go <- dat_separated %>%
  left_join(goterms, by = c("annotation" = "gene"))

#Group back and concatenate GO terms
dat_final <- dat_with_go %>%
  group_by(across(-c(annotation, go_terms))) %>%
  summarise(
    annotation = paste(unique(annotation), collapse = ";"),
    go_terms = paste(unique(go_terms[!is.na(go_terms)]), collapse = ","),
    .groups = "drop"
  ) %>%
   mutate(go_terms = na_if(go_terms, ""))


# get a "gene universe" that is in the dataset:
gene_universe <- goterms[(goterms$gene %in% unique(dat_final$annotation)),]
nrow(gene_universe)
#[1] 1840

write.table(file="~/tonsa_epigenetics/analysis/go_enrich/meth_gene_universe.txt", gene_universe,sep="\t", row.names=F, quote=FALSE, col.names=F)

# 72511 total genes in the transcriptome.
# 4043 genes in the methylation dataset
# 35014 genes annotated
# 1840 methylation genes annotated. basically the exact same ratio as the entire dataset
sum(is.na(goterms$go_terms))
length(unique(dat$annotation))

# I also want to add in the description of the genes for each. 

description <- read.csv("/data/copepods/tonsa_transcriptome/tonsa_trinotate_report.txt", sep="\t")

# result columns are: sprot_Top_BLASTP_hit, Pfam, EggNM.Description


!!!!!!!!!!!!!!!!!! NOT DONE !!!!!!!!!!!!!!!!!!!!!!!
X.gene_id

# Function to map Description based on annotation
map_description <- function(annotation) {
  genes <- str_split(annotation, ";")[[1]]
  descriptions <- lookup_description %>% filter(Gene %in% genes) %>% pull(Description)
  descriptions <- paste(descriptions, collapse = ";")
  return(descriptions)
}

# Apply the functions to add the GO_Terms and Description columns to dat
dat2 <- dat %>%
  mutate(Description = sapply(annotation, map_description),
        GO_Terms = sapply(annotation, map_go_terms))




head(as.data.frame(dat_final))

sum(!is.na(dat_final$go_terms))


genes_all <- (unique(unlist(strsplit(as.character(dat_final$annotation),split = ";"))))
# 4125

```



## TopGO


Run topGO

```r

library(topGO)
library(dplyr)

#go_terms <- read.table("~/tonsa_epigenetics/analysis/go_enrich/meth_GOterms.corrected.out", header=T)
go_terms <- read.table("~/tonsa_epigenetics/analysis/go_enrich/meth_gene_universe.txt", header=F)
colnames(go_terms) <- c("gene", "go_terms")
geneID2GO <- readMappings(file = "~/tonsa_epigenetics/analysis/go_enrich/meth_gene_universe.txt")

dat <- read.table("~/tonsa_genomics/analysis/gene_level_analysis/methylation_results_table.txt", header=T)

#dat <- read.table("~/tonsa_genomics/analysis/gene_level_analysis/methylation_results_table.txt", header=T)
#geneID2GO <- readMappings(file = "~/tonsa_epigenetics/analysis/go_enrich/meth_GOterms.topgo.out")


# identify significant genes
dat$ah_sig <- FALSE
dat$ha_sig <- FALSE
dat$hh_sig <- FALSE
dat$ah_sig[which(dat$ah_fdr_meth < 0.05 & dat$ah_large == "True") ] <- TRUE
dat$ha_sig[which(dat$ha_fdr_meth < 0.05 & dat$ha_large == "True") ] <- TRUE
dat$hh_sig[which(dat$hh_fdr_meth < 0.05 & dat$hh_large == "True") ] <- TRUE

#dat$ah_sig[which(dat$ah_fdr_meth < 0.01) ] <- TRUE
#dat$ha_sig[which(dat$ha_fdr_meth < 0.01) ] <- TRUE
#dat$hh_sig[which(dat$hh_fdr_meth < 0.01) ] <- TRUE

sum(dat$ah_sig)
# 128
sum(dat$ha_sig)
# 161
sum(dat$hh_sig)
# 753


result_ah <- dat %>%
  group_by(annotation) %>%
  summarise(has_true = any(ah_sig == "TRUE"))
result_ha <- dat %>%
  group_by(annotation) %>%
  summarise(has_true = any(ha_sig == "TRUE"))
result_hh <- dat %>%
  group_by(annotation) %>%
  summarise(has_true = any(hh_sig == "TRUE"))

sum(result_ah$has_true)
#36
sum(result_ha$has_true)
#57
sum(result_hh$has_true)
#185

# Split the annotation column by ';' and unnest the data frame
dat$annotation <- as.character(dat$annotation)
dat2 <- dat %>%
  tidyr::separate_rows(annotation, sep = ";") %>%
  filter(annotation != "-")
nrow(dat2)
#60398

sum(dat2$ah_sig)
# 77
sum(dat2$ha_sig)
# 89
sum(dat2$hh_sig)
#487

# only those with go terms:
dat3 <- dat2[dat2$annotation %in% go_terms$gene,]
nrow(dat3)
# 27023
length(unique(dat3$annotation))
# 1840
dat2 <- dat3

sum(dat2$ah_sig)
# 13
sum(dat2$ha_sig)
# 20
sum(dat2$hh_sig)
# 129
# but the above are still redundant by gene

# Group by 'annotation', summarize to find if 'TRUE' exists
result_ah <- dat2 %>%
  group_by(annotation) %>%
  summarise(has_true = any(ah_sig == "TRUE"))
result_ha <- dat2 %>%
  group_by(annotation) %>%
  summarise(has_true = any(ha_sig == "TRUE"))
result_hh <- dat2 %>%
  group_by(annotation) %>%
  summarise(has_true = any(hh_sig == "TRUE"))
sum(result_ah$has_true)
#9
sum(result_ha$has_true)
#14
sum(result_hh$has_true)
#49

head(result_ah)

toTest <- list(
  ah = result_ah,
  ha = result_ha,
  hh = result_hh
)
for(test_name in names(toTest)){
    #out <- as.data.frame(matrix(nrow=length(genes), ncol=2))
    test_set <- toTest[[test_name]]
    ofInterest <- test_set$annotation[test_set$has_true]
    out.save <- as.data.frame(matrix(ncol=8, nrow=0))
    colnames(out.save) <- c("GO.ID", "Term","Annotated","Significant","Expected","classicFisher","weight", "ontology")
    
    # set gene background
    geneUniverse <- names(geneID2GO)
    genesOfInterest <- names(geneID2GO)[names(geneID2GO) %in% ofInterest]
    
    #show genes of interest in universe vector
    geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
    names(geneList) <- geneUniverse

 for(j in c("BP", "CC", "MF")){

        myGOdata <- (new("topGOdata", description="My project",
            ontology=j, allGenes=geneList,
            annot = annFUN.gene2GO, gene2GO = geneID2GO,
            nodeSize = 5 ))
        resultWeight <- runTest(myGOdata, algorithm="weight01", statistic="fisher")
        resultClassic <- runTest(myGOdata, algorithm="classic", statistic="fisher")
        allRes <- GenTable(myGOdata, classicFisher = resultClassic,  weight = resultWeight, orderBy = "weight", ranksOf = "weight"  ,
            topNodes = length(resultWeight@score))
        allRes$ontology <- j
        all_filt <- allRes[which(allRes$weight < 0.05 & allRes$Significant > 1 ),]
        if(nrow(all_filt) > 0){


            out.save <- rbind(out.save,all_filt)

            #### get genes in go terms

            myterms =all_filt$GO.ID # change it to results.table.bh$GO.ID if working with BH corrected values
           mygenes = genesInTerm(myGOdata, myterms)

          genes_out <- as.data.frame(matrix(nrow=0, ncol=2))
          colnames(genes_out) <- c("goterm", "genes")
          for (i in 1:length(myterms)){
             myterm=myterms[i]
             mygenesforterm= mygenes[myterm][[1]]
             mygenesforterm=paste(mygenesforterm, collapse=',')
             tmp_term <- data.frame(goterm = myterm, genes = mygenesforterm)
             genes_out <- rbind(genes_out, tmp_term)
           }
        }
        write.table(genes_out,
            paste("~/tonsa_epigenetics/analysis/go_enrich/",test_name, "_", j, "_genetoGOmapping.txt", sep=""),
                sep="\t",quote=F, row.names=F)
    }
    
    write.table(file=paste("~/tonsa_epigenetics/analysis/go_enrich/", test_name, "_epigenetics_GO.txt", sep=""),
             out.save, col.names=TRUE,
            row.names=FALSE, quote=FALSE,sep="\t")
}

```
