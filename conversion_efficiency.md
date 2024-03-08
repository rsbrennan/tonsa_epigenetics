# Lambda analysis

Checking conversion rate for the bisulfite sequencing. 

1. align to lambda genome.
2. calculate conversion rate at CpG sites.
3. % of methylated sites out of total corresponds to conversion rate. Want the % of methylated to be low, because none were methylated.

### align

Use raw reads from all samples, align to lambda genome.

NCBI reference: NC_001422.1; https://www.ncbi.nlm.nih.gov/nucleotide/NC_001422


```bash
cd ~/tonsa_epigenetics/analysis/lambda/aligned

for sample in `ls ~/tonsa_epigenetics/data/trimmed | grep '.fq.gz$' | cut -f 1-3 -d "_"| uniq`

do

    echo "starting sample ${sample}"

        ~/bin/Bismark-0.24.0/bismark --bowtie2 --multicore 1 \
        --genome ~/tonsa_epigenetics/analysis/lambda/ \
        --output_dir ~/tonsa_epigenetics/analysis/lambda/aligned \
        -1 ~/tonsa_epigenetics/data/trimmed/${sample}_1.fq.gz \
        -2 ~/tonsa_epigenetics/data/trimmed/${sample}_2.fq.gz \
        --rg_tag --rg_id ${sample} --rg_sample ${sample} --un --gzip --local --maxins 1000

     echo "finished sample ${sample}"

done

```

### extract lambda methylation

```bash

cd ~/tonsa_epigenetics/analysis/lambda/aligned

~/bin/Bismark-0.24.0/bismark_methylation_extractor --gzip --bedGraph --scaffolds \
  --cytosine_report --comprehensive \
  --ignore 2 --ignore_r2 2 \
  --no_header \
  --parallel 4 \
  --output ~/tonsa_epigenetics/analysis/lambda/methylation_extract/ \
  --genome_folder ~/tonsa_epigenetics/analysis/lambda/Bisulfite_Genome/ \
  ~/tonsa_epigenetics/analysis/lambda/aligned/*pe.bam


```

```r

library(methylKit)

setwd("~/tonsa_epigenetics/analysis/lambda/methylation_extract")

dir <- "/users/r/b/rbrennan/tonsa_epigenetics/analysis/lambda/methylation_extract"

samples <- read.table("~/tonsa_epigenetics/analysis/lambda/methylation_extract/sample_id.txt", header=FALSE)

files <- file.path(dir, samples$V1)
all(file.exists(files)) # checking that all files exist

file.list <- as.list(files)

# get the names only for naming our samples
nmlist <- as.list(gsub("_1_bismark_bt2_pe.bismark.cov.gz","",samples$V1))

# use methRead to read in the coverage files
myobj <- methRead(location= file.list,
        sample.id =   nmlist,
                      assembly = "lambda", # this is just a string. no actual database
                      dbtype = "tabix",
                      context = "CpG",
                      resolution = "base",
                      mincov = 10,
                            treatment = 
                              c(0,0,0,0,
                                0,0,0,0,
                                1,1,1,1,
                                2,2,2,2,
                                3,3,3,3,
                                4,4,4,4),
                      pipeline = "bismarkCoverage",
                      dbdir = "~/tonsa_epigenetics/analysis/lambda/methylation_extract")

myobj.filt <- filterByCoverage(myobj,
                      lo.count=10,
                      lo.perc=NULL,
                      hi.count=NULL,
                      hi.perc=99.9)

getMethylationStats(myobj[[1]],plot=TRUE,both.strands=FALSE)
str(getMethylationStats(myobj[[1]],plot=FALSE,both.strands=FALSE))
meth <- unite(myobj, destrand=FALSE)
pm=percMethylation(meth)
colMeans(pm)

out <- as.data.frame(matrix(nrow=nrow(samples), ncol=2))
colnames(out) <- c("sample", "conversion_efficiency")
out$sample <- gsub("_1_bismark_bt2_pe.bismark.cov.gz","",samples$V1)

for(i in 1:nrow(samples)){

  out$conversion_efficiency[i] <- round(mean(getData(myobj[[i]])$numTs/getData(myobj[[i]])$coverage),3)

}
mean(out$conversion_efficiency)
# 0.9841667
sd(out$conversion_efficiency)
# 0.003559026
min(out$conversion_efficiency)
# 0.976

write.table(file = "~/tonsa_epigenetics/analysis/conversion_efficiency.txt", out, row.names=F, quote=F)

```
