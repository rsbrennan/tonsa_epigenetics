# tonsa epigenetics

This repository contains code to run analyses for:

`Epigenetic and evolutionary mechanisms uniquely contribute to rescue from global change. Brennan... Pespeni.`

### Data availability

Raw reads for all analyses can be found on NCBI BioProject `PRJNA590963`.

Other relevant files for the analysis can be found on Zenodo: https://doi.org/10.5281/zenodo.10797734


### scripts to run analyses

It is necessary to change paths to run these correctly. Point them to wherever you have saved the files. 

quality check: `fastqc.sh`

trim: `trim.sh`

align: `align.sh`

calculate alignment rates: `alignment_rates.md`

check for bias in methylation calls: `bias_check.sh`

Plot bias results: `bias_plot.R`

Extraction methylation calls: `methylation_extract.sh`

Merge and convert to bedgraph: `to_bedgraph.sh`

Differential methylation analysis: `differential_methylation.R`

Conversion efficiency from the phix spike in: `conversion_efficiency.md`

Compare allele frequency and methylation divergence: `methylation_vs_fst.R` 


Gene level analysis to compare methylation to DGE and allele freq divergence at the gene level: `gene_level_analysis.md`




### Figures:

- Fig. 2: `differential_methylation.R`
- Fig. 3: `methylation_vs_fst.R`
- Fig. 4: `gene_level_analysis.md` 
- Fig. S1: `differential_methylation.R`
- Fig. S2: `gene_level_analysis.md`
- Fig. S3: `gene_level_analysis.md`
- Fig. S4: `methylation_vs_fst.R`
- Fig. S5: `methylation_vs_fst.R`
- Fig. S6: `methylation_vs_fst.R`
- Fig. S7: `methylation_vs_fst.R`
- Fig. S8: `pi.md`
- Fig. S9: `gene_level_analysis.md`
- Fig. S10: `differential_methylation.R`
