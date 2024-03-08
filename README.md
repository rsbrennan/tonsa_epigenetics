# tonsa epigenetics

This repository contains code to run analyses for:

`Epigenetic and evolutionary mechanisms uniquely contribute to rescue from global change. Brennan... Pespeni.`

### Data availability

Raw reads for all analyses can be found on NCBI BioProject `PRJNA590963`.

Other relevant files for the analysis can be found on Zenodo: ADD LINK HERE


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

Compare allele frequency and methylation divergence: 


Gene level analysis to compare methylation to DGE and allele freq divergence at the gene level: `gene_level_analysis.md`




### Figures:

- Fig. 1: `differential_methylation.R`
- Fig. 2: `methylation_vs_fst.R`
- Fig. 3: `gene_level_analysis.md` 
- Fig. S1: `differential_methylation.R`
- Fig. S2: `gene_level_analysis.md`
- Fig. S3: `sliding_win_fst_epi.R`
- Fig. S4: `sliding_win_fst_epi.R`
- Fig. S5: `pi.md`
- Fig. S6: `gene_level_analysis.md`
- Fig. S7: `differential_methylation.R`
