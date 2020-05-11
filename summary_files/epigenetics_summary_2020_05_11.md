
# Conversion Rate



# mapping rates

![](figures/map_rates.png)

between the F25 samples, relatively similar mapping rates. The righer rates in F0 and F3 likely have to do with higher input dna. Those samples had more than 30 indivs.

# mbias

![](figures/mbias_all.png)

There is bias in the first 2 bases. Remove them from the methylation calls. This is standard practice in RRBS data.

# methylation statistics

loci present in all samples > 20x coverage: 
- loci: 55,297
- scaffolds: 8,816

# visualization, broad patterns

Percent methylation where boxes indicate 95% CI of the mean  
![](figures/meth_percent.png)



![](figures/mds_edgeR.png)


# differential methylation

We're broadly interseted in identifying loci where there is differential methylation between AA and the other groups. We have AA at both F0 and F25. Any locus that is truly affected by the experimental condition should be different from AA at both of these timepoints, where AA_F25 is a control for lab acclimation/selection. There are two main ways to identify these loci. We incorporate both into edgeR where all line and generations are included (methylation ~ line*generation)

1. Pairwise comparisons between specific groups and look for overlaps.
2. An averaging approach.

Of these two, the first is more stringent. We ask, where are loci significantly different from both? For second is more lenient because it is asking if, on average, a locus is different from the AA treatment. If we have a large effect in AA at one generation, but not the other, we would likely detect a significant difference.

An additional layer, which is typical in methylation studies, is requiring a minimum difference in methylation levels between the comparisons. I'll add this to both approaches, below.

### pairwise differences

#### Pairwise difference between AA F0 and other treatments:

|    Comparison   | Number significant |
|:---------------:|:------------------:|
| AA_F0 vs AA_F25 |        5970       |
| AA_F0 vs AH_F25 |       8328       |
| AA_F0 vs HA_F25 |        6172       |
| AA_F0 vs HH_F25 |       8783       |

#### Pairwise difference between AA F25 and other treatments:

|    Comparison    | Number significant |
|:----------------:|:------------------:|
| AA_F25 vs AH_F25 |         26         |
| AA_F25 vs HA_F25 |         861        |
| AA_F25 vs HH_F25 |         905       |


#### Pairwise difference between AA F25 and F0 and other treatments:

|    Comparison    | Number significant |
|:----------------:|:------------------:|
| AA_F25 and AA_F0 vs AH_F25 |           9       |
| AA_F25 and AA_F0 vs HA_F25 |           298      |
| AA_F25 and AA_F0 vs HH_F25 |           440     |

#### Pairwise difference between AA F25 and F0 and other treatments with 10 % change required:

|    Comparison    | Number significant |
|:----------------:|:------------------:|
| AA_F25 and AA_F0 vs AH_F25 |         2     |
| AA_F25 and AA_F0 vs HA_F25 |       116     |
| AA_F25 and AA_F0 vs HH_F25 |       244     |

##### overlap between the treatments

|    Comparison    | overlap |
|:----------------:|:------------------:|
| HH vs HA |         39     |
| HH vs AH |         0     |
| AH vs HA |         0     |


Note that all plots are normalized to AA_F25 mean methylation levels, so show relative change from this.

##### HH plot
![](figures/heatmap_HH_pw.png)

##### HA plot
![](figures/heatmap_HA_pw.png)

##### AH plot
![](figures/heatmap_AH_pw.png)

### Average difference

| Comparison | Number significant |
|:----------:|:------------------:|
|   AH_F25   |        2190        |
|   HA_F25   |        1988        |
|   HH_F25   |        4011        |

#### Average difference with 10 % change required:

| Comparison | Number significant |
|:----------:|:------------------:|
|   AH_F25   |        283        |
|   HA_F25   |        381        |
|   HH_F25   |        1286       |

##### overlap between the treatments

|    Comparison    | overlap |
|:----------------:|:------------------:|
| HH vs HA |         147     |
| HH vs AH |         119     |
| AH vs HA |         58     |


##### HH plot
![](figures/heatmap_HH_av.png)

##### HA plot
![](figures/heatmap_HA_av.png)

##### AH plot
![](figures/heatmap_AH_av.png)


