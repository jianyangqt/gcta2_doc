## Genomic risk prediction

### SBLUP
This is a method to run BLUP analysis using summary data from a GWAS/meta-analysis with LD from a reference sample with individual-level data. Similar methods have been proposed in recent studies to predict complex traits and diseases using GWAS summary data (Vilhjálmsson et al. 2015 AJHG; Robinson et al. 2017 Nat Hum Behav) or age using summary data from transcriptome-wise association studies (Peters et al. 2015 Nat Comm).

> Example
```bash
gcta64 --bfile test --cojo-file test.ma --cojo-sblup 1.33e6 --cojo-wind 1000 --thread-num 20
```

--cojo-file test.ma  
Input file in COJO format (see COJO)

--cojo-sblup  
Perform COJO-SBLUP analysis. The input parameter lambda = *m* \* (1 / *h*<sup>2</sup><sub>SNP</sub> - 1) where *m* is the total number of SNPs used in this analysis (i.e. the number of SNPs in common between the summary data and the reference set), and *h*<sup>2</sup><sub>SNP</sub> is the proportion of variance in the phenotype explained by all SNPs. *h*<sup>2</sup><sub>SNP</sub> can be estimated from GCTA-GREML if individual-level data are available or from LD score regression analysis of the summary data. 

**Note:** for the ease of computation, the analysis can be performed separately on individual chromosomes with the same lambda value calculated from all the genome-wide SNPs. 

--cojo-wind 1000  
Specify a distance d (in Kb unit). LD between SNPs more than d Kb away from each other are ignored. The default value is 10000 Kb (i.e. 10 Mb) if not specified. We recommend a window size of 1 Mb for the ease of computation.

> Output file format
```nohighlight
rs10055084    T     0.0094     0.00046175
rs10076494    T     0.0074     0.000361283
rs10057531    T     -0.0075    -0.0003681
rs10039735    T     -0.0095    -0.000469503
rs1507712     A     0.0091     0.000438602
rs6869386     T     -0.0068    -3.35496e-05
rs7734346     T     0.0106     0.000133413
rs10065373    T     0.0093     -1.55927e-05
rs10070362    T     -0.0108    -0.000150901
```
> Columns are SNP, the coded allele, effect size in the original GWAS summary data, and BLUP estimate of the SNP effect (all SNPs are fitted jointly).

**Note:** Let b = per-allele effect size and u = effec size per standardized SNP genotype. GCTA-SBLUP reads b, fits the model based on u, and output the BLUP SNP effect in the scale of b. So, the output can directly be used to compute the profile score using PLINK --score.


#### Citations
**COJO-SBLUP method**: Robinson et al. (2017) Genetic evidence of assortative mating in humans. Nat Hum Behav, 1:0016.

**GCTA software**: Yang J, Lee SH, Goddard ME and Visscher PM (2011) GCTA: a tool for Genome-wide Complex Trait Analysis. Am J Hum Genet, 88: 76-82. [PubMed ID: 21167468]

