## Mendelian randomisation

### GSMR

**GSMR: Generalised Summary-data-based Mendelian Randomisation**

The GSMR method tests for putative causal association between a risk factor and a disease using summary-level data from genome-wide association studies (GWAS). Details of the method can be found in Zhu et al. ([2018 Nat. Commun.](https://www.nature.com/articles/s41467-017-02317-2)). The corresponding R package is at [http://cnsgenomics.com/software/gsmr/](http://cnsgenomics.com/software/gsmr/). This GCTA module runs the same analysis as in the R version of GSMR but is meant to be faster and more flexible.

> Example
```bash
gcta64 --mbfile gsmr_ref_data.txt --gsmr-file gsmr_exposure.txt gsmr_outcome.txt --gsmr-direction 0 --out test_gsmr_result
```

--mbfile gsmr\_ref\_data.txt  

The GSMR analysis requires a reference sample with individual level genotypes (in PLINK binary format) for LD estimation. If the genotype data are very large, the data are often saved in separate PLINK files (e.g. one for each chromosome). Here, we provide an option to read GWAS genotype data saved in multiple PLINK files. The input is a text file with each row representing a PLINK binary file.
> Input file format
```nohighlight
gsmr_ref_data_chr1
gsmr_ref_data_chr2
…
```
**Note:** the option --bfile is still valid when there is only a single PLINK file.

--gsmr-file gsmr\_exposure.txt gsmr\_outcome.txt  
The inputs are two text or compressed text (gz format) files containing the filepaths of the GWAS summary data. The first one is for exposures and the second one is for outcomes. 
> Input file format 
 
gsmr\_exposure.txt
```nohighlight
hdl hdl_test.raw
bmi bmi_test.raw
```
gsmr\_outcome.txt
```nohighlight
t2d t2d_test.raw
```
Columns are the trait name, filepath of the GWAS summary data. Each row represents a trait.
> Format of the GWAS summary data (i.e. the [GCTA-COJO format](#COJO))  

bmi\_test.raw  
```nohighlight
SNP A1  A2  freq    b   se  p   N
rs1000000   G   A   0.781838245 1.00E-04    0.0044  0.9819  231410
rs10000010  T   C   0.513760872 -0.0029 0.003   0.3374  322079
rs10000012  G   C   0.137219265 -0.0095 0.0054  0.07853 233933
rs10000013  A   C   0.775931455 -0.0095 0.0044  0.03084 233886
```
Columns are SNP, the effect allele, the other allele, frequency of the effect allele, effect size, standard error, p-value and sample size.

--gsmr-direction 0   
There are 3 GSMR analyses, forward-GSMR analysis (coded as 0), reverse-GSMR analysis (coded as 1) and bi-GSMR analysis (both forward- and reverse-GSMR analyses, coded as 2). 

--out test\_gsmr\_result  
> Output file format  

test\_gsmr\_result.gsmr  
```nohighlight
exposure  outcome      bxy          se          pval       nsnps
bmi         t2d      0.798596    0.086785     3.51315e-20   77
hdl         t2d     -0.125651    0.0431615    0.00360073    130
```
Columns are exposure, outcome, GSMR estimates of *b*<sub>xy</sub>, standard error, p-value and number of SNPs.

#### Visualization
--effect-plot
This flag will save the information required to generate a plot of SNP effects as in Figure 3d of Zhu et al. ([2018 Nature Communications](https://www.nature.com/articles/s41467-017-02317-2); also see the figure below) in a compressed text file (\*.eff\_plot.gz). We provide an R script ([gsmr\_plot.r](./res/gsmr_plot.r)) to generate the effect size plot based on the \*.eff\_plot.gz file (see the example below).

> Example

```r
source("gsmr_plot.r")
gsmr_data = read_gsmr_data("test_gsmr_result.eff_plot.gz")
gsmr_summary(gsmr_data)      # show a summary of the data
plot_gsmr_effect(gsmr_data, "bmi", "t2d", colors()[75])           
```

![effect_size_plot](./res/gsmr_toy_bmi_t2d.jpg)

#### Optional flags  
--gsmr2-beta to use the new HEIDI-outlier method.   
Note: We included a new HEIDI-outlier method (as part of [the GSMR analysis](#Download)) in GCTA v1.91.7. However, the new HEIDI-outlier method is currently under development and subject to changes during the method development. From GCTA version v1.92.0, we changed the default back to the original HEIDI-outlier method described in Zhu et al. (2018 Nature Communications) and added this temporary flag to test the new method. The new HEIDI-outlier method in GCTA v1.92.0 has been tested by extensive simulations and real data analyses. We will make a formal release in our next GSMR paper.

*Quality control*

--diff-freq 0.2   
To check the difference in allele frequency of each SNP between the GWAS summary datasets and the LD reference sample. SNPs with allele frequency differences greater than the specified threshold value will be excluded from the analysis. The default value is 0.2.

*Clumping analysis*  

--gwas-thresh 5e-8   
To specify a threshold p-value to select SNPs for clumping. Note that the SNP instruments used in the GSMR analysis ([Zhu et al. 2018 Nat. Commun.](https://www.nature.com/articles/s41467-017-02317-2)) are a set of near-independent SNPs from the clumping analysis at this threshold. The default value is 5e-8.

--clump-r2 0.05  
LD *r*<sup>2</sup> threshold for clumping analysis. The default value is 0.05. 

*GSMR analysis*  

--heidi-thresh 0.01  
The HEIDI-outlier method described in Zhu et al. ([2018 Nature Communications](https://www.nature.com/articles/s41467-017-02317-2))   
To specific a p-value threshold for the HEIDI-outlier analysis to remove horizontal pleiotropic SNPs. The default threshold is 0.01.

--heidi-thresh 0.01 0.01 (to be used together with --gsmr2-beta)   
The input parameters are two p-value thresholds used in the new HEIDI-outlier method.   
The new HEIDI-outlier method involves two steps: 1) single-SNP-based HEIDI-outlier analysis, excluding a SNP with the smallest p-value from a single-SNP-based HEIDI-outlier analysis iteratively until the HEIDI-outlier p-values of all the remaining SNPs are not smaller than the specified p-value threshold (the first input parameter); 2) multi-SNP-based HEIDI-outlier, excluding a SNP with the smallest p-value (from a single-SNP-based HEIDI-outlier analysis) iteratively until the p-value from the multi-SNP-based HEIDI-outlier test is not smaller than the specified p-value threshold (the second input parameter). The default value is 0.01 for both.

--gsmr-snp-min 10   
To specify the minimum number of genome-wide significant and near-independent SNPs required for the GSMR analysis. Note that the SNP instruments will be pruned for LD by a clump analysis and filtered for horizontal pleiotropy by the HEIDI-outlier analysis. This option will count the number of SNPs after clumping and HEIDI-outlier filtering. The default value is 10.

--gsmr-ld-fdr 0.05  
FDR threshold to shrink the chance correlations between the SNP instruments to zero. The default value is 0.05. If the reference sample is independent from the GWAS samples, it is not valid to approximate the chance correlations between SNPs in the GWAS data by those estimated from the reference sample. Under the null that the SNPs are not correlated, *n**r*<sup>2</sup> follows a chi-squared distribution with df = 1, where n is the sample size and df is the degrees of freedom.  

#### Citation  
Zhu, Z. et al. (2018) Causal associations between risk factors and common diseases inferred from GWAS summary data. Nat. Commun. 9, 224.  



