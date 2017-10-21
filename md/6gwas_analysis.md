
## GWAS analysis

### COJO

**GCTA-COJO: multi-SNP-based conditional & joint association analysis using GWAS summary data**

--cojo-file test.ma  
Input the summary-level statistics from a meta-analysis GWAS (or a single GWAS).
> Input file format  
> test.ma
```nohighlight
SNP A1 A2 freq b se p N 
rs1001 A G 0.8493 0.0024 0.0055 0.6653 129850 
rs1002 C G 0.0306 0.0034 0.0115 0.7659 129799 
rs1003 A C 0.5128 0.0045 0.0038 0.2319 129830
...
```

Columns are SNP, the effect allele, the other allele, frequency of the effect allele, effect size, standard error, p-value and sample size. The headers are not keywords and will be omitted by the program. Important: "A1" needs to be the effect allele with "A2" being the other allele and "freq" should be the frequency of "A1".

Note: 1) For a case-control study, the effect size should be log(odds ratio) with its corresponding standard error. 2) Please always input the summary statistics of all the SNPs even if your analysis only focuses on a subset of SNPs because the program needs the summary data of all SNPs to calculate the phenotypic variance.

--cojo-slct  
Perform a stepwise model selection procedure to select independently associated SNPs. Results will be saved in a *.jma file with additional file 
*.jma.ldr showing the LD correlations between the SNPs.

--cojo-top-SNPs 10  
Perform a stepwise model selection procedure to select a fixed number of independently associated SNPs without a p-value threshold. The output format is the same as that from --cojo-slct.

--cojo-joint  
Fit all the included SNPs to estimate their joint effects without model selection. Results will be saved in a *.jma file with additional file *.jma.ldr showing the LD correlations between the SNPs.

--cojo-cond cond.snplist  
Perform association analysis of the included SNPs conditional on the given list of SNPs. Results will be saved in a *.cma.
> Input file format  
> cond.snplist
```nohighlight
rs1001
rs1002
...
```

--cojo-p 5e-8  
Threshold p-value to declare a genome-wide significant hit. The default value is 5e-8 if not specified. This option is only valid in conjunction with the option --cojo-slct. 

Note: it will be extremely time-consuming if you set a very low significance level, e.g. 5e-3.

--cojo-wind 10000  
Specify a distance d (in Kb unit). It is assumed that SNPs more than d Kb away from each other are in complete linkage equilibrium. The default value is 10000 Kb (i.e. 10 Mb) if not specified.

--cojo-collinear 0.9  
During the model selection procedure, the program will check the collinearity between the SNPs that have already been selected and a SNP to be tested. The testing SNP will not be selected if its multiple regression *R<sup>2</sup>* on the selected SNPs is greater than the cutoff value. By default, the cutoff value is 0.9 if not specified.

--cojo-gc  
If this option is specified, p-values will be adjusted by the genomic control method. By default, the genomic inflation factor will be calculated from the summary-level statistics of all the SNPs unless you specify a value, e.g. --cojo-gc 1.05.

--cojo-actual-geno  
If the individual-level genotype data of the discovery set are available (e.g. a single-cohort GWAS), you can use the discovery set as the reference sample. In this case, the analysis will be equivalent to a multiple regression analysis with the actual genotype and phenotype data. Once this option is specified, GCTA will take all pairwise LD correlations between all SNPs into account, which overrides the –cojo-wind option. This option also allows GCTA to calculate the variance taken out from the residual variance by all the significant SNPs in the model, otherwise the residual variance will be fixed constant at the same level of the phenotypic variance.

Examples (Individual-level genotype data of the discovery set is NOT available) - Robust and recommended
```bash
# Select multiple associated SNPs through a stepwise selection procedure
gcta64  --bfile test  --chr 1 --maf 0.01 --cojo-file test.ma --cojo-slct --out test_chr1

# Select a fixed number of of top associated SNPs through a stepwise selection procedure
gcta64  --bfile test  --chr 1 --maf 0.01 --cojo-file test.ma --cojo-top-SNPs 10 --out test_chr1

# Estimate the joint effects of a subset of SNPs (given in the file test.snplist) without model selection
gcta64  --bfile test  --chr 1 --extract test.snplist  --cojo-file test.ma --cojo-joint --out test_chr1

# Perform single-SNP association analyses conditional on a set of SNPs (given in the file cond.snplist) without model selection
gcta64  --bfile test  --chr 1 --maf 0.01 --cojo-file test.ma --cojo-cond cond.snplist --out test_chr1
```

It should be more efficient to separate the analysis onto individual chromosomes or even some particular genomic regions. Please refer to the Data management section for some other options, e.g. including or excluding a list of SNPs and individuals or filtering SNPs based on the imputation quality score.

Examples (Individual-level genotype data of the discovery set is available)
```bash
# Select multiple associated SNPs through a stepwise selection procedure
gcta64  --bfile test  --maf 0.01 --cojo-file test.ma --cojo-slct --cojo-actual-geno --out test

# In this case, it is recommended to perform the analysis using the data of all the genome-wide SNPs rather than separate the analysis onto individual chromosomes because GCTA needs to calculate the variance taken out from the residual variance by all the significant SNPs in the model, which could give you a bit more power.
# Estimate the joint effects of a subset of SNPs (given in the file test.snplist) without model selection
gcta64  --bfile test  --extract test.snplist  --cojo-file test.ma --cojo-actual-geno  --cojo-joint --out test

# Perform single-SNP association analyses conditional on a set of SNPs (given in the file cond.snplist) without model selection
gcta64  --bfile test  --maf 0.01 --cojo-file test.ma --cojo-actual-geno  --cojo-cond cond.snplist --out test
```

> Output file format  
> test.jma (generate by the option --cojo-slct or --cojo-joint)
```nohighlight
Chr SNP bp freq refA b se p n freq_geno bJ bJ_se pJ LD_r
1 rs2001 172585028 0.6105 A 0.0377 0.0042 6.38e-19 121056 0.614 0.0379 0.0042 1.74e-19 -0.345
1 rs2002 174763990 0.4294 C 0.0287 0.0041 3.65e-12 124061 0.418 0.0289 0.0041 1.58e-12 0.012
1 rs2003 196696685 0.5863 T 0.0237 0.0042 1.38e-08 116314 0.589 0.0237 0.0042 1.67e-08 0.0 
...
```

Columns are chromosome; SNP; physical position; frequency of the effect allele in the original data; the effect allele; effect size, standard error and p-value from the original GWAS or meta-analysis; estimated effective sample size; frequency of the effect allele in the reference sample; effect size, standard error and p-value from a joint analysis of all the selected SNPs; LD correlation between the SNP i and SNP i + 1 for the SNPs on the list.

> LD correlation matrix between all pairwise SNPs listed in test.jma.  
> test.jma.ldr (generate by the option --cojo-slct or --cojo-joint)
```nohighlight
SNP rs2001 rs2002 rs2003 ...
rs2001 1 0.0525 -0.0672 ...
rs2002 0.0525 1 0.0045 ...
rs2003 -0.0672 0.0045 1 ...
...
```

> test.cma (generate by the option --cojo-slct or --cojo-cond)
```nohighlight
Chr SNP bp freq refA b se p n freq_geno bC bC_se pC
1 rs2001 172585028 0.6105 A 0.0377 0.0042 6.38e-19 121056 0.614 0.0379 0.0042 1.74e-19
1 rs2002 174763990 0.4294 C 0.0287 0.0041 3.65e-12 124061 0.418 0.0289 0.0041 1.58e-12
1 rs2003 196696685 0.5863 T 0.0237 0.0042 1.38e-08 116314 0.589 0.0237 0.0042 1.67e-08 
...
```

Columns are chromosome; SNP; physical position; frequency of the effect allele in the original data; the effect allele; effect size, standard error and p-value from the original GWAS or meta-analysis; estimated effective sample size; frequency of the effect allele in the reference sample; effect size, standard error and p-value from conditional analyses.

#### The choice of reference sample for GCTA-COJO analysis

1) If the summary data are from a single cohort based GWAS, the best reference sample is the GWAS sample itself.  
2) For a meta-analysis where individual-level genotype data are not available, you could use one of the large participating cohorts. For example, Wood et al. 2014 Nat Genet used the ARIC cohort (data available from dbGaP).  
3) We suggest you use a reference sample with a sample size > 4000 (see Supplementary Figure 4 of Yang et al. 2012 Nat Genet).  
4) We do NOT suggest you use HapMap or 1000G panels as the reference sample. The sample sizes of HapMap and 1000G are not large enough.

#### GCTA-COJO analysis conditioning on a single SNP

1) create a file including the SNP ID.  
For example, cond.snplist)
```nohighlight
rs1001
```

2) then run 
```bash
gcta64  --bfile test  --cojo-file test.ma --cojo-cond cond.snplist --out test
```

#### References

**Conditional and joint analysis method**: Yang et al. (2012) Conditional and joint multiple-SNP analysis of GWAS summary statistics identifies additional variants influencing complex traits. Nat Genet 44(4):369-375. [PubMed ID: 22426310]

**GCTA software**: Yang J, Lee SH, Goddard ME and Visscher PM. GCTA: a tool for Genome-wide Complex Trait Analysis. Am J Hum Genet. 2011 Jan 88(1): 76-82. [PubMed ID: 21167468]

### mtCOJO

**mtCOJO: multi-trait-based conditional & joint analysis using GWAS summary data**

If you have two phenotypes (y and x, which can be measured on two different samples), and you want to run GWAS analysis for y conditioning on x, this analysis can be achieved by a mtCOJO analysis using GWAS summary-level data for y and x. Details of the method can be found in the Zhu et al. ([2017 bioRxiv](https://www.biorxiv.org/content/early/2017/07/26/168674)) paper. In the mtCOJO analysis, we used the GSMR method to estimate the effect of x on y, the mtCOJO estimate is free of the collider bias as described in Aschard et al. ([2015 AJHG](http://www.sciencedirect.com/science/article/pii/S0002929714005278?via%3Dihub)). It should be noted that mtCOJO allows for an analysis of y conditioning on multiple x.

> Example
```bash
gcta64 --bfile mtcojo_ref_data --mtcojo-file mtcojo_summary_data.list --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out test_mtcojo_result
```

--mtcojo-file mtcojo\_summary\_data.list  
Reading a list that contains filepaths of the GWAS summary data and prevalence of diseases.

> Input file format
```nohighlight
t2d t2d_test.raw 0.176306984 0.09
bmi bmi_test.raw
```

Columns are the trait name, filepath of the GWAS summary data, sample prevalence and population prevalence. Each row represents a trait. The first row is for the target trait (i.e. y), and the remaining rows are for covariate traits.

**Note:** If the sample prevalence and the population prevalence are not specified, the estimate of the SNP-based h<sup>2</sup> will be on the observed scale.

> Format of the GWAS summary data (i.e. the [GCTA-COJO format](http://cnsgenomics.com/software/gcta/#COJO))

bmi\_test.raw
```nohighlight
SNP A1  A2  freq    b   se  p   N
rs1000000   G   A   0.781838245 1.00E-04    0.0044  0.9819  231410
rs10000010  T   C   0.513760872 -0.0029 0.003   0.3374  322079
rs10000012  G   C   0.137219265 -0.0095 0.0054  0.07853 233933
rs10000013  A   C   0.775931455 -0.0095 0.0044  0.03084 233886
```
Columns are SNP, the effect allele, the other allele, frequency of the effect allele, effect size, standard error, p-value and sample size.

--ref-ld-chr eur\_w\_ld\_chr/  
The directory of LD score files (the same format as in [LDSC software tool](https://github.com/bulik/ldsc)). Note that the LD scores will be used in the LD score regression analysis to estimate the SNP-based heritability for a trait based on summary data, genetic correlation between two traits, and to estimate the degree of overlap between two samples (see [Zhu et al. 2017 bioRxiv](https://www.biorxiv.org/content/early/2017/07/26/168674)).

--w-ld-chr eur\_w\_ld\_chr/  
The directory of LD scores for the regression weights (the same format as in [LDSC software tool](https://github.com/bulik/ldsc)).

--out test\_mtcojo\_result    
> Output file format  

test\_mtcojo\_result.mtcojo.cma
```nohighlight
SNP A1  A2  freq    b   se  p   N   bC  bC_se   bC_pval
rs4040617   G   A   0.143126    -0.0295588  0.0397075   0.50    28944   -0.0343525  0.0399332   0.389652
rs6687776   T   C   0.133909    0.0295588   0.0343927   0.35    38288   0.0217897   0.0344916   0.527557
```
Columns are SNP, effect allele, the other allele, frequency of the effect allele, effect size, standard error, p-value, sample size from the original GWAS summary data, mtCOJO effect size, mtCOJO standard error and mtCOJO p-value.

#### Optional flags

*Clumping analysis*

--clump-p1 5e-8  
P-value threshold for index SNPs. The default threshold is 5e-8.

--clump-kb 10000  
Window size for clumping analysis. The default distance is 10,000 kb (i.e. 10 Mb).

--clump-r2 0.05  
LD r<sup>2</sup> threshold for clumping analysis. The default value is 0.05.

*GSMR analysis*

--gwas-thresh 5e-8  
P-value threshold to select instruments for the GSMR analysis (see [Zhu et al. 2017 bioRxiv](https://www.biorxiv.org/content/early/2017/07/26/168674)). Instruments are filtered from the index SNPs. The default threshold is 5e-8.

--heidi-thresh 0.01  
P-value threshold for the HEIDI-outlier test to filter instruments. The default threshold is 0.01.

--heidi-snp 10  
The minimum number of SNP instruments for the HEIDI-outlier test. The default number is 10.

--gsmr-snp 10  
The minimum number of SNP instruments for the GSMR analysis. The default number is 10.


### Mixed model association

**GCTA-MLMA and GCTA-LOCO: mixed linear model based association analysis** 

The following options are designed to perform an MLM based association analysis. Previous data management options such as --keep, --extract and --maf, REML analysis options such as --reml-priors, --reml-maxit and --reml-no-constrain and multi-threading option --thread-num are still valid for this analysis.

--mlma  
This option will initiate an MLM based association analysis including the candidate SNP  
*y = a + bx + g + e*  
where *y* is the phenotype, *a* is the mean term, *b* is the additive effect (fixed effect) of the candidate SNP to be tested for association, *x* is the SNP genotype indicator variable coded as 0, 1 or 2, g is the polygenic effect (random effect) i.e. the accumulated effect of all SNPs (as captured by the GRM calculated using all SNPs) and *e* is the residual. For the ease of computation, the genetic variance, *var(g)*, is estimated based on the null model i.e. *y = a  + g + e* and then fixed while testing for the association between each SNP and the trait. This analysis would be similar as that implemented in other software tools such as EMMAX, FaST-LMM and GEMMA.
The results will be saved in the *.mlma file.
 
--mlma-loco  
This option will implement an MLM based association analysis with the chromosome, on which the candidate SNP is located, excluded from calculating the GRM. We call it MLM leaving-one-chromosome-out (LOCO) analysis. The model is  
*y = a + bx + g<sup>-</sup> + e*  
where *g<sup>-</sup>* is the accumulated effect of all SNPs except those on the chromosome where the candidate SNP is located. The *var(g<sup>-</sup>)* will be re-estimated each time when a chromosome is excluded from calculating the GRM. The MLM-LOCO analysis is computationally less efficient but more powerful as compared with the MLM analysis including the candidate (--mlma).
The results will be saved in the *.loco.mlma file. **Note: not recommended for data with related individuals.** 
 
--mlma-no-adj-covar  
If there are covariates included in the analysis, the covariates will be fitted in the null model, a model including the mean term (fixed effect), covariates (fixed effects), polygenic effects (random effects) and residuals (random effects). By default, in order to improve computational efficiency, the phenotype will be adjusted by the mean and covariates, and the adjusted phenotype will subsequently be used for testing SNP association. However, if SNPs are correlated with the covariates, pre-adjusting the phenotype by the covariates will probably cause loss of power. If this option is specified, the covariates will be fitted together with the SNP for association test. However, this will significantly reduce computational efficiency.

--mlma-subtract-grm  
Subtract a GRM for a subset of SNPs (e.g. calculated from SNPs on one chromosome) from that for all the SNPs. This option is designed to parallelise the MLMA-LOCO analysis for large data set. Please see the example below.
 
> Examples
```bash
# MLMA analysis - If you have already computed the GRM
gcta64 --mlma --bfile test --grm test --pheno test.phen --out test --thread-num 10

# MLMA analysis using multiple GRMs - If you have already computed the GRM
gcta64 --mlma --bfile test --mgrm multi_grm.txt --pheno test.phen --out test --thread-num 10 
 
# MLMA analysis including the candidate SNP (MLMi)
gcta64 --mlma --bfile test --pheno test.phen --out test --thread-num 10

# MLMA leaving-one-chromosome-out (LOCO) analysis
gcta64 --mlma-loco --bfile test --pheno test.phen --out test --thread-num 10

# MLMA-LOCO analysis for large data sets
gcta64 --mlma --grm test_all --mlma-subtract-grm test_chr1 --bfile test --chr 1 --pheno test.phen --out test_loco_chr1 --thread-num 10
gcta64 --mlma --grm test_all --mlma-subtract-grm test_chr2 --bfile test --chr 2 --pheno test.phen --out test_loco_chr2 --thread-num 10
...
gcta64 --mlma --grm test_all --mlma-subtract-grm test_chr22 --bfile test --chr 22 --pheno test.phen --out test_loco_chr22 --thread-num 10

```

> Note: test\_all is the GRM calculated from all SNPs; test_chr1 is the GRM calculated from SNPs on chromosome 1.
 
> Output file format  
> test.mlma or test.loco.mlma (columns are chromosome, SNP, physical position, reference allele (the coded effect allele), the other allele, frequency of the reference allele, SNP effect, standard error and p-value).
```nohighlight
Chr SNP     bp    ReferenceAllele  OtherAllele Freq   b           se          p
1   qtl2_1  1001  L                H           0.366  0.0143857   0.0411682   0.726761
1   qtl2_2  1002  H                L           0.326  -0.0240756  0.0421248   0.56764
1   qtl2_3  1003  H                L           0.146  -0.0921772  0.0565541   0.103124
1   qtl2_4  1004  H                L           0.3865 -0.0771376  0.0394826   0.0507357
1   qtl2_5  1005  H                L           0.1665 0.00251276  0.0526821   0.961958
1   qtl2_6  1006  L                H           0.119  -0.0153568  0.059891    0.797632
1   qtl2_7  1007  L                H           0.1675 -0.0487809  0.0512279   0.340979
...
``` 
 
#### References
 
**An overview of the MLM based association methods**: Yang J, Zaitlen NA, Goddard ME, Visscher PM and Price AL (2014) Mixed model association methods: advantages and pitfalls. Nat Genet. 2014 Feb;46(2):100-6. [Pubmed ID: 24473328]

**REML analysis and GCTA Software**: Yang J, Lee SH, Goddard ME and Visscher PM. GCTA: a tool for Genome-wide Complex Trait Analysis. Am J Hum Genet. 2011 Jan 88(1): 76-82. [PubMed ID: 21167468]

### Gene-based test

**GCTA-fastBAT: a fast and flexible gene- or set-Based Association Test using GWAS summary data**

This method performs a fast set-based association analysis for human complex traits using summary-level data from genome-wide association studies (GWAS) and linkage disequilibrium (LD) data from a reference sample with individual-level genotypes. Please see [Bakshi et al. (2016 Scientific Reports)](https://www.nature.com/articles/srep32894) for details about the method. This module is developed by Andrew Bakshi and Jian Yang.

**Note**: most other GCTA options are also valid in this analysis.

> Examples
```bash
# Gene-based test
gcta64 --bfile test --maf 0.01 --fastBAT assoc.txt --fastBAT-gene-list gene_list.txt --out test --thread-num 10

# Segment-based test (size of a segment = 100Kb)
gcta64 --bfile test --maf 0.01 --fastBAT assoc.txt --fastBAT-seg 100 --out test --thread-num 10

# Set-based test with a customized set file (note that this can be used to test all SNPs involved in a pathway)
gcta64 --bfile test --maf 0.01 --fastBAT assoc.txt --fastBAT-set-list set.txt --out test --thread-num 10
```

#### Options 

--bfile test  
Input SNP genotype data (in PLINK binary PED format) as the reference set for LD estimation. For a single-cohort based GWAS, the GWAS cohort itself can be used as the reference set. For a meta-analysis, you can use one of the largest participating cohorts as the reference set. If none of them are available, you might use data from the 1000 Genomes Project (you will need PLINK2 --vcf option to convert the data into PLINK binary PED format). Please see Figure 1 of Bakshi et al. 2016 for a comparison of results using different reference sets for LD. 

--fastBAT assoc.txt  
Input association p-values of a list of SNPs. This can be from a standard GWAS or from a meta-analysis.
> Input file format  
> assoc.txt
```nohighlight
SNP         p
rs1001      0.0055
rs1002      0.0115
......
```

--fastBAT-gene-list gene_list.txt  
Input gene list with gene start and end positions.
> Input file format  
> gene_list.txt (columns are gene ID, chromosome, left and right side boundary of the gene region)  
```nohighlight
Chr   Start       End         Gene
1     19774       19899       Gene1
1     34627       35558       Gene2
......
```

Please click the link below to download the gene list file.
> Gene list (hg18): [glist-hg18.txt](./glist-hg18.txt)  
> Gene list (hg19): [glist-hg19.txt](./glist-hg19.txt)  

--fastBAT-set-list set.txt  
Input set list with name and list of SNPs in the set.
> Input file format  
> set.txt (set ID, followed by SNPs, then END, then blank space before next set)  
```nohighlight
Set1
rs1234534
rs5827743
rs9737542
END

Set2    
rs1252514
...
```

This option provides an opportunity for you to customize your own sets of SNPs. For example, you can create a SNP set which contains all the 1KGP SNPs in genes involved in a pathway listed in the file below.
> pathway list: [c2.cp.v5.1.symbols.gmt](./c2.cp.v5.1.symbols.gmt) (downloaded from [Broad GSEA](http://software.broadinstitute.org/gsea/index.jsp))

--fastBAT-seg 100  
Perform fastBAT analysis based on segments of size 100Kb (default).

#### Other options 

--fastBAT-wind 50  
Used in conjunction with --fastBAT-gene-list to define a gene region. By default, a gene region is defined as +-50kb of UTRs of a gene.

--fastBAT-ld-cutoff 0.9  
Threshold LD r-squared value for LD pruning. The default value is 0.9. You can turn off LD pruning by setting this value to 1.

--fastBAT-write-snpset   
Write the sets of SNPs included in the analysis. The SNP sets will be saved in a text file in the same format as the input file of --fastBAT-set-list. 

> Output file format  
> Possible output file names:  
```nohighlight
test.fastbat (set-based test)
test.seg.fastbat (segment-based test)
test.gene.fastbat (gene-based test)
```
> test.gene.fbat (columns are)
```nohighlight
Gene: gene ID  
Chr: chromosome  
Start and End: left and right side boundaries of the gene region
No.SNPs: number of SNPs in the gene region
SNP_start and SNP_end: the SNP at the left and right side boundary of the gene region
Chisq(Obs): sum of chi-squared test-statstics of all SNPs in the gene region
Pvalue: gene-based test p-value
TopSNP.Pvalue: smallest single-SNP GWAS p-value in the gene region
TopSNP: the top associated GWAS SNP
```

> test.seg.fbat (columns are)
```nohighlight
Chr: chromosome
Start and End: left and right side boundaries of the segment
No.SNPs: number of SNPs in the gene region
SNP_start and SNP_end: the SNP at the left and right side boundary of the gene region
Chisq(Obs): sum of chi-squared test-statstics of all SNPs in the segment
Pvalue: segment-based test p-value
TopSNP.Pvalue: smallest single-SNP GWAS p-value in the segment
TopSNP: the top associated GWAS SNP
```

> test.fbat (columns are)
```nohighlight
Set: set ID
No.SNPs: number of SNPs in the gene region
SNP_start and SNP_end: the SNP at the left and right side boundary of the gene region
Chisq(Obs): sum of chi-squared test-statstics of all SNPs in the set
Pvalue: segment-based test p-value
TopSNP.Pvalue: smallest single-SNP GWAS p-value in the segment
TopSNP: the top associated GWAS SNP
```

#### References:

**fastBAT method**: Bakshi A., Zhu Z., Vinkhuyzen A.A.E., Hill W.D., McRae A.F., Visscher P.M., and Yang J. (2016). Fast set-based association analysis using summary data from GWAS identifies novel gene loci for human complex traits. Scientific Reports 6, 32894.

**GCTA Software**: Yang J, Lee SH, Goddard ME and Visscher PM. GCTA: a tool for Genome-wide Complex Trait Analysis. Am J Hum Genet. 2011 Jan 88(1): 76-82. [PubMed ID: 21167468]


### GWAS Simulation
GCTA can simulate a GWAS based on real genotype data. The phenotypes are simulated based on a set of real genotype data and a simple additive genetic model *y<sub>j</sub> = sum(w<sub>ij</sup>*u<sub>i</sub>) + e<sub>j</sub>*, where *w<sub>ij</sub> = (x<sub>ij</sub> - 2p<sub>i</sub>) / sqrt[2p<sub>i</sub>(1 - p<sub>i</sub>)]* with *x<sub>ij</sub>* being the number of reference alleles for the i-th causal variant of the j-th individual and *p<sub>i</sub>* being the frequency of the i-th causal variant, *u<sub>i</sub>* is the allelic effect of the i-th causal variant and *e<sub>j</sub>* is the residual effect generated from a normal distribution with mean of 0 and variance of *va(sum(w<sub>ij</sub>\*u<sub>i</sub>))(1 / h<sup>2</sup> - 1)*. For a case-control study, under the assumption of threshold model, cases are sampled from the individuals with disease liabilities (*y*) exceeding the threshold of normal distribution truncating the proportion of *K* (disease prevalence) and controls are sampled from the remaining individuals.
 
--simu-qt  
Simulate a quantitative trait.

--simu-cc   100   200  
Simulate a case-control study. Specify the number of cases and the number of controls, e.g. 100 cases and 200 controls. Since the simulation is based on the actual genotype data, the maximum numbers of cases and controls are restricted to be *n \* K* and *n \* (1-K)*, respectively, where n is the sample size of the genotype data.
 
--simu-causal-loci   causal.snplist  
Assign a list of SNPs as causal variants. If the effect sizes are not specified in the file, they will be generated from a standard normal distribution.
> Input file format  
> causal.snplist (columns are SNP ID and effect size)
```nohighlight
rs113645    0.025
rs185292   -0.021
...
```

--simu-hsq   0.8  
Specify the heritability (or heritability of liability), e.g. 0.8. The default value is 0.1 if this option is not specified.
 
--simu-k   0.01  
Specify the disease prevalence, e.g. 0.01. The default value is 0.1 if this option is not specified.
 
--simu-rep   100  
Number of simulation replicates.  The default value is 1 if this option is not specified.
 
#### Examples
```bash
# Simulate a quantitative trait with the heritability of 0.5 for a subset of individuals for 3 times
gcta64  --bfile test  --simu-qt  --simu-causal-loci causal.snplist  --simu-hsq 0.5 --simu-rep 3  --keep test.indi.list --out test
# Simulate 500 cases and 500 controls with the heritability of liability of 0.5 and disease prevalence of 0.1 for 3 times
gcta64  --bfile test  --simu-cc 500 500  --simu-causal-loci causal.snplist  --simu-hsq 0.5  --simu-k 0.1  --simu-rep 3  --out test
```
 
> Output file format  
> test.par (one header line; columns are the name of the causal variant, reference allele, frequency of the reference allele, and effect size).
```nohighlight
QTL             RefAllele     Frequency       Effect       
rs13626255      C             0.136           -0.0837      
rs779725        G             0.204           -0.0677      
...
```

> test.phen (no header line; columns are family ID, individual ID and the simulated phenotypes). For the simulation of a case-control study, all the individuals involved in the simulation will be outputted in the file and the phenotypes for the indivdiuals neither sampled as cases nor as controls are treated as missing, i.e. -9.
```nohighlight
011     0101     1     -9    1
012     0102     2      2    -9    
013     0103     1      1    1
...
```

### GCTA-SBLUP
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

#### References
**COJO-SBLUP method**: Robinson et al. (2017) Genetic evidence of assortative mating in humans. Nat Hum Behav, 1:0016.

**GCTA software**: Yang J, Lee SH, Goddard ME and Visscher PM (2011) GCTA: a tool for Genome-wide Complex Trait Analysis. Am J Hum Genet, 88: 76-82. [PubMed ID: 21167468]

