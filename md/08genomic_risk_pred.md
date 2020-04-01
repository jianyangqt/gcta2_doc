## Genomic risk prediction

### BLUP
BLUP (best linear unbiased prediction) is a method frequently used to estimate the random effects in mixed linear models (MLMs). In GCTA, we have implemented options to estimate the BLUP effects of individuals’ genetic values, and then transform them to the BLUP effects of all SNPs. The BLUP SNP effects can be used to predict the PRS (polygenic risk scores) in an independent sample. We have also implemented the cross validated BLUP (cvBLUP) approach to perform a leave-one-out BLUP analysis in one run ([Mefford et al. 2019 bioRxiv](https://www.biorxiv.org/content/10.1101/517821v1.full)). The BLUP and cvBLUP functions share multiple flags with --reml.

--reml-pred-rand  
To estimate random effects in a MLM by BLUP. It can be used to predict the genetic value (called "breeding value" in animal genetics) of each individual attributed to the aggregative effect of all SNPs used to compute the GRM. The genetic values of all the individuals will be saved in a plain text file \*.indi.blp.

> Output file format
>> test.indi.blp (no headers). The columns are family ID, individual ID, an intermediate variable, the genetic value, another intermediate variable and the residual effect. If there are multiple GRMs fitted in the model, each GRM will add two additional columns, i.e., an intermediate variable (= **Py**, see Yang et al. 2011 AJHG for the definitions of matrix **P** and vector **y**) and a genetic value, before the last two columns.

```nohighlight
01       0101    -0.012    -0.014   -0.010    -0.035
02       0203    0.021     0.031    -0.027    -0.031
03       0305    0.097     0.102    -0.026    -0.041
...
```

--blup-snp test.indi.blp  
To compute the BLUP solutions for the SNPs (you will need to specify the option --bfile to read the genotype data). This option takes the output from --reml-pred-rand as the input (\*.indi.blp file) and transforms the BLUP solutions for individuals to those for SNPs. The SNP BLUP effects can then be used to predict the genetic values (so called PRS) of individuals in an independent sample by the PLINK --score option. Note that for the ease of using the BLUP SNP effects in PLINK --score, the BLUP SNPs effects are scaled by sqrt[2p(1-p)]\(see pages 77 and 78 of Yang et al. 2011 AJHG for details).

> Output file format
>> test.snp.blp (no heades). The columns are SNP ID, reference allele and BLUP of SNP effect. If there are multiple GRMs fitted in the model, each GRM will add an additional column to the file. The last column is for the residual effect.

```nohighlight
rs103645   A     0.00312    0.00451
rs175292   G    -0.00021    0.00139
...
```

--cvblup  
cvBLUP is a computationally efficient approach that uses a mixed linear model (MLM) to generate polygenic risk scores (PRSs), which are equivalent to those from a leave-one-individual-out BLUP analysis, for all the individuals in a sample (see [Mefford et al. 2020 J Comput. Biol.](https://www.liebertpub.com/doi/full/10.1089/cmb.2019.0325) for details of the method).

> Output file format
>> test.indi.cvblp (no headers). The columns are family ID, individual ID, an intermediate variable, the genetic effect, another intermediate variable and the residual effect.

```nohighlight
01       0101    -0.012    -0.014   -0.010    -0.035
02       0203    0.021     0.031    -0.027    -0.031
03       0305    0.097     0.102    -0.026    -0.041
...
```

#### Examples
```bash
# To obtain BLUP solutions for the genetic values of individuals
gcta64 --reml --grm test --pheno test.phen --reml-pred-rand –qcovar test_10PCs.txt  --out test
# To obtain BLUP solutions for the SNP effects
gcta64 --bfile geno_to_predict --blup-snp test.indi.blp --out test
# To compute the polygenic risk score (PRS) in an independent sample
plink --bfile geno_to_predict --score test.snp.blp 1 2 3

# To obtain cvBLUP solutions for the genetic values of individuals
gcta64 --reml --grm test --pheno test.phen --cvblup –qcovar test_10PCs.txt --out test
# To obtain cvBLUP solutions for the SNP effects
gcta64 --bfile geno --blup-snp test.indi.cvblp --out test
# To compute the polygenic risk score (PRS)
plink --bfile geno --score test.snp.blp 1 2 3
```

#### Citation
BLUP method or GCTA software: Yang J, Lee SH, Goddard ME and Visscher PM. GCTA: a tool for Genome-wide Complex Trait Analysis. Am J Hum Genet. 2011 Jan 88(1): 76-82. \[PubMed ID: [21167468](https://www.ncbi.nlm.nih.gov/pubmed/21167468)]

cvBLUP method: Mefford JA, Park D, Zheng Z, et al. Efficient Estimation and Applications of Cross-Validated Genetic Predictions to Polygenic Risk Scores and Linear Mixed Models. Journal of Computational Biology. 2020 Feb online; doi: https://doi.org/10.1089/cmb.2019.0325



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

