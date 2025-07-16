
## Data Resource

### UK Biobank GWAS results

We developed two resource-efficient tools (called [fastGWA](#fastGWA) and [fastGWA-GLMM](#fastGWA-GLMM)) for mixed model-based association analysis in large-scale data (Jiang et al. 2019 Nat Genet and Jiang et al. 2021 Preprint). 1. We first applied fastGWA to 2,173 traits on 456,422 array-genotyped as well as 49,960 whole-exome-sequenced individuals of European ancestry in the UK Biobank (UKB) (Jiang et al. 2019 Nat Genet). 2. We then applied fastGWA-GLMM to 2,989 binary traits on 456,348 array-genotyped individuals of European ancestry in the UKB (Jiang et al. 2021 Preprint). See below for detailed instructions to query and download the data. One can also query or visualize the GWAS summary data using the online tool for [fastGWA](http://fastgwa.info/ukbimp) and [fastGWA-GLMM](http://fastgwa.info/ukbimpbin).



* fastGWA summary statistics using the imputed data (Jiang et al. 2019 Nat Genet): 456,422 individuals of European ancestry; 8,531,416 variants (MAF > 0.01 and missingness rate < 0.1); 2,173 traits.
    * Summary table: [UKB\_impute\_v1.1.csv](./res/UKB_impute_v1.1.csv)
    * Online tool: [http://fastgwa.info/ukbimp/phenotypes](http://fastgwa.info/ukbimp/phenotypes) 
    * Linux command to download all the summary statistics (2,173 files; 454 GB in total):
```bash
mkdir ukb && cd ukb && wget http://cnsgenomics.com/software/gcta/res/UKB_impute_v1.1.list && wget -i UKB_impute_v1.1.list
```
* fastGWA summary statistics using the whole-exome sequence (WES) data (Jiang et al. 2019 Nat Genet): 46,191 individuals of European ancestry; 152,327 variants (MAF > 0.01 and missingness rate < 0.1); 2,048 valid traits.
    * Summary table: [UKB\_WES\_v1.1.csv](./res/UKB_WES_v1.1.csv)
    * Online tool: [http://fastgwa.info/ukbwes/phenotypes](http://fastgwa.info/ukbwes/phenotypes) 
    * Linux command to download all the summary statistics (2,048 files; 8 GB in total):
```bash
mkdir wes && cd wes && wget http://cnsgenomics.com/software/gcta/res/UKB_WES_v1.1.list && wget -i UKB_WES_v1.1.list
```
* fastGWA-GLMM summary statistics using the imputed data (Jiang et al. 2021 Preprint): 456,422 individuals of European ancestry; 11,842,647 variants (MAF > 0.0001 and missingness rate < 0.1); 2,989 binary traits.
    * Summary table: [UKB\_binary\_v1.11.csv](./res/UKB_binary_v1.11.csv)
    * Online tool: [http://fastgwa.info/ukbimpbin/phenotypes](http://fastgwa.info/ukbimpbin/phenotypes) 
    * Linux command to download all the summary statistics (2,989 files; 1.2 TB in total):
```bash
mkdir ukb_binary && cd ukb_binary && wget http://cnsgenomics.com/software/gcta/res/UKB_binary_v1.11.list && wget -i UKB_binary_v1.11.list
```


#### Data format
Columns in the summary table for results from fastGWA:
```nohighlight
ID: the trait ID.
Description: trait description.
Data_type:  the type of phenotype (Continuous: quantitative trait; Ordered_Categorical: ordered categorical trait; Binary: binary trait)
Method: LR - Linear Regression; MLM - Mixed Linear Model. Note that the program will switch to use LR for analysis if the estimated genetic variance from an fastGWA is not significant (p > 0.05).
N: sample size.
Ncase: number of affected individuals for a binary disease trait.
Gender_specific:  to indicate if it is a gender-specific trait.
URL: the link to download the summary statistics.
```

Columns in the summary table for results from fastGWA-GLMM:
```nohighlight
ID: the trait ID.
Description: trait description.
N_case: number of affected individuals for a binary disease trait.
N_control: number of unaffected individuals for a binary disease trait.
ratio: case-control ratio.
URL:  the link to download the summary statistics.
```

Format of the summary statistics from fastGWA:
```nohightlight
CHR:  chromosome
SNP:  SNP ID
POS:  SNP position
A1:   effect allele
A2:   the other allele
N:    per allele sample size
AF1:  the allele frequency of A1
BETA: SNP effect
SE:   standard error
P:    p value
```

Format of the summary statistics from fastGWA-GLMM:
```nohightlight
CHR:  chromosome
SNP:  SNP ID
POS:  SNP position
A1:   effect allele
A2:   the other allele
N:    per allele sample size
AF1:  the allele frequency of A1
T:    GLMM score statistic
SE_T: standard error of the score statistic
P_noSPA: raw p-value
BETA: SNP effect or log(odds ratio)
SE:   standard error for the estimated effect size after the SPA correction
P:    p-value after the SPA correction
CONVERGE: to indicate whether the SPA correction is converged for that variant
```


Note: the names of the variants were kept the same as provided (the coordinates of the variants were based on GRCh37).

#### Credits and Acknowledgements
Zhili Zheng (online tool development and data analysis), Longda Jiang (data analysis), Jian Yang (overseeing). The online tool was developed based on the source code modified from Pheweb. We thank Alibaba Cloud - Australia and New Zealand for hosting the online tool.

#### Questions and Help Requests
If you have any question, please send an email to Jian Yang [jian.yang@westlake.edu.cn](mailto:jian.yang@westlake.edu.cn)

#### Citation
Jiang L, Zheng Z, Qi T, Kemper KE, Wray NR, Visscher PM, Yang J (2019) A resource-efficient tool for mixed model association analysis of large-scale data. Nature Genetics, 51: 1749–1755.

Jiang L, Zheng Z, Yang J (2021) FastGWA-GLMM: a generalized linear mixed model association tool for biobank-scale data, February 2021, PREPRINT (Version 1) available at Research Square https://doi.org/10.21203/rs.3.rs-128758/v1

#### Update log
* Version v1.11 (14 May 2021): added GWAS summary statistics from GCTA-fastGWA-GLMM v1.93.3 for 2,989 binary traits from the UKB.
* Version v1.1 (9 Aug 2019): reran all the analyses using GCTA-fastGWA v1.92.3; removed variants with MAF < 0.01 or missingness rate > 0.1 per trait, and removed binary traits with case fraction < 0.01.
* Version v1.0 (Apr 2019): first release.
