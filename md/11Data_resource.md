
## Data Resource

### UK Biobank GWAS results
We developed a resource-efficient tool (called [fastGWA](#fastGWA)) for mixed model association analysis, and applied it to 2,173 traits on 456,422 array-genotyped as well as 49,960 whole-exome-sequenced individuals of European ancestry in the UK Biobank (UKB). All the summary data can be downloaded using the Linux commands below. One can also query or visualize the summary data using the [online tool](http://fastgwa.info).

* GWAS summary statistics from the imputed data: 456,422 individuals of European ancestry; 8,531,416 variants (MAF > 0.01 and missingness rate < 0.1); 2,173 traits.
    * Summary table: [UKB\_impute\_v1.1.csv](./res/UKB_impute_v1.1.csv)
    * Online tool: [http://fastgwa.info/ukbimp/phenotypes](http://fastgwa.info/ukbimp/phenotypes) 
    * Linux command to download all the summary statistics (2,173 files; 454 GB in total):
```bash
mkdir ukb && cd ukb && wget http://cnsgenomics.com/software/gcta/res/UKB_impute_v1.1.list && wget -i UKB_impute_v1.1.list
```
* GWAS summary statistics from the whole-exome sequence (WES) data: 46,191 individuals of European ancestry; 152,327 variants (MAF > 0.01 and missingness rate < 0.1); 2,048 valid traits.
    * Summary table: [UKB\_WES\_v1.1.csv](./res/UKB_WES_v1.1.csv)
    * Online tool: [http://fastgwa.info/ukbwes/phenotypes](http://fastgwa.info/ukbwes/phenotypes) 
    * Linux command to download all the summary statistics (2,048 files; 8 GB in total):
```bash
mkdir wes && cd wes && wget http://cnsgenomics.com/software/gcta/res/UKB_WES_v1.1.list && wget -i UKB_WES_v1.1.list
```

#### Data format
Columns in the summary table:
```nohighlight
ID: the trait ID
Description: trait description
Data_type:  the type of phenotype (Continuous: quantitative traits; Ordered_Categorical: ordered categorical traits; Binary: binary trait)
Method: LR - Linear Regression; MLM - Mixed Linear Model. Note that fastGWA will switch to use LR for analysis if the estimated genetic variance from an MLM is not significant (p > 0.05) .
N:  sample size
Ncase: number of affected (unaffected) individuals for binary trait.
Gender_specific:  is it a gender-specific trait?
URL:  the link to download the summary statistics.
```

Association results:
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

Note: the names of the variants were kept the same as provided (the coordinates of the variants were based on GRCh37).

#### Credits and Acknowledgements
Zhili Zheng (online tool development and data analysis), Longda Jiang (data analysis), Jian Yang (overseeing). The online tool was developed based on the source code modified from Pheweb. We thank Alibaba Cloud - Australia and New Zealand for hosting the online tool.

#### Questions and Help Requests
If you have any question, please send an email to Jian Yang [jian.yang@uq.edu.au](mailto:jian.yang@uq.edu.au)

#### Citation
Jiang et al. (2019) A resource-efficient tool for mixed model association analysis of large-scale data. [bioRxiv 598110; doi:10.1101/598110](https://www.biorxiv.org/content/10.1101/598110v1).

#### Update log
* Version v1.1 (9 Aug 2019): reran all the analyses using GCTA-fastGWA v1.92.3; removed variants with MAF < 0.01 or missingness rate > 0.1 per trait, and removed binary traits with case fraction < 0.01.
* Version v1.0 (Apr 2019): first release.
