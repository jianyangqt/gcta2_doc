
## Data Resource

### UK Biobank GWAS results
We developed a resource-efficient tool (called [fastGWA](#fastGWA)) for mixed model association analysis, applied it to 3,613 traits on 456,422 array-genotyped and 46,191 whole-exome-sequenced individuals of European ancestry in the UK Biobank (UKB).

* GWAS summary statistics from imputed/genotyped data: 456,422 individuals of European ancestry; 12,602,502 variants (MAF > 0.0001); 3613 traits (summary table [UKB\_impute.csv](./static/UKB_impute.csv), variants information [IMPUTE\_SNP\_info.txt.gz](http://data.qld.edu.au/public/Q1031/IMPUTE_SNP_info.txt.gz)).
* GWAS summary statistics from whole-exome sequence (WES) data: 46,191 individuals of European ancestry; 3,264,503 variants; 3610 valid traits (summary table [UKB\_WES.csv](./static/UKB_WES.csv), variants information [WES\_SNP\_info.txt.gz](http://data.qld.edu.au/public/Q1031/WES_SNP_info.txt.gz)).

Note:  the names of the variants were kept the same as provided. Thus, the coordinates of the UKB imputed data were based on GRCh37, whereas the UKB WES data were based on GRCh38.

Columns in the summary table:
```nohighlight
ID: the trait ID in the UKB
Description: trait description
Data_type:  the type of phenotype (Continuous: quantitative traits; Ordered_Categorical: ordered categorical traits; Binary: binary trait)
Method: LR - Linear Regression; MLM - Mixed Linear Model. Note that the fastGWA program will switch to use LR for analysis if the estimated genetic variance from an MLM is not significant (p > 0.05) .
N:  sample size
N_Case (N_control): number of affected (unaffected) individuals for binary trait.
Gener_specific:  Is the trait limited to one gender?
URL:  the link to download the association summary statistics.
```

Variants information:
```nohighlight
CHR:  chromosome
SNP:  SNP ID
POS:  SNP position
A1:   effect allele
A2:   the other allele
```

Association Results:
```nohightlight
SNP:  SNP ID
AF1:  the allele frequency of A1
beta: SNP effect
se:   standard error
p:    p value
```

#### Citation
Jiang L, Zheng Z, Qi T, Kemper KE, Wray NR, Visscher PM, Yang J. (2019) A resource-efficient tool for mixed model association analysis of large-scale data. bioRxiv.


