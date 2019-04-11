
## Data Resource

### UK Biobank GWAS results
We developed a resource-efficient tool (called [fastGWA](#fastGWA)) for mixed model association analysis, and applied it to 3,613 traits on 456,422 array-genotyped and 46,191 whole-exome-sequenced individuals of European ancestry in the UK Biobank (UKB).

* GWAS summary statistics from imputed/genotyped data: 456,422 individuals of European ancestry; 12,602,502 variants (MAF > 0.0001); 3613 traits.
    * Summary table: [UKB\_impute\_v1.csv](./static/UKB_impute_v1.csv)
    * Variants information: [IMPUTE\_SNP\_info\_v1.txt.gz](http://data.qld.edu.au/public/Q1031/IMPUTE_SNP_info_v1.txt.gz)) 
    * Summary statistics: [http://data.qld.edu.au/public/Q1031/UKB\_impute\_v1/](http://data.qld.edu.au/public/Q1031/UKB_impute_v1/)
* GWAS summary statistics from whole-exome sequence (WES) data: 46,191 individuals of European ancestry; 3,264,503 variants; 3610 valid traits.
    * Summary table: [UKB\_WES\_v1.csv](./static/UKB_WES_v1.csv)
    * Variants information: [WES\_SNP\_info\_v1.txt.gz](http://data.qld.edu.au/public/Q1031/WES_SNP_info_v1.txt.gz)
    * Summary statistics: [http://data.qld.edu.au/public/Q1031/UKB\_WES\_v1/](http://data.qld.edu.au/public/Q1031/UKB_WES_v1/)

Note:   
1) The names of the variants were kept the same as provided. Thus, the coordinates of the variants in the UKB imputed data were based on GRCh37, whereas those in the UKB WES data were based on GRCh38.  
2) You can download all the summary data files by the following Linux commands.
```bash
wget http://data.qld.edu.au/public/Q1031/UKB_impute_v1.list && wget -i UKB_impute_v1.list
wget http://data.qld.edu.au/public/Q1031/UKB_WES_v1.list && wget -i UKB_WES_v1.list
```

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
N:    sample size
AF1:  the allele frequency of A1
beta: SNP effect
se:   standard error
p:    p value
```

#### Citation
Jiang L, Zheng Z, Qi T, Kemper KE, Wray NR, Visscher PM, Yang J. (2019) A resource-efficient tool for mixed model association analysis of large-scale data. bioRxiv.

