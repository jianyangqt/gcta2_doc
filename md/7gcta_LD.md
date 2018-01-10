
## LD

### Computing LD scores

**GCTA-LDS: calculating LD score for each SNP**

LD score is defined as the sum of LD *r<sup>2</sup>* between a variant and all the variants in a region.

Example - calculating LD score to stratify SNPs
```bash
gcta64 --bfile test --ld-score --ld-wind 1000 --ld-rsq-cutoff 0.01 --out test
```

--ld-wind 10000  
The default value is L = 10000 (in Kb unit), i.e. L = 10Mb. The genome is chopped into segments with length of L for LD calculation. Two adjacent segments are overlapped. The size of the overlap is L/2.

--ld-rsq-cutoff 0  
The default value is 0. LD r<sup>2</sup> smaller than this value will be ignored in the calculated. If the threshold is > 0, the LD score estimate is biased since it's alway positive. The LD score generated from this option can be used for stratifying SNPs (see GREML-LDMS).

> Output
```nohighlight
SNP chr bp MAF mean_rsq snp_num max_rsq ldscore
rs12260013 10 66326 0.0709329 0.0475853 2211 0.807478 106.211
...
```
mean_rsq: mean LD r<sup>2</sup> between the target SNP and all other SNPs in the window.  
snp_num: number of SNPs used in the calculation  
max_rsq: maximum LD r<sup>2</sup> between the target SNP and its best tagging SNP in the window.  
ldscore: LD score

LD score is calculated as  
```r
1 + mean_rsq * snp_num
```

Example - calculating LD score for LDSC regression
```bash
gcta64 --bfile test --ld-score --ld-wind 1000 --ld-score-adj --out test
```

--ld-score-adj  
LD r<sup>2</sup> is alway positive which is not an unbiased estimate of squared correlation (rho2). The adjustment is r<sup>2</sup>adj = r<sup>2</sup> - [(1 - r<sup>2</sup>) / (n -2)], where n is the sample size (Bulik-Sullivan et al. 2015 Nat Genet). 
The output from this analysis can be used for LDSC regression analysis. We do not recommend using the --ld-rsq-cutoff option in this analysis. Otherwise, the LD score estimate is biased.

--ld-score-multi  
Creating LD score of each SNP against multiple SNP set. This option can be used to perform multi-component LDSC regression analysis following Finucane et al. (2015 Nat Genet). Note that the --ld-score-adj option also applies to this analysis.

> Output - the same as above.

Example - calculating LD scores for multi-component LDSC regression
```bash
gcta64 --bfile test --ld-score-multi test_multi_snplist.txt --ld-wind 1000 --out test
```

> Note that this is an analysis of calculating the LD score for each SNP against multiple SNP sets, e.g. the LD score of each SNP against all SNPs in exons and that against all SNPs in introns.

> Input format 
```nohighlight
test_multi_snplist.txt
test_snp_set1.snplist
test_snp_set2.snplist
...
```

> Output - an example of calculating LD score of SNP against two SNP sets
```nohighlight
SNP chr bp MAF mean_rsq_1 snp_num_1 max_rsq1 ldscore1 mean_rsq_2 snp_num_2 max_rsq2 ldscore2
rs4475691 1 836671 0.197698 0.000867814 499 0.216874 0.000308932 500 0.0022564
rs28705211 1 890368 0.278112 0.000911328 499 0.216874 0.000237098 500 0.00254858
rs9777703 1 918699 0.0301614 0.00240581 499 0.854464 0.000222185 500 0.00222427
...
```


#### Citations

**LD score regression analysis**: Bulik-Sullivan BK, Loh PR, Finucane HK, Ripke S, Yang J, Schizophrenia Working Group of the Psychiatric Genomics Consortium, Patterson N, Daly MJ, Price AL, Neale BM (2015) LD Score regression distinguishes confounding from polygenicity in genome-wide association studies. Nat Genet, 47: 291-295.

**GCTA software**: Yang J, Lee SH, Goddard ME and Visscher PM. (2011) GCTA: a tool for Genome-wide Complex Trait Analysis. Am J Hum Genet, 88: 76-82. [PubMed ID: 21167468]


### Searching for LD friends 

**GCTA-LDF: finding LD friends for each of the SNPs given on the list**
	
For each target SNP, GCTA uses the simple regression approach to search for SNPs that are in significant LD with the target SNP.

--ld ld.snplist  
Specify a list of SNPs.

--ld-wind 5000  
Search for SNPs in LD with a target SNP within d Kb (e.g. 5000 Kb) region in either direction by simple regression test.

--ld-sig 0.05  
Threshold p-value for regression test, e.g. 0.05.

Example
```bash
gcta64  --bfile test  --ld ld.snplist  --ld-wind 5000  --ld-sig 0.05  --out test
```

> Output files

1) test.rsq.ld, summary of LD structure with each row corresponding to each target SNP. The columns are target SNP
```nohighlight
length of LD block
two flanking SNPs of the LD block
total number of SNPs within the LD block
mean r2 (r squared)
median r2
maximum r2
SNP in highest LD with the target SNP
```
2) test.r.ld, the correlations (r) between the target SNP and all the SNPs in the LD block.  
3) test.snp.ld, the names of all the SNPs in the LD with the target SNP.  
> Note: LD block is defined as a region where SNPs outside this region are not in significant LD with the target SNP. According to this definition, the length of LD block depends on user-specified window size and significance level.


#### Citations

**Method to search for LD friends**: Yang et al. (2011) Genomic inflation factors under polygenic inheritance. Eur J Hum Genet. 19(7): 807-812. [Pubmed ID: 21407268]

**GCTA software**: Yang J, Lee SH, Goddard ME and Visscher PM. GCTA: a tool for Genome-wide Complex Trait Analysis. Am J Hum Genet. 2011 Jan 88(1): 76-82. [PubMed ID: 21167468]


