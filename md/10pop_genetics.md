
## Population genetics

### Fst

**GCTA-Fst: calculating Fst values of SNPs using GWAS data**

It follows the Fst method described in Weir (1996).

> Example
```bash
gcta64 --bfile test --fst --sub-popu subpopu.txt --out test
```

> Format of input file ("subpopu.txt")
```nohighlight
1	11	Popu1
2	21	Popu1
3	31	Popu2
4	41	Popu2
5	51	Popu1
6	61	Popu2
...
```

> Output  
> Results are saved in *.fst file.
```nohighlight
Chr	SNP	      bp	 refA freq_Popu1(n=1000) freq_Popu2(n=2925)	Fst	
1	rs4475691	836671	T	0.208561	0.193984	0.000508832	
1	rs28705211	890368	C	0.287427	0.274928	0.000295543	
1	rs9777703	918699	C	0.0265765	0.0313871	0.000300492	
...
```


#### Citations

**Fst method**: Weir, B.S. 1996. Genetic Data Analysis II: Methods for Discrete Population Genetic Data. Sinauer Associates, Inc. Sunderland, Massachusetts.

**GCTA software**: Yang J, Lee SH, Goddard ME and Visscher PM. (2011) GCTA: a tool for Genome-wide Complex Trait Analysis. Am J Hum Genet, 88: 76-82. [PubMed ID: 21167468]


### PCA

**GCTA-PCA: Principal component analysis**

--pca 20  
Input the GRM and output the first n (n = 20, by default) eigenvectors (saved as *.eigenvec, plain text file) and all the eigenvalues (saved as *.eigenval, plain text file), which are equivalent to those calculated by the program EIGENSTRAT. The only purpose of this option is to calculate the first m eigenvectors, and subsequently include them as covariates in the model when estimating the variance explained by all the SNPs (see below for the option of estimating the variance explained by genome-wide SNPs). Please find the EIGENSTRAT software if you need more sophisticated principal component analysis of the population structure. 
> Output file format  
> test.eigenval (no header line; the first m eigenvalues)
```nohighlight
20.436  
7.1293  
6.7267  
......
```

> test.eigenvec (no header line; the first m eigenvectors; columns are family ID, individual ID and the first m eigenvectors)  
```nohighlight
011      0101       0.00466824      -0.000947       0.00467529      -0.00923534  
012      0102       0.00139304      -0.00686406     -0.0129945      0.00681755  
013      0103       0.00457615      -0.00287646     0.00420995      -0.0169046  
......
```

Examples
```bash
# Input the GRM file and output the first 20 eigenvectors for a subset of individuals
gcta64  --grm test --keep test.indi.list  --pca 20  --out test
```

#### Citations

**PCA method**: Price AL, Patterson NJ, Plenge RM, Weinblatt ME, Shadick NA and Reich D (2006) Principal components analysis corrects for stratification in genome-wide association studies. Nat Genet, 38: 904-909. [PubMed ID: 16862161]

**GCTA software**: Yang J, Lee SH, Goddard ME and Visscher PM. (2011) GCTA: a tool for Genome-wide Complex Trait Analysis. Am J Hum Genet, 88: 76-82. [PubMed ID: 21167468]

### PC loading and projection

--pc-loading ref_pca  
To generate the loadings of each SNP on the PCs computed from a reference sample (or a subset of the whole sample)

--project-loading ref\_snp\_loading 20  
To project the SNP genotypes of a target sample onto a reference sample (or a subset of the target sample) via the precomputed SNP loadings (see the option above), and generate PCs in the target sample (20 PCs by default).

> Note: 
> * Only the SNPs in common between the target and reference samples will be used in the PC projection analysis.
> * By default, PC loadings of the SNPs will be computed based on allele frequencies estimated in the reference sample (similar to the PC projection approach implemented in EIGENSOFT), which can be changed by the option --update-freq with allele frequencies computed from the target sample or a meta-analysis of the target and reference samples.

**Example**  
REF: SNP genotype data of the reference sample; 
TAR: SNP genotype data of the target sample;
```bash
# To make a GRM
gcta64 --bfile REF --maf 0.01 --autosome --make-grm --out REF
# PCA analysis
gcta64 --grm REF --pca 20 --out REF_pca20

# To use the PCs generated above to produce PC loadings of each SNP
gcta64 --bfile REF --pc-loading REF_pca20 --out REF_snp_loading

# To compute the PCs of the target sample using the PC loading generated above
# Note that the analysis can be performed with one chromosome at a time 
gcta64 --bfile TAR --project-loading REF_snp_loading 20 --out TAR_pca20
```

#### Citation
Jiang et al. (2019) A resource-efficient tool for mixed model association analysis of large-scale data. Nature Genet. 51, 1749–1755 (2019) [doi:10.1038/s41588-019-0530-8](https://www.nature.com/articles/s41588-019-0530-8).
