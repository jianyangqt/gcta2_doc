
## Population genetics

### GCTA-Fst: calculating Fst using GWAS data

Calculating Fst following the method described in Weir (1996).

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

#### References

**Fst method**: Weir, B.S. 1996. Genetic Data Analysis II: Methods for Discrete Population Genetic Data. Sinauer Associates, Inc. Sunderland, Massachusetts.

**GCTA software**: Yang J, Lee SH, Goddard ME and Visscher PM. (2011) GCTA: a tool for Genome-wide Complex Trait Analysis. Am J Hum Genet, 88: 76-82. [PubMed ID: 21167468]


### GCTA-PCA: Principal component analysis

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

### PC loading and projection

--pc-loading  pca
Generate the SNP loading according the the GCTA-PCA results

--project-loading snp_loading 20  
Project the genotype into the SNP loading, and generate 20 principal component

> Note: 
> * The SNPs mismatch will obsolete in the PC loadings, so it might be biased if lots of SNPs missing in the target genotype.
> * The default mode is using the MAF in reference genotype to caculate the projection, which is coincided with EIGENSOFT projection.

**Example**  
SRC: source genotype, the source genotype to generate PC  
TAR: target genotype to project the PC to SRC
```bash
# make GRM
gcta64 --bfile SRC --maf 0.01 --autosome --make-grm --out SRC
# PCA analysis
gcta64 --grm SRC --pca 20 --out SRC_pca20

# use the pca generated above to produce the SNP loading
gcta64 --bfile SRC --pc-loading SRC_pca20 --out SRC_snp_loading

# project the TAR to the SNP loading. The number 20 is the PC want to project.
# It is the same if do chromosome by chromosome caculation add them up after each caculation 
gcta64 --bfile TAR --project-loading SRC_snp_loading 20 --out TAR_pca20
```
