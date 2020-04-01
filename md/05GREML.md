
## GREML

### Tutorial

If you have used PLINK before, you will find it easy to use GCTA. In this tutorial, all the options used are not detailed. Please refer to the documentation of GCTA for details of the options and formats of the input or output files.

#### GCTA-GRM: calculating the genetic relationship matrix (GRM) from all the autosomal SNPs
Suppose you have a GWAS data set in PLINK binary PED format, e.g. test.bed, test.bim and test.fam. You can type this command to calculate the genetic relationships between pairwise individuals from all the autosomal SNPs
```bash
gcta64 --bfile test --autosome --maf 0.01 --make-grm --out test --thread-num 10
```
The genetic relationship matrix will be saved in the files test.grm.bin, test.grm.N.bin and test.grm.id .

For datasets with an extremely large number of SNPs and large sample size (e.g. 1000G imputed data, you can use the following commands:

```bash
gcta64 --bfile test --chr 1 --maf 0.01 --make-grm --out test_chr1 --thread-num 10
gcta64 --bfile test --chr 2 --maf 0.01 --make-grm --out test_chr2 --thread-num 10
...
gcta64 --bfile test --chr 22 --maf 0.01 --make-grm --out test_chr22 --thread-num 10
```
which calculate the GRM for each autosome and then merge the 22 GRMs by the following command:
```bash
gcta64 --mgrm grm_chrs.txt --make-grm --out test
```
You can use this command to remove cryptic relatedness

```bash
gcta64 --grm test --grm-cutoff 0.025 --make-grm --out test_rm025
```

which creates a new GRM of "unrelated" individuals. Please be aware that the cutoff value 0.025 is quite arbitrary.

#### GCTA-GREML analysis: estimating the variance explained by the SNPs
```bash
gcta64 --grm test --pheno test.phen --reml --out test --thread-num 10
```
The results will be saved in the file test.hsq.  

You can also include the first 4 or 10 eigenvectos from principal component analysis (PCA) as covariates by the command
```bash
gcta64 --grm test --pheno test.phen --reml --qcovar test_10PCs.txt --out test --thread-num 10
```
You can also estimate the variance explained by the SNPs on each chromosome by fitting one chromosome at a time
```bash
gcta64 --grm test_chr1 --pheno test.phen --reml --out test_chr1 --thread-num 10
gcta64 --grm test_chr2 --pheno test.phen --reml --out test_chr2 --thread-num 10
......
gcta64 --grm test_chr22 --pheno test.phen --reml --out test_chr22 --thread-num 10
```
or fitting all the 22 autosomes simultaneously by
```bash
gcta64 --mgrm grm_chrs.txt --pheno test.phen --reml --out test_all_chrs --thread-num 10
```
You are also allowed to include the first 4 or 10 eigenvectors from PCA as covariates in any of these analyses.

#### GCTA-GREML analysis for a case-control study
For a case-control study, the phenotypic values of cases and controls should be specified as 1 and 0, respectively. Suppose you have prepared a phenotype file test_cc.phen. You can type the following command to estimate the variance explained by all the autosomal SNPs on the observed 0-1 scale and transform the estimate to that on the underlying liability scale (assuming the disease prevalence is 0.01 in this example)
```bash
gcta64 --grm test --pheno test_cc.phen --reml --prevalence 0.01 --out test --thread-num 10
```

### Making a GRM

**GCTA-GRM: estimating genetic relatedness from SNPs**

--make-grm  
or  
--make-grm-bin  
Estimate the genetic relationship matrix (GRM) between pairs of individuals from a set of SNPs and save the lower triangle elements of the GRM to binary files, e.g. test.grm.bin, test.grm.N.bin, test.grm.id.

> Output file 
>> test.grm.bin (it is a binary file which contains the lower triangle elements of the GRM).  
>> test.grm.N.bin (it is a binary file which contains the number of SNPs used to calculate the GRM).  
>> test.grm.id (no header line; columns are family ID and individual ID, see above).  
>> You can not open test.grm.bin or test.grm.N.bin by a text editor but you can use the following R script to read them in R)  
```r
# R script to read the GRM binary file
ReadGRMBin=function(prefix, AllN=F, size=4){
  sum_i=function(i){
    return(sum(1:i))
  }
  BinFileName=paste(prefix,".grm.bin",sep="")
  NFileName=paste(prefix,".grm.N.bin",sep="")
  IDFileName=paste(prefix,".grm.id",sep="")
  id = read.table(IDFileName)
  n=dim(id)[1]
  BinFile=file(BinFileName, "rb");
  grm=readBin(BinFile, n=n*(n+1)/2, what=numeric(0), size=size)
  NFile=file(NFileName, "rb");
  if(AllN==T){
    N=readBin(NFile, n=n*(n+1)/2, what=numeric(0), size=size)
  }
  else N=readBin(NFile, n=1, what=numeric(0), size=size)
  i=sapply(1:n, sum_i)
  return(list(diag=grm[i], off=grm[-i], id=id, N=N))
}
```

**Note**: --make-grm has been rewritten with orders of magnitude improvement in speed and memory usage. Currently, It can only used in combination with a limited number of other flags, i.e., --keep, --remove, --chr, --autosome-num, --autosome, --extract, --exclude, --maf, --max-maf, --thread-num, --update-ref-allele, --update-sex, --update-freq. You can use --make-grm-part to reduce the memory usage further.

Make GRM function can combine with --mbfile to calculate GRMs in multiple PLINK files without merge them together. 

--mbfile chrs.txt  
If the genotype data is very large, the data is often saved in separate PLINK files (e.g. one for each chromosome). Use --mbfile to specify multiple PLINK files. The input is a text file with each row representing a PLINK binary file (without file name suffix).
> Input file format
```nohighlight
data_chr1
data_chr2
…
```

**Note**: All these files shall have same sample size and order, the program will prompt an error if not. 

--make-grm-part m i  
Partition the GRM into m parts (by row), and compute the i-th part in the current run. 

**Note**: This option is designed to compute the GRM in a very large sample (e.g. the UK Biobank data). The memory usage of each run is the total memory required divided by m. Thus partitioning a large number of parts can reduce the memory usage significantly. The total memory required is approximately [n * (n + 1) / 2 * 12] / 1024<sup>3</sup> GB + 0.5GB, where n is the sample size. As some computer clusters limit the virtual memory, allocating 1 to 2GB more memory to each job will be safer. In our computation of the GRM in the UKB data, we partitioned the whole data set (n = 456,426) into 250 parts and allocated 6700MB memory to each job.

> Example:
```bash
# Partition the GRM into 3 parts
gcta64 --bfile test --make-grm-part 3 1 --thread-num 5 --out test
gcta64 --bfile test --make-grm-part 3 2 --thread-num 5 --out test
gcta64 --bfile test --make-grm-part 3 3 --thread-num 5 --out test
# Merge all the parts together (Linux, Mac)
cat test.part_3_*.grm.id > test.grm.id
cat test.part_3_*.grm.bin > test.grm.bin
cat test.part_3_*.grm.N.bin > test.grm.N.bin
# Windows alternative
copy /b test.part_3_*.grm.id test.grm.id
copy /b test.part_3_*.grm.bin test.grm.bin
copy /b test.part_3_*.grm.N.bin test.grm.N.bin
```
 
--make-grm-alg 0    
The default value is 0, and the GRM is calculated using the equation *sum{[(x<sub>ij</sub> - 2p<sub>i</sub>)\*(x<sub>ik</sub> - 2p<sub>i</sub>)] / [2p<sub>i</sub>(1-p<sub>i</sub>)]}* as described in Yang et al. 2010 Nat Genet. If the value = 1, the GRM will be calculated using the equation *sum[(x<sub>ij</sub> - 2p<sub>i</sub>)*(x<sub>ik</sub> - 2p<sub>i</sub>)] / sum[2p<sub>i</sub>(1-p<sub>i</sub>)]*. 

--make-grm-gz  
Estimate the GRM, save the lower triangle elements to a compressed text file (e.g. test.grm.gz) and save the IDs in a plain text file (e.g. test.grm.id).
> Output file format  
> test.grm.gz (no header line; columns are indices of pairs of individuals (row numbers of the test.grm.id), number of non-missing SNPs and the estimate of genetic relatedness)
```nohighlight
1    1    1000    1.0021
2    1    998     0.0231
2    2    999     0.9998
3    1    1000    -0.0031
...
```

> test.grm.id (no header line; columns are family ID and individual ID)
```nohighlight
011      0101
012      0102
013      0103
...
```

--make-grm-xchr  
Estimate the GRM from SNPs on the X-chromosome. The GRM will be saved in the same binary format as above (\*.grm.bin, \*.grm.N.bin and \*.grm.id). Due to the speciality of the GRM for the X-chromosome, it is not recommended to manipulate the matrix by --grm-cutoff or --grm-adj, or merge it with the GRMs for autosomes (see below for the options of manipulating the GRM). 

Note 1: this flag has been re-implemented in GCTA 1.91.4, it has same performance and memory consumption as --make-grm.

Note 2: the function treats X chr as non-pseudoautosomal region (nPAR) with genotype coding for male as 0, 2. For pseudoautosomal region (PAR), we can alter the chromosome number in bim file to autosome and use --make-grm to run. Don't put nPAR and PAR together as X chr, GCTA will give weird results.

--make-grm-xchr-part m i  
Partition the GRM of X chromosome into m parts (by row), and compute the i-th part in the current run.

See the document of [--make-grm-part](#MakingaGRM)


--make-grm-xchr-gz  
Same as --make-grm-xchr but the GRM will be in compressed text files (see --make-grm-gz for the format of the output files).

--make-grm-inbred or --make-grm-inbred-gz  
Make a GRM for an inbred population such as inbred mice or inbred crops.

--ibc  
Estimate the inbreeding coefficient from the SNPs by 3 different methods.
> Output file format  
> test.ibc (one header line; columns are family ID, individual ID, number of nonmissing SNPs, estimator 1, estimator 2 and estimator 3)
```nohighlight
FID      IID        NOMISS       Fhat1       Fhat2         Fhat3
011      0101       999          0.00210     0.00198       0.00229
012      0102       1000         -0.0033     -0.0029       -0.0031
013      0103       988          0.00120     0.00118       0.00134
```
See [Yang et al. 2011 AJHG](http://www.cell.com/ajhg/abstract/S0002-9297(10)00598-7) for the definitions of Fhat1, Fhat2 and Fhat3.

> Examples
```bash
# Estimate the GRM from all the autosomal SNPs
gcta64  --bfile test  --autosome  --make-grm  --out test

# Estimate the GRM from the SNPs on the X-chromosome
gcta64  --bfile test  --make-grm-xchr  --out test_xchr

# Estimate the GRM from the SNPs on chromosome 1 with MAF from 0.1 to 0.4
gcta64  --bfile test  --chr 1  --maf 0.1  --max-maf 0.4  --make-grm  --out test

# Estimate the GRM using a subset of individuals and a subset of autosomal SNPs with MAF < 0.01
gcta64  --bfile test  --keep test.indi.list  --extract test.snp.list  --autosome  --maf 0.01 --make-grm  --out test

# Estimate the GRM from the imputed dosage scores for the SNPs with MAF > 0.01 and imputation R2 > 0.3
gcta64  --dosage-mach  test.mldose.gz  test.mlinfo.gz  --imput-rsq  0.3  --maf 0.01  --make-grm --out test

# Estimate the GRM from the imputed dosage scores for a subset of individuals and a subset of SNPs
gcta64  --dosage-mach  test.mldose.gz  test.mlinfo.gz  --keep test.indi.list  --extract test.snp.list  --make-grm --out test

# Estimate the inbreeding coefficient from all the autosomal SNPs
gcta64  --bfile test  --autosome  --ibc  --out test

# Calculate the GRM using the alternative method
gcta64  --bfile test  --autosome --make-grm  --make-grm-alg 1  --out test_alg1
```
 

#### Citations
 
**Method for estimating the GRM**: Yang et al. (2010) Common SNPs explain a large proportion of the heritability for human height. Nat Genet. 42(7): 565-9. [PubMed ID: 20562875]
 
**Method for estimating the inbreeding coefficients and GCTA software**: Yang J, Lee SH, Goddard ME and Visscher PM. GCTA: a tool for Genome-wide Complex Trait Analysis. Am J Hum Genet. 2011 Jan 88(1): 76-82. [PubMed ID: 21167468]

### Manipulating the GRM

Manipulation of the genetic relationship matrix

--grm test  
or  
--grm-bin test  
Input the GRM generated by --make-grm option. This option actually tells GCTA to read three files, e.g. test.grm.bin, test.grm.N.bin and test.grm.id (See the option --make-grm). GCTA automatically adds suffix ".grm.bin", ".grm.N.bin" or ".grm.id" to the specified root filename. If the test.grm.N.bin file (which contains the number of SNPs used to calculate GRM) is missing, the program will still be running because all the analysis except --grm do not actually need the the number of SNPs used to calculate the GRM.

--grm-gz test  
To be compatible with the previous version of GCTA. Same as --grm but read the GRM files in compressed text format generated by --make-grm-gz option. This option actually tells GCTA to read two files, e.g. test.grm.gz and test.grm.id (See the option --make-grm-gz). GCTA automatically adds suffix ".grm.gz" and ".grm.id" to the specified root filename. 

Examples: converting the two formats from each other
```bash
# From *.grm.gz to *.grm.bin
gcta64  --grm-gz test  --make-grm  --out test
# From *.grm.bin to *.grm.gz
gcta64  --grm test  --make-grm-gz --out test
```


--mgrm multi_grm.txt  
or  
--mgrm-bin multi_grm.txt  
Input multiple GRMs in binary format (See the option --make-grm). The root filenames of multiple GRMs are given in a file, e.g. multi_grm.txt  
> Input file format  
> multi_grm.txt (full paths can be specified if the GRM files are in different directories)
```nohighlight
test_chr1  
test_chr2  
test_chr3  
......  
test_chr22  
```

--unify-grm  
This option is designed to unify the individual IDs (as well as the order) of multiple GRMs (with --mgrm option) used in analyses such as REML and HE regression.

Examples  
```bash
gcta64 --mgrm multi_grm.txt --keep sample_list.txt --remove sample_rm.txt --unify-grm --out common
# The output are GRMs with the same individual IDs (i.e. the individuals in common among all the GRM files) in the same order.
```

--mgrm-gz multi_grm.txt  
To be compatible with the previous version of GCTA. Same as --mgrm but read the GRM files in compressed text format generated by --make-grm-gz.

Examples  
```bash
# This option is very useful to deal with large dataset. You can firstly run the jobs (split one job into 22 pieces)
gcta64  --bfile test  --chr 1  --make-grm  --out test_chr1
gcta64  --bfile test  --chr 2  --make-grm --out test_chr2
...
gcta64  --bfile test  --chr 22  --make-grm  --out test_chr22

# To estimate the GRMs from the SNPs on each chromosome, then merge them by the command
gcta64  --mgrm multi_grm.txt  --make-grm  --out test
```

--grm-cutoff 0.05  
Remove one of a pair of individuals with estimated relatedness larger than the specified cut-off value (e.g. 0.05). GCTA selectively removes individuals to maximize the remaining sample size rather than doing it at random. 
 
**Note**: 1) This flag has been rewritten to save memory usage. Currently, it can only be used in combination with other three flags, i.e., --grm --keep --remove and --make-grm.  
2) When merging multiple GRMs with --mgrm flag, this option does not apply to each single GRM but to the final merged GRM.

--grm-singleton 0.05  
Output IDs of individuals who do not have any relatives in sample given the relatedness threshold. This option will lead to two output files: \*.singleton.txt and \*.family.txt. It can be used in combination with --keep and --remove to manupulate the subjects.

> Format for \*.singleton.txt (FID IID)
```nohighlight
17      171
295     2951
429     4291
827     8271
2585    25851
...
```

> Format for \*.family.txt (FID1 IID1 FID2 IID2 GRM)
```nohighlight
5       51      3       31      0.129183
7       71      1       11      0.0732403
9       91      1       11      0.0618603
9       91      7       71      0.0703791
15      151     5       51      0.0623071
...
```

--grm-adj 0  
When using the SNPs to predict the genetic relationship at causal loci, we have to adjust the prediction errors due to imperfect LD because of two reasons: 1) the use of only a finite number of SNPs; 2) causal loci tend to have lower MAF than the genotyped SNPs (input 0 if you assume that the causal loci have similar distribution of allele frequencies as the genotyped SNPs) (see Yang et al. 2010 Nat Genet for details).

--dc 1  
By default, the GRM, especially for the X-chromosome, is parameterized under the assumption of equal variance for males and females, unless the option --dc is specified (1 and 0 for full and no dosage compensation, respectively). You need to use the option --update-sex to read sex information of the individuals from a file (see the --update-sex option above).

**NOTE**: you can add the option --make-grm or --make-grm-gz afterwards to save the modified GRM. You can also use the option --keep and/or --remove in combination with these five commands. It is also possible to use these five commands in the REML analysis (see the section below).

Examples  
```bash
# Prune the GRM for relatedness by a cutoff of 0.05
gcta64 --grm test --grm-cutoff 0.05 --make-grm --out test
# Extract the GRM subject id of all the singletons by a cutoff of 0.05
gcta64 --grm test --grm-singleton 0.05  --out test

# Use --keep or --remove option
gcta64  --grm test  --keep test.indi.list  --grm-cutoff  0.05  --make-grm  --out test_adj
gcta64  --grm test  --remove test.indi.list  --grm-adj 0  --make-grm  --out test_adj

# Assume full and no dosage compensation for the X chromosome
gcta64  --grm test_xchr  --dosage-compen 1  --update-sex test.indi.sex.list  --make-grm  --out test_xchr_fdc
gcta64  --grm test_xchr  --dosage-compen 0  --update-sex test.indi.sex.list  --make-grm  --out test_xchr_ndc
```


#### Citations

**Method for estimating the GRM**: Yang et al. (2010) Common SNPs explain a large proportion of the heritability for human height. Nat Genet. 42(7): 565-9. [PubMed ID: 20562875]

**Method for estimating the GRM for the X chromosome and GCTA software**: Yang J, Lee SH, Goddard ME and Visscher PM. GCTA: a tool for Genome-wide Complex Trait Analysis. Am J Hum Genet. 2011 Jan 88(1): 76-82. [PubMed ID: 21167468]

**A demonstration of estimating variance explained by the X chromosome for height and BMI**: Yang et al. (2011) Genome partitioning of genetic variation for complex traits using common SNPs. Nat Genet. 43(6): 519-525. [PubMed ID: 21552263]


### GREML analysis 

**GCTA-GREML: Estimate variance explained by all the SNPs**

--reml  
Perform a REML (restricted maximum likelihood) analysis. This option is usually followed by the option --grm (one GRM) or --mgrm (multiple GRMs) to estimate the variance explained by the SNPs that were used to estimate the GRM.

--reml-priors 0.45 0.55  
Specify the starting values for REML iterations. The number of starting values specified should NOT be smaller than the number of variance components in the model. By default, GCTA will use equal variances of all the components as the starting values if this option is not specified.

--reml-alg 0  
Specify the algorithm to run REML iterations, 0 for average information (AI), 1 for Fisher-scoring and 2 for EM. The default option is 0, i.e. AI-REML, if this option is not specified.

--reml-no-constrain  
By default, if an estimate of variance component escapes from the parameter space (i.e. negative value), it will be set to be a small positive value i.e. Vp * 10-6 with Vp being the phenotypic variance. If the estimate keeps escaping from the parameter space, the estimate will be constrained to be Vp * 10-6. If the option --reml-no-constrain is specified, the program will allow an estimate of variance component to be negative, which may result in the estimate of proportion variance explained by all the SNPs > 100%.

--reml-maxit 100  
Specify the maximum number of iterations. The default number is 100 if this option is not specified.

--pheno test.phen  
Input phenotype data from a plain text file, e.g. test.phen. If the phenotypic value is coded as 0 or 1, then it will be recognized as a case-control study (0 for controls and 1 for cases). Missing value should be represented by "-9" or "NA".  
> Input file format  
> test.phen (no header line; columns are family ID, individual ID and phenotypes)  
```nohighlight
011       0101    0.98   
012       0102    -0.76  
013       0103    -0.06  
......
```

--mpheno 2  
If the phenotype file contains more than one trait, by default, GCTA takes the first trait for analysis (the third column of the file) unless this option is specified. For example, --mpheno 2 tells GCTA to take the second trait for analysis (the fourth column of the file).

--gxe test.gxe  
Input an environmental factor from a plain text file, e.g. test.gxe. Apart from estimating the genetic variance, this command tells GCTA to estimate the variance of genotype-environment (GE) interaction. You can fit multiple environmental factors simultaneously. The main effects of an environmental factor will be included in the model as fixed effects and the GE interaction effects will be treated as random effects. NOTE: the design matrix of the overall mean in the model (which is a vector of all ones) is always a linear combination of the design matrix of a discrete environmental factor so that not all the main effects (fixed effects) are estimable. GCTA will always constrain the main effect of the first level to be zero and the main effect of any other level represents its difference in effect compared to the first level. For example, if you fit sex as an environmental factor, GCTA will fit only one main effect in the model, i.e. the mean difference between males and females.
> Input file format  
> test.gxe (no header line; columns are family ID, individual ID and environmental factors)  
```nohighlight
01       0101    F       smoker  
02       0203    M       nonsmoker  
03       0305    F       smoker  
......
```

--covar test.covar  
Input discrete covariates from a plain text file, e.g. test.covar. Each discrete covariate is recognized as a categorical factor with several levels. The levels of each factor can be represented by a single character, word or numerical number. NOTE: the design matrix of the mean in the model (which is a vector of all ones) is always a linear combination of the design matrix of a discrete covariate so that not all the effects of the levels (or classes, e.g. male and female) of a discrete covariate are estimable. GCTA will always constrain the effect of the first level to be zero and the effect of any other level represents its difference in effect compared to the first level.
> Input file format  
> test.covar (no header line; columns are family ID, individual ID and discrete covariates)  
```nohighlight
01       0101    F       Adult             0  
02       0203    M       Adult             0  
03       0305    F       Adolescent        1  
......
```

--qcovar test.qcovar  
Input quantitative covariates from a plain text file, e.g. test.qcovar. Each quantitative covariate is recognized as a continuous variable.  
> Input file format  
> test.qcovar (no header line; columns are family ID, individual ID and quantitative covariates)  
```nohighlight
01       0101    -0.024    0.012  
02       0203    0.032     0.106  
03       0305    0.143     -0.056  
......
```

--reml-res-diag  res\_diag\_file   
To specify the diagonal elements of the residual correlation matrix in REML (note that the original correlation matrix is an identity matrix). All the diagonal elements need to be positive.

> res\_diag\_file (no header line; Columns are family ID, individual ID and diagonal element)
```nohighlight
sub1 sub1 1.1
sub2 sub2 0.9
```

--reml-lrt 1  
Calculate the log likelihood of a reduce model with one or multiple genetic variance components dropped from the full model and calculate the LRT and p-value. By default, GCTA will always calculate and report the LRT for the first genetic variance component, i.e. --reml-lrt 1, unless you re-specify this option, e.g. --reml-lrt 2 assuming there are a least two genetic variance components included in the analysis. You can also test multiple components simultaneously, e.g. --reml-lrt 1 2 4. See FAQ #1 for more details. 

--reml-no-lrt  
Turn off the LRT.

--prevalence 0.01  
Specify the disease prevalence for a case-control study. Once this option is specified, GCTA will transform the estimate of variance explained, *V(1)/Vp*, on the observed scale to that on the underlying scale, *V(1)/Vp_L*. The prevalence should be estimated from a general population in literatures rather than that estimated from the sample. 

**NOTE**:  
1. You do not have to have exactly the same individuals in these files. GCTA will find the individuals in common in the files and sort the order of the individuals.  
2. Please be aware that if the GRM is estimated from the imputed SNPs (either "best guess" or "dosage score"), the estimate of variance explained by the SNPs will depend on the imputation-R<sup>2</sup> cutoff used to select SNPs because the imputation-R<sup>2</sup> is correlated with MAF, so that selection on imputation-R<sup>2</sup> will affect the MAF spectrum and thus affect the estimate of variance explained by the SNPs.  
3. For a case-control study, the phenotypic values of cases and controls should be specified as 1 and 0 (or 2 and 1, compatible with PLINK), respectively.  
4. Any missing value (either phenotype or covariate) should be represented by "-9" or "NA".  
5. The summary result of REML analysis will be saved in a plain text file (*.hsq).  

> Output file format  
> test.hsq (rows are
>> header line;  
>> name of genetic variance, estimate and standard error (SE);  
>> residual variance, estimate and SE;  
>> phenotypic variance, estimate and SE;  
>> ratio of genetic variance to phenotypic variance, estimate and SE;  
>> log-likelihood;  
>> sample size). If there are multiple GRMs included in the REML analysis, there will be multiple rows for the genetic variance (as well as their ratios to phenotypic variance) with the names of *V(1)*, *V(2)*, ... .  
```nohightlight
Source  Variance        SE
V(1)    0.389350        0.161719
V(e)    0.582633        0.160044
Vp      0.971984        0.031341
V(1)/Vp 0.400573        0.164937
The estimate of variance explained on the observed scale is transformed to that on the underlying scale:
(Proportion of cases in the sample = 0.5; User-specified disease prevalence = 0.1)
V(1)/Vp_L       0.657621      0.189123
logL    -945.65
logL0   -940.12
LRT     11.06
Pval    4.41e-4
n       2000
```


--reml-est-fix  
Output the estimates of fixed effects on the screen.

--reml-pred-rand  
Predict the random effects by the BLUP (best linear unbiased prediction) method. This option is actually to predict the total genetic effect (called "breeding value" in animal genetics) of each individual attributed by the aggregative effect of the SNPs used to estimate the GRM. The total genetic effects of all the individuals will be saved in a plain ext file *.indi.blp.
> Output file format  
> test.indi.blp (no header line; columns are family ID, individual ID, an intermediate variable, the total genetic effect, another intermediate variable and the residual effect.  
> If there are multiple GRMs fitted in the model, each GRM will insert additional two columns, , i.e. an intermediate variable (the intermediate variable = Py, please see Yang et al. 2011 AJHG for the definitions of P and y) and a total genetic effect, in front of the last two columns)  
```nohighlight
01       0101    -0.012    -0.014   -0.010    -0.035
02       0203    0.021     0.031    -0.027    -0.031
03       0305    0.097     0.102    -0.026    -0.041
......
```

--blup-snp test.indi.blp  
Calculate the BLUP solutions for the SNP effects (you have to specify the option --bfile to read the genotype data). This option takes the output of the option --reml-pred-rand as input (\*.indi.blp file) and transforms the BLUP solutions for individuals to the BLUP solutions for the SNPs, which can subsequently be used to predict the total genetic effect of individuals in an independent sample by PLINK --score option. Note that for the ease of using the BLUP solutions in a PLINK-score analysis, the BLUP effects are scaled by *sqrt[2p(1-p)]* (please see pages 77 and 78 of Yang et al. 2011 AJHG for details). 
> Output file format  
> test.snp.blp (columns are SNP ID, reference allele and BLUP of SNP effect; if there are multiple GRMs fitted in the model, each GRM will add an additional column to the file; the last column is for the residual effect)  
```nohighlight
rs103645   A     0.00312    0.00451
rs175292   G    -0.00021    0.00139
......
```

#### Examples

NOTE: if your GRMs files were generated by the --grm-bin option (i.e. saved in binary format, *.grm.bin), you could simply replace the --grm option by the --grm-bin option in the examples below.
```bash
# Without GRM (fitting the model under the null hypothesis that the additive genetic variance is zero)
gcta64  --reml  --pheno test.phen  --out test_null
gcta64  --reml  --pheno test.phen  --keep test.indi.list  --out test_null

# One GRM (quantitative traits)
gcta64  --reml  --grm test  --pheno test.phen  --reml-pred-rand –qcovar test_10PCs.txt  --out test
gcta64  --reml  --grm test  --pheno test.phen  --grm-adj 0  --grm-cutoff 0.05  --out test
gcta64  --reml  --grm test  --pheno test.phen  --keep test.indi.list  --grm-adj 0  --out test

# One GRM (case-control studies)
gcta64  --reml  --grm test  --pheno test_cc.phen  --prevalence 0.01  --out test_cc
gcta64  --reml  --grm test  --pheno test_cc.phen  --prevalence 0.01  --qcovar test_10PCs.txt  --out test_cc

# GxE interaction (LRT test for the significance of GxE)
gcta64  --reml  --grm test  --pheno test.phen   --gxe test.gxe  --reml-lrt 2  --out test

# Multiple GRMs
gcta64  --reml  --mgrm multi_grm.txt  --pheno test.phen  --reml-no-lrt  --out test_mgrm 
gcta64  --reml  --mgrm multi_grm.txt  --pheno test.phen  --keep test.indi.list  --reml-no-lrt  --out test_mgrm 

# BLUP solutions for the SNP effects
gcta64  --bfile   test   --blup-snp test.indi.blp  --out test
# Then use plink --score test.snp.blp 1 2 3 
```

--reml-bendV  
The GREML method uses REML for variance estimation, which requires the inverse of the variance-covariance matrix V. If V is not positive definite, the inverse of V does not exist. We therefore could not estimate the variance component. This usually happens when one (or more) of the variance components are negative or constrained at zero. It might also indicate there is something wrong with the GRM or the data which you might need to check carefully.

Unfortunately, there has not been an ultimate solution. Tricks such as adding a small number of to the diagonal elements of V also do not guarantee the modified V being invertible. In some cases, you might be able to get around the problem by using alternative REML algorithms e.g. the Fisher scoring approach (--reml-alg 1).

We have implemented the "bending" approach (Hayes and Hill 1981 Biometrics) in GCTA to invert V if V is not positive definite (you could add the --reml-bendV option to a REML or MLMA analysis to activate this approach). The "bending" approach guarantees to get an approximate of V-1 but it does not guarantee the REML analysis being converged.

Note that the --reml-bendV option only provides an approximate inverse of V and has not been tested extensively. The results from analyses using this option might not be reliable.


#### Citations

**Method for estimating the variance explained by all SNPs**: Yang et al. (2010) Common SNPs explain a large proportion of the heritability for human height. Nat Genet. 42(7): 565-9. [PubMed ID: 20562875]

**Method for estimating the variance explained by all SNPs using case-control data**: Lee et al. (2011) Estimating Missing Heritability for Disease from Genome-wide Association Studies. Am J Hum Genet. 88(3): 294-305. [PubMed ID: 21376301]

**Method for partitioning the genetic variance captured by all SNPs onto chromosomes and genomic segments**: Yang et al. (2011) Genome partitioning of genetic variation for complex traits using common SNPs. Nat Genet. 43(6): 519-525. [PubMed ID: 21552263]


### GREML in family data

**GCTA-GREML analysis in family data**

Zaitlen et al. [2013 PLoS Genetics](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1003520) proposed a method to estimate pedigree-based and SNP-based *h<sup>2</sup>* simultaneously in one model using family data. The main advantage of this method is that it allows us to estimate SNP-based *h<sup>2</sup>* in family data without having to remove related individuals. The method has now been implemented in GCTA.

--make-bK 0.05  
The default value is 0.05. This option will set the GRM off-diagonal elements that are below the threshold to 0. It has been updated to save the memory (Memory usage less than 500MB).

#### Examples
```bash
# Making a GRM from all SNPs in a family data set
gcta64 --bfile test --make-grm --out test

# Creating an additional GRM from the GRM above (setting the off-diagonals that are < 0.05 to 0)
gcta64 --grm test --make-bK 0.05 --out test_bK
```

> An example of the mgrm.txt file
```nohighlight
test  
test_bK
```

```bash
# Running a REML analysis with two GRMs
gcta64 --reml --mgrm mgrm.txt --pheno test.phen --out test_bKsK
```

> Here is an example of the output file (test_bKsK.hsq)
```nohighlight
Source	Variance	SE
V(G1)	0.294615	0.102976
V(G2)	0.322424	0.144884
V(e)	0.377467	0.104458
Vp	0.994506	0.027059
V(G1)/Vp	0.296242	0.102655
V(G2)/Vp	0.324205	0.145112

Sum of V(G)/Vp	0.620447	0.105741
logL	-1357.892
n	2753
```

where "*V(G<sub>1</sub>) / V<sub>p</sub>*" provides an estimate of SNP-based *h<sup>2</sup>* (h<sup>2</sup>SNP), "Sum of *V(G) / V<sub>p</sub>*" provides an estimate of pedigree-based *h<sup>2</sup>* (*h<sup>2</sup>ped*), and *V(G<sub>2</sub>) / V<sub>p</sub> = h<sup>2</sup>ped - h<sup>2</sup>SNP*.


#### Citations

**Method for estimating the GRM**: Yang et al. (2010) Common SNPs explain a large proportion of the heritability for human height. Nat Genet. 42(7): 565-9. [PubMed ID: 20562875]

**The Zaitlen et al. method**: Zaitlen N, Kraft P, Patterson N, Pasaniuc B, Bhatia G, Pollack S, Price AL (2013) Using extended genealogy to estimate components of heritability for 23 quantitative and dichotomous traits. PLoS Genet. 2013 May;9(5):e1003520. PubMed ID: 23737753]

**REML analysis and GCTA Software**: Yang J, Lee SH, Goddard ME and Visscher PM. GCTA: a tool for Genome-wide Complex Trait Analysis. Am J Hum Genet. 2011 Jan 88(1): 76-82. [PubMed ID: 21167468]


### GREML in WGS or imputed data

**GCTA-LDMS: estimating heritability from WGS or imputed WGS data**

The GREML-LDMS method is proposed to estimate heritability using whole genome sequence (WGS) data ([Yang et al. 2015 Nature Genetics](http://www.nature.com/ng/journal/v47/n10/full/ng.3390.html)). **It corrects for the LD bias in the estimated SNP-based heritability.** It can also be applied to (imputed) GWAS data. The method is unbiased regardless the properties (e.g. MAF and LD) of the underlying causal variants. The analysis involves four steps.  
1) calculating segment-based LD score;  
2) stratifying SNPs based on the segment-based LD score (this is done in R);  
3) computing GRMs using the stratified SNPs;  
4) performing REML analysis using the multiple GRMs.  

**Tutorial**:

#### Step 1: segment based LD score
```bash
gcta64 --bfile test --ld-score-region 200 --out test
```

--ld-score-region 200  
The default value is 200Kb, i.e. the length of the segment is 200Kb (with 100Kb overlap between two adjacent segments). Results are save a *.score.ld file.

> Output file format
test.score.ld(Columns are SNP ID, chromosome, physical position, allele frequency, mean LD rsq between the target SNP and other SNPs in a window (specified by the --ld-wind option), number of SNPs in LD with the target SNP passing the threshold (specified by the --ld-rsq-cutoff option), maximum rsq between the target SNP and its best tagging SNP within the window, LD score of the SNP, and LD score of the region). 
```nohighlight
SNP chr bp freq mean_rsq snp_num max_rsq ldscore_SNP ldscore_region
rs4475691 1 836671 0.197698 0.000588093 999 0.216874 1.5875 2.75538
rs28705211 1 890368 0.278112 0.000573876 999 0.216874 1.5733 2.75538
rs9777703 1 918699 0.0301614 0.00131291 999 0.854464 2.31159 2.75538
.... 
```

#### Step 2 (option #1): stratify the SNPs by LD scores of individual SNPs in R
Below is an example of R script to stratify the SNPs by the LD scores of individual SNPs.
```r
lds_seg = read.table("test.score.ld",header=T,colClasses=c("character",rep("numeric",8)))
quartiles=summary(lds_seg$ldscore_SNP)

lb1 = which(lds_seg$ldscore_SNP <= quartiles[2])
lb2 = which(lds_seg$ldscore_SNP > quartiles[2] & lds_seg$ldscore_SNP <= quartiles[3])
lb3 = which(lds_seg$ldscore_SNP > quartiles[3] & lds_seg$ldscore_SNP <= quartiles[5])
lb4 = which(lds_seg$ldscore_SNP > quartiles[5])

lb1_snp = lds_seg$SNP[lb1]
lb2_snp = lds_seg$SNP[lb2]
lb3_snp = lds_seg$SNP[lb3]
lb4_snp = lds_seg$SNP[lb4]

write.table(lb1_snp, "snp_group1.txt", row.names=F, quote=F, col.names=F)
write.table(lb2_snp, "snp_group2.txt", row.names=F, quote=F, col.names=F)
write.table(lb3_snp, "snp_group3.txt", row.names=F, quote=F, col.names=F)
write.table(lb4_snp, "snp_group4.txt", row.names=F, quote=F, col.names=F)
```

In each LD group, you can use the --maf and --max-maf options GCTA to further stratify the SNPs into MAF groups.

Option 1 is preferred because it is shown in Evans et al. ([2018 Nature Genetics](https://www.nature.com/articles/s41588-018-0108-x)) stratifying SNPs by LD scores of individual SNPs as opposed to regional LD scores significantly improves the performance of GREML-LDMS.

#### Step 2 (option #2): stratify the SNPs by segment-based LD scores in R
This assumes that functional SNPs are clustered in regions with higher or lower LD.  
Below is an example of R script to stratify the SNPs by the segment-based mean LD scores.
```r
lds_seg = read.table("test.score.ld",header=T,colClasses=c("character",rep("numeric",8)))
quartiles=summary(lds_seg$ldscore_region)

lb1 = which(lds_seg$ldscore_region <= quartiles[2])
lb2 = which(lds_seg$ldscore_region > quartiles[2] & lds_seg$ldscore_region <= quartiles[3])
lb3 = which(lds_seg$ldscore_region > quartiles[3] & lds_seg$ldscore_region <= quartiles[5])
lb4 = which(lds_seg$ldscore_region > quartiles[5])

lb1_snp = lds_seg$SNP[lb1]
lb2_snp = lds_seg$SNP[lb2]
lb3_snp = lds_seg$SNP[lb3]
lb4_snp = lds_seg$SNP[lb4]

write.table(lb1_snp, "snp_group1.txt", row.names=F, quote=F, col.names=F)
write.table(lb2_snp, "snp_group2.txt", row.names=F, quote=F, col.names=F)
write.table(lb3_snp, "snp_group3.txt", row.names=F, quote=F, col.names=F)
write.table(lb4_snp, "snp_group4.txt", row.names=F, quote=F, col.names=F)
```

#### Step 3: making GRMs using SNPs stratified into different groups
```bash
gcta64 --bfile test --extract snp_group1.txt --make-grm --out test_group1
gcta64 --bfile test --extract snp_group2.txt --make-grm --out test_group2
...
```

#### Step 4: REML analysis with multiple GRMs
```bash
gcta64 --reml --mgrm multi_GRMs.txt --pheno phen.txt --out test
```

> format of multi_grm.txt (no headline; each line represents the prefix of a GRM file)
```nohighlight
test_group1
test_group2
...
```


#### Citations:

**Method paper**: Yang et al. (2015) Genetic variance estimation with imputed variants finds negligible missing heritability for human height and body mass index. Nature Genetics, 47:1114–1120.

**GCTA software**: Yang J, Lee SH, Goddard ME and Visscher PM. GCTA: a tool for Genome-wide Complex Trait Analysis. Am J Hum Genet. 2011 Jan 88(1): 76-82. [PubMed ID: 21167468]


### GREML for dominance variance

**GCTA-GREMLd: estimating dominance variance in unrelated individuals using SNP data**

Details of the method can be found in [Zhu et al. (2015 AJHG)](http://www.sciencedirect.com/science/article/pii/S0002929715000099).

--make-grm-d  
or  
--make-grm-d-bin  
Estimate the dominance genetic relationship matrix (GRM) between pairs of individuals from a set of SNPs and save the lower triangle elements of the dominance GRM to binary files. eg. test.grm.d.bin, test.grm.d.N.bin, test.grm.d.id.

Note: the memory usage of --make-grm-d is same with --make-grm. --make-grm-d-bin takes much more memory.

> Output file format:
test.grm.d.bin Binary file which contains the lower triangle elements of the dominance GRM).  
test.grm.d.N.bin Binary file which contains the number of SNPs used to calculate the dominance GRM).  
test.grm.d.id No header line; columns are family ID and individual ID

--make-grm-d-part m i  
Partition the dominance GRM into m parts (by row), and compute the i-th part in the current run.

See the document of [--make-grm-part](#MakingaGRM)

--make-grm-d-gz  
Estimate the dominance GRM, save the lower triangle elements to a compressed text file (e.g. test.grm.d.gz) and save the IDs in a plain text file (e.g. test.grm.d.id). 
> Output format:
> test.grm.d.gz (No header line; columns are indices of pairs of individuals (row numbers of the test.grm.d.id), number of non-missing SNPs and the estimate of dominance genetic relatedness)
```nohighlight
1    1     1000       0.0021       
2    1      998       0.0231       
2    2      999       0.0238       
3    1     1000       0.0031       
..... 
```

> test.grm.d.id (no header line; columns are family ID and individual ID)
```nohighlight
011    0101     
012    0102     
013    0103     
..... 
```

> Examples:
```bash
# Calculating the additive GRM from all the autosomal SNPs
gcta64 --bfile test --autosome --make-grm --thread-num 10 --out test_add

# Calculating the dominance GRM from all the autosomal SNPs
gcta64 --bfile test --autosome --make-grm-d --thread-num 10 --out test_domi

# Estimating additive and dominance genetic variance by fitting an AD model
gcta64 --reml --mgrm add_domi_grm.txt --pheno test.phen --thread-num 10 --out test_add_domi

# format of add_domi_grm.txt (no headline; each line represents the prefix of a GRM file)
test_add
test_domi

# Note that most of the other GCTA options (e.g. --extract and --keep) are also valid for these analyses
```


#### Citations:

**Method paper**: Zhu Z, Bakshi A, Vinkhuyzen AA, Hemani G, Lee SH, Nolte IM, van Vliet-Ostaptchouk JV, Snieder H, The LifeLines Cohort Study, Esko T, Milani L, Mägi R, Metspalu A, Hill WG, Weir BS, Goddard ME, Visscher PM, Yang J (2015) Dominance Genetic Variation Contributes Little to the Missing Heritability for Human Complex Traits. Am J Hum Genet, 96: 1-9. [PubMed ID: 25683123]

**GCTA software**: Yang J, Lee SH, Goddard ME and Visscher PM. GCTA: a tool for Genome-wide Complex Trait Analysis. Am J Hum Genet. 2011 Jan 88(1): 76-82. [PubMed ID: 21167468]


### Bivariate GREML analysis

These options are designed to perform a bivariate GREML analysis to estimate the genetic correlation between two quantitative traits, between two disease traits (binary) from case control studies, and between a quantitative trait and a binary disease trait. The analysis will output the estimated genetic variance for each trait (captured by the SNPs) and the genetic covariance betwen traits.  

--reml-bivar 1 2  
By default, GCTA will take the first two traits in the phenotype file for analysis. The phenotype file is specified by the option --pheno as described in univariate REML analysis. All the options for univariate REML analysis are still valid here except --mpheno, --gxe, --prevalence, --reml-lrt, --reml-no-lrt and --blup-snp. All the input files are in the same format as in univariate REML analysis.

--reml-bivar-nocove  
By default, GCTA will model the residual covariance between two traits. However, if the traits were measured on different individuals (e.g. two diseases), the residual covariance will be automatically dropped from the model. You could also specify this option to exclude the residual covariance at all time.

--reml-bivar-lrt-rg 0  
To test for the hypothesis of fixing the genetic correlation at a particular value, e.g. fixing genetic correlation at -1, 0 and 1. By default bivariate GCTA-GREML does not perform a log likelihood test unless this option is specified.

--reml-bivar-prevalence 0.1 0.05  
For a bivariate analysis of two disease traits, you can specify the prevalence rates of the two diseases in the general population so that GCTA will transform the estimate of variance explained by the SNPs from the observed 0-1 scale to that on the underlying scale for both diseases.

--reml-bivar-no-constrain  
By default, the genetic correlation estimate is constrained between -1 and 1. This option will allow the estimate of genetic correlation > 1 or < -1. Note that not all the analyses can converge with this option.

> Examples
```bash
# With residual covariance
gcta64  --reml-bivar --grm test  --pheno test.phen  --out test

# Without residual covariance
gcta64  --reml-bivar --reml-bivar-nocove --grm test  --pheno test.phen  --out test

# To test for genetic correlation = 0 or 1
gcta64  --reml-bivar --reml-bivar-nocove --grm test  --pheno test.phen --reml-bivar-lrt-rg 0  --out test

gcta64  --reml-bivar --reml-bivar-nocove --grm test  --pheno test.phen --reml-bivar-lrt-rg 1  --out test

# Case-control data for two diseases (the residual covariance will be automatically dropped from the model if there are not too many samples affected by both diseases)
gcta64  --reml-bivar  --grm test_CC  --pheno test_CC.phen --reml-bivar-prevalence 0.1 0.05  --out test_CC

# Bivariate GREML analysis with multiple GRMs
gcta64  --reml-bivar --mgrm multi_grm.txt  --pheno test.phen  --out test
```

See <a href="#ManipulationoftheGRM" target="_blank">Manipulation of the GRM</a> for the format of multi_grm.txt.

> Output file format  
> test.hsq (rows are
> * header line;  
> * genetic variance for trait 1, estimate and standard error (SE);  
> * genetic variance for trait 2, estimate and SE;  
> * genetic covariance between traits 1 and 2, estimate and SE;  
> * residual variance for trait 1, estimate and SE;  
> * residual variance for trait 2, estimate and SE;  
> * residual covariance between traits 1 and 2, estimate and SE;  
> * proportion of variance explained by all SNPs for trait 1, estimate and SE;  
> * proportion of variance explained by all SNPs for trait 2, estimate and SE;  
> * genetic correlation;  
> * sample size).  
```nohighlight
Source      Variance    SE
V(G)_tr1    0.479647    0.179078
V(G)_tr2    0.286330    0.181329
C(G)_tr12   0.230828    0.147958
V(e)_tr1    0.524264    0.176650
V(e)_tr2    0.734654    0.181146
C(e)_tr12   0.404298    0.146863
Vp_tr1      1.003911    0.033202
Vp_tr2      1.020984    0.033800
V(G)/Vp_tr1  0.477779    0.176457
V(G)/Vp_tr2  0.280445    0.176928
rG    0.622864    0.217458
n    3669
```


#### Citations

**The first bivariate GREML example**: Deary et al. (2012) Genetic contributions to stability and change in intelligence from childhood to old age. Nature, 482: 212-215. [Pubmed ID: 22258510]

**Bivariate GREML analysis method**: Lee et al. (2012) Estimation of pleiotropy between complex diseases using SNP-derived genomic relationships and restricted maximum likelihood. Bioinformatics, 28: 2540-2542. [PubMed ID: 22843982]

**GCTA software**: Yang J, Lee SH, Goddard ME and Visscher PM. (2011) GCTA: a tool for Genome-wide Complex Trait Analysis. Am J Hum Genet, 88: 76-82. [PubMed ID: 21167468]


### GREML power calculator
R functions for GCTA power calculation

Although the online version of [GCTA power calculator](http://shiny.cnsgenomics.com/gctaPower/) is available, I think it would still be useful for some users to have the R functions below.

Reference: [Visscher et al. (2014) Statistical power to detect genetic (co)variance of complex traits using SNP data in unrelated samples. PLoS Genetics, 10(4): e1004269.](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004269)
```r
# Function for a quantitative trait
# n = sample size 
# hsq = variance explained by all SNPs
# alpha = significance level
# var_pi = variance of the off-diagonal elements of the GRM
# The output are: se (standard error), ncp (non-centrality parameter) and power
calcUniQt <- function(
	n     =1000, 
	hsq   =0.5, 
	alpha =0.05,
  var_pi=2e-5
){
	l <- list()
	var_vg <- var_vg_func(n, var_pi)
	l$se <- sqrt(var_vg)
	l$ncp <- hsq^2/var_vg;
	l$power <- power_func(l$ncp, alpha)
	return(l)
}
```
```r
# Function for case-control study
# ncase = number of cases
# ncontrol = number of controls
# K = disease prevalence in the population
calcUniCc <- function(
	ncase    = 1000, 
	ncontrol = 1000, 
	hsq      = 0.5, 
	K        = 0.1, 
	alpha    = 0.05,
	var_pi=2e-5
){
	h <- h2O_func(ncase, ncontrol, K, hsq, var_pi)
	l <- list()
	l$se <- sqrt(h$var_h2L)
	l$ncp <- h$h2L^2/h$var_h2L
	l$power <- power_func(l$ncp, alpha)
	return(l)
}
```
```r
# Function for bivariate analysis of two quantitative traits
# rg = genetic correlation
# rp = phenotypic correlation
# overlap = whether or not the traits are measured on the same samples
calcBiQt <- function(
	n1      = 1000, 
	n2      = 1000, 
	hsq1    = 0.5, 
	hsq2    = 0.5, 
	rg      = 0.5, 
	rp      = 0.5, 
	overlap = FALSE, 
	alpha   = 0.05,
	var_pi=2e-5
){
	var_rg <- var_rg_func(n1, n2, hsq1, hsq2, rg, rp, overlap, var_pi)
	l <- list()
	l$se <- sqrt(var_rg)
	l$ncp <- rg^2/var_rg;
	l$power <- power_func(l$ncp, alpha)
	return(l)
}
```
```r
# Function for bivariate analysis of two case-control studies
calcBiCc <- function(
	ncase1    = 1000, 
	ncase2    = 1000, 
	ncontrol1 = 1000, 
	ncontrol2 = 1000, 
	hsq1      = 0.5, 
	hsq2      = 0.5, 
	K1        = 0.1, 
	K2        = 0.1, 
	rg        = 0.5, 
	overlap   = FALSE, 
	alpha     = 0.05,
	var_pi=2e-5
){
	h1 <- h2O_func(ncase1, ncontrol1, K1, hsq1, var_pi)
	h2 <- h2O_func(ncase2, ncontrol2, K2, hsq2, var_pi)
	n1 <- ncase1+ncontrol1
	n2 <- ncase2+ncontrol2
	var_rg <- var_rg_func(n1, n2, h1$h2O, h2$h2O, rg, rg, overlap, var_pi)
	l <- list()
	l$se <- sqrt(var_rg)
	l$ncp <- rg^2/var_rg;
	l$power <- power_func(l$ncp, alpha)
	return(l)
}
```
```r
# Function for bivariate analysis of a quantitative trait and a binary trait (case-control study)
calcBiQtCc <- function(
	n        = 1000, 
	ncase    = 1000, 
	ncontrol = 1000, 
	hsq1     = 0.5, 
	hsq2     = 0.5, 
	K        = 0.1, 
	rg       = 0.5, 
	overlap  = FALSE, 
	alpha    = 0.05,
	var_pi=2e-5
){
	h2=h2O_func(ncase, ncontrol, K, hsq2, var_pi)
	n2=ncase+ncontrol
	var_rg=var_rg_func(n, n2, hsq1, h2$h2O, rg, rg, overlap, var_pi)
	l <- list()
	l$se <- sqrt(var_rg)
	l$ncp <- rg^2/var_rg;
	l$power <- power_func(l$ncp, alpha)
	return(l)
}
```
```r
# Functions used in the functions above
var_vg_func <- function(N, var_pi=2e-5){
	return(2/(N^2*var_pi))
}

var_rg_func <- function(N1, N2, hsq1, hsq2, rg, rp, overlap=TRUE, var_pi=2e-5){
	if(overlap==T) var_rg=((1-rg*rp)^2+(rg-rp)^2)/(hsq1*hsq2*N1^2*var_pi)
	if(overlap==F) var_rg=(rg^2*(N1^2*hsq1^2+N2^2*hsq2^2)+2*hsq1*hsq2*N1*N2)/(2*hsq1^2*hsq2^2*N1^2*N2^2*var_pi)
	return(var_rg)
}

power_func <- function(ncp, alpha){
	pchisq(qchisq(alpha, df=1,lower.tail=F), ncp=ncp, df=1, lower.tail=F)
}

h2O_func <- function(ncase, ncontrol, K, h2L, var_pi=2e-5){
	n=ncase+ncontrol
	v=ncase/(ncase+ncontrol)
	z=dnorm(qnorm(K))
	c=(K*(1-K))^2/(v*(1-v)*z^2)
	h2O=h2L/c
	var_h2O=var_vg_func(n, var_pi)
	var_h2L=c^2*var_h2O
	return(list(h2L=h2L, var_h2L=var_h2L, h2O=h2O, var_h2O=var_h2O))
}
```

### Haseman-Elston regression

**GCTA-HEreg: Haseman-Elston regression analysis**

Haseman-Elston (HE) regression is a moment-based method for estimating the heritability. It is computationally much more efficient but slightly less powerful than REML as the SE of the estimate from HE regression is larger than that from REML. We implemented a HE regression that allows fitting multiple GRMs and facilitates bivariate analysis as in the GREML analysis, and only requires a small amount of memory (e.g. <2GB for n=120,000). The bivariate analysis is essentially three sets of independent HE regression for the variance and covariance components, where the sampling variance/covariance of the estimates (including the genetic correlation rG) are computed using leave-one-individual-out Jackknife technique.

> Example
```bash
# One GRM
gcta64 --HEreg --grm test --pheno test.phen --out test

# Multiple GRMs
gcta64 --HEreg --mgrm multi_grm.txt --pheno test.phen --out test

# Bivariate analysis with one GRM
gcta64 --HEreg-bivar 1 2 --grm test --pheno test.phen --out test

# Bivariate analysis with multiple GRMs
gcta64 --HEreg-bivar 1 2 --mgrm multi_grm.txt --pheno test.phen --out test
```

> Output results are saved in *.HEreg file.

> Univariate analysis with one GRM
> * HE-SD: HE regression using the square difference of the phenotypes for pairwise individuals  
> * HE-CP: HE regression using the cross product of the phenotypes for pairwise individuals
> * SE_OLS: standard error estimated from the ordinary least squares, which is likely to be an underestimation in a large sample
> * SE_Jackknife: standard error estimated using leave-one-individual-out Jackknife technique
```nohighlight
HE-CP
Coefficient     Estimate        SE_OLS          SE_Jackknife    P_OLS           P_Jackknife     
Intercept       -2.8813e-06     1.11965e-05     8.16333e-08     0.79692         6.83316e-273    
V(G)/Vp         0.710794        0.00260404      0.0118354       0               0               

HE-SD
Coefficient     Estimate        SE_OLS          SE_Jackknife    P_OLS           P_Jackknife     
Intercept       -0.999995       1.62064e-05     0.00416409      0               0               
V(G)/Vp         0.709868        0.0037692       0.0121754       0               0               
```

> Bivariate analysis with two GRMs
```nohighlight
HE-CP
Coefficient         Estimate        SE_OLS          SE_Jackknife    P_OLS           P_Jackknife     
Intercept_tr1       -1.73059e-05    2.77668e-05     2.7531e-07      0.533114        0               
Intercept_tr2       -1.4924e-05     2.50997e-05     2.38696e-07     0.552119        0               
Intercept_tr12      2.94888e-06     1.86676e-05     8.62647e-06     0.874481        0.732469        
V(G1)/Vp_tr1        0.65882         0.0064606       0.0175264       0               3.12495e-309    
V(G1)/Vp_tr2        0.696183        0.00584055      0.0174994       0               0               
C(G1)/Vp_tr12       0.660164        0.00434369      0.0115293       0               0               
V(G2)/Vp_tr1        0.0140144       0.00116567      0.00228244      2.70237e-33     8.24779e-10     
V(G2)/Vp_tr2        0.0112924       0.00105554      0.00201111      1.03724e-26     1.96549e-08     
C(G2)/Vp_tr12       0.012702        0.000784421     0.0014176       5.66913e-59     3.24039e-19     
Sum of V(G)/Vp_tr1  0.672834        0.00656473      0.0178002       0               1.16416e-312    
Sum of V(G)/Vp_tr2  0.707476        0.0059348       0.0177104       0               0               
Sum of C(G)/Vp_tr12 0.672866        0.00441375      0.011687        0               0               
rG1                 0.97478         0.00898327      0.0125365       
rG2                 1.0097          0.0887597       0.114919        
Total rG            0.975256        0.0089607       0.0124629       
N_tr1               50930           
N_tr2               56342           
```


#### Citation

Yang J, Zeng J, Goddard ME, Wray NR, Visscher PM  (2017) Concepts, estimation and interpretation of SNP-based heritability. Nature Genetics, 49: 1304-1310.
