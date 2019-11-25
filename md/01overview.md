## Overview {: .expand}

### About
GCTA (Genome-wide Complex Trait Analysis) was initially designed to estimate the proportion of phenotypic variance explained by all genome-wide SNPs for complex traits (i.e., the GREML method). It has been subsequently extended for many other analyses to better understand the genetic architecture of complex traits. GCTA currently supports the following analyses.
* GRM: estimating the genetic relationships among individuals in GWAS data;
* Estimating the inbreeding coefficients of individuals in GWAS data;
* GREML: estimating the proportion of variance in a phenotype explained by all GWAS SNPs (i.e. the SNP-based heritability);
* Partitioning the genetic variance onto individual chromosomes, MAF bins or functional categories;
* Estimating the genetic variance attributed to the X chromosome, and testing for the effect of dosage compensation;
* GREMLd: estimating the dominance variance in unrelated individuals using GWAS data;
* Bivariate GREML: estimating the genetic correlation between two traits (diseases) using GWAS data;
* PCA analysis and estimation of Fst in GWAS data;
* Computing LD scores and searching for LD friends for a list of target SNPs;
* Simulating a phenotype based on GWAS data;
* Conditional & joint (COJO) analysis of GWAS summary statistics without individual-level genotype data;
* mtCOJO: multi-trait-based conditional & joint analysis using GWAS summary data;
* GSMR: generalised summary-data-based mendelian randomisaion;
* MLMA and MLMA-LOCO: mixed linear model association analysis;
* fastBAT: gene- or set-based association analysis;
* sBLUP: summary-data based BLUP analysis for genomic risk prediction;
* Haseman-Elston regression to estimate the the SNP-based heritability for a trait and the genetic correlation between two traits;
* fastGWA: an extremely source-efficient (mixed) linear model association tool.

**Latest release [v1.92.4beta2](#Download), click to download or view update log (25 Nov 2019)**

### Credits 

[Jian Yang](http://scholar.google.com.au/citations?user=aLuqQs8AAAAJ&hl=en) developed the original version of the software with supports from [Peter Visscher](mailto:peter.visscher@uq.edu.au), [Mike Goddard](mailto:Mike.Goddard@dpi.vic.gov.au) and [Hong Lee](mailto:hong.lee@uq.edu.au). 

[Zhili Zheng](mailto:zhili.zheng@uq.edu.au) programmed the fastGWA module, rewrote the I/O and GRM modules, improved the GREML and bivariate GREML modules, extended the GCTA-PCA module, improved the SBLUP module, developed the website, and is currently maintaining the software.

[Zhihong Zhu](mailto:z.zhu1@uq.edu.au) programmed the GCTA-mtCOJO, GCTA-GSMR module and improved the GCTA-COJO module. 

[Jian Zeng](mailto:j.zeng@imb.uq.edu.au) rewrote the GCTA-HEreg module. 

[Andrew Bakshi](mailto:andrew.bakshi@gmail.com) contributed to the GCTA-fastBAT module. 

[Robert Maier](mailto:rmaier@broadinstitute.org) improved the GCTA-SBLUP module.

Contributions to the developments of methods included in GCTA (e.g., GREML methods, COJO, MLMA-LOCO, fastBAT and fastGWA) can be found in the papers cited in the corresponding webpages.

### Questions and Help Requests 
If you have any bug reports or questions please send an email to [Jian Yang](http://scholar.google.com.au/citations?user=aLuqQs8AAAAJ&hl=en) at [jian.yang@uq.edu.au](mailto:jian.yang@uq.edu.au)

### Citations 
**GCTA Software tool:**  
Yang et al. (2011) GCTA: a tool for Genome-wide Complex Trait Analysis. Am J Hum Genet. 88(1): 76-82. \[[PubMed ID: 21167468](http://www.ncbi.nlm.nih.gov/pubmed/21167468)\]

**Method for estimating the variance explained by all SNPs (GREML method) with its application in human height:**  
Yang et al. (2010) Common SNPs explain a large proportion of the heritability for human height. Nat Genet. 42(7): 565-9. \[[PubMed ID: 20562875](http://www.ncbi.nlm.nih.gov/pubmed/20562875)\]

**GREML method being extended for case-control design with its application to the WTCCC data:**  
Lee et al. (2011) Estimating Missing Heritability for Disease from Genome-wide Association Studies. Am J Hum Genet. 88(3): 294-305. \[[PubMed ID: 21376301](http://www.ncbi.nlm.nih.gov/pubmed?term=Estimating%20Missing%20Heritability%20for%20Disease%20from%20Genome-wide%20Association%20Studies)\]

**Extension of GREML method to partition the genetic variance into individual chromosomes and genomic segments with its applications in height, BMI, vWF and QT interval:**  
Yang et al. (2011) Genome partitioning of genetic variation for complex traits using common SNPs. Nat Genet. 43(6): 519-525. \[[PubMed ID: 21552263](http://www.ncbi.nlm.nih.gov/pubmed/21552263)\]

**Method for conditional and joint analysis using summary statistics from GWAS with its application to the GIANT meta-analysis data for height and BMI:**  
Yang et al. (2012) Conditional and joint multiple-SNP analysis of GWAS summary statistics identifies additional variants influencing complex traits. Nat Genet 44(4):369-375. \[[PubMed ID: 22426310](http://www.ncbi.nlm.nih.gov/pubmed/22426310)\]

**Bivariate GREML method:**  
Lee et al. (2012) Estimation of pleiotropy between complex diseases using SNP-derived genomic relationships and restricted maximum likelihood. Bioinformatics. 28(19): 2540-2542. \[[PubMed ID: 22843982](http://www.ncbi.nlm.nih.gov/pubmed/22843982)\]

**Mixed linear model based association analysis:**  
Yang et al. (2014) Mixed model association methods: advantages and pitfalls. Nat Genet. 2014 Feb;46(2):100-6. \[[Pubmed ID: 24473328](http://www.ncbi.nlm.nih.gov/pubmed/24473328)\]

**GREML-LDMS method and LD-score calculation:**  
Yang et al. (2015) Genetic variance estimation with imputed variants finds negligible missing heritability for human height and body mass index. Nat Genet. 47(10):1114-20.\[[PMID: 26323059](http://www.nature.com/ng/journal/vaop/ncurrent/full/ng.3390.html)\]

**Method to search for LD friends:**  
Yang et al. (2011) Genomic inflation factors under polygenic inheritance. Eur J Hum Genet. 19(7): 807-812. \[[Pubmed ID: 21407268](http://www.nature.com/ejhg/journal/v19/n7/full/ejhg201139a.html)\]

**fastBAT method:**  
Bakshi et al. (2016) Fast set-based association analysis using summary data from GWAS identifies novel gene loci for human complex traits. Scientific Reports 6, 32894. \[[PMID: 27604177](https://www.nature.com/articles/srep32894)\]

**mtCOJO and GSMR methods:**  
Zhu et al. (2018) Causal associations between risk factors and common diseases inferred from GWAS summary data. Nat Commun. 9, 224.\[[PMID: 29335400](https://www.ncbi.nlm.nih.gov/pubmed/?term=29335400)\]  

**fastGWA method:**  
Jiang et al. (2019) A resource-efficient tool for mixed model association analysis of large-scale data. bioRxiv 598110; doi:10.1101/598110.


<p style="color: rgb(204,0,0);font-weight:bold;">Last update: 25 Nov 2019</p>

