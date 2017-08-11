## Overview {: .expand}

### About
GCTA (Genome-wide Complex Trait Analysis) was originally designed to estimate the proportion of phenotypic variance explained by genome- or chromosome-wide SNPs for complex traits (the GREML method), and has subsequently extended for many other analyses to better understand the genetic architecture of complex traits. GCTA currently supports the analyses as follows.
* GRM: estimating the genetic relationships among individuals in GWAS data;
* Estimating the inbreeding coefficients of individuals in GWAS data;
* GREML: estimating the proportion of variance in a phenotype explained by all GWAS SNPs (i.e. the SNP-based heritability);
* Partitioning the genetic variance onto individual chromosomes, MAF bins or functional categories;
* Estimating the genetic variance attributed to chromosome X, and testing for the effect of dosage compensation;
* GREMLd: estimating the dominance variance in unrelated individuals using GWAS data;
* Bivariate GREML: estimating the genetic correlation between two traits (diseases) using GWAS data;
* PCA analysis and estimation of Fst in GWAS data
* Computing LD scores and searching for LD friends for a list of target SNPs;
* Simulating a phenotype based on GWAS data;
* Conditional & joint (COJO) analysis of GWAS summary statistics without individual level genotype data;
* Mixed linear model association analysis;
* Gene- or set-based association analysis;
* Sumamry-data based BLUP;
* Haseman-Elston regression to estimate the the SNP-based heritability for a trait and the genetic correlation between two traits.

**Latest release [v1.90.0beta](#Download), click to download or view update log (8 Aug. 2017)**

### Credits 

[Jian Yang](http://scholar.google.com.au/citations?user=aLuqQs8AAAAJ&hl=en) developed the original version of the software with assistance and guidance from [Hong Lee](http://researchers.uq.edu.au/researcher/2703), [Mike Goddard](mailto:Mike.Goddard@dpi.vic.gov.au) and [Peter Visscher](mailto:peter.visscher@uq.edu.au). [Zhili Zheng](mailto:zhilizheng@outlook.com) rewrote the GRM module, improved the bivariate GREML module, extended the GCTA-PCA module, and designed the website. [Jian Zeng](j.zeng@imb.uq.edu.au) contributed to the GCTA-HEreg module. [Andrew Bakshi](mailto:andrew.bakshi@gmail.com) contributed to the GCTA-fastBAT module. [Zhihong Zhu](mailto:z.zhu1@uq.edu.au) improved the GCTA-COJO module. Robert Marie contributed to the GCTA-sBLUP module. 

### Questions and Help Requests 
If you have any bug reports or questions please send an email to [Jian Yang](http://scholar.google.com.au/citations?user=aLuqQs8AAAAJ&hl=en) at [jian.yang@uq.edu.au](mailto:jian.yang@uq.edu.au)

### Citations 
**Software tool:**  
Yang J, Lee SH, Goddard ME and Visscher PM. GCTA: a tool for Genome-wide Complex Trait Analysis. Am J Hum Genet. 2011 Jan 88(1): 76-82. \[[PubMed ID: 21167468](http://www.ncbi.nlm.nih.gov/pubmed/21167468)\]

**Method for estimating the variance explained by all SNPs (GREML method) with its application in human height:**  
Yang J, Benyamin B, McEvoy BP, Gordon S, Henders AK, Nyholt DR, Madden PA, Heath AC, Martin NG, Montgomery GW, Goddard ME, Visscher PM. Common SNPs explain a large proportion of the heritability for human height. Nat Genet. 2010 Jul 42(7): 565-9. \[[PubMed ID: 20562875](http://www.ncbi.nlm.nih.gov/pubmed/20562875)\]

**GREML method being extended for case-control design with its application to the WTCCC data:**  
Lee SH, Wray NR, Goddard ME and Visscher PM. Estimating Missing Heritability for Disease from Genome-wide Association Studies. Am J Hum Genet. 2011 Mar 88(3): 294-305. \[[PubMed ID: 21376301](http://www.ncbi.nlm.nih.gov/pubmed?term=Estimating%20Missing%20Heritability%20for%20Disease%20from%20Genome-wide%20Association%20Studies)\]

**Extension of GREML method to partition the genetic variance into individual chromosomes and genomic segments with its applications in height, BMI, vWF and QT interval:**  
Yang J, Manolio TA, Pasquale LR, Boerwinkle E, Caporaso N, Cunningham JM, de Andrade M, Feenstra B, Feingold E, Hayes MG, Hill WG, Landi MT, Alonso A, Lettre G, Lin P, Ling H, Lowe W, Mathias RA, Melbye M, Pugh E, Cornelis MC, Weir BS, Goddard ME, Visscher PM: Genome partitioning of genetic variation for complex traits using common SNPs. Nat Genet. 2011 Jun 43(6): 519-525. \[[PubMed ID: 21552263](http://www.ncbi.nlm.nih.gov/pubmed/21552263)\]

**Method for conditional and joint analysis using summary statistics from GWAS with its application to the GIANT meta-analysis data for height and BMI:**  
Yang J, Ferreira T, Morris AP, Medland SE; Genetic Investigation of ANthropometric Traits (GIANT) Consortium; DIAbetes Genetics Replication And Meta-analysis (DIAGRAM) Consortium, Madden PA, Heath AC, Martin NG, Montgomery GW, Weedon MN, Loos RJ, Frayling TM, McCarthy MI, Hirschhorn JN, Goddard ME, Visscher PM (2012) Conditional and joint multiple-SNP analysis of GWAS summary statistics identifies additional variants influencing complex traits. Nat Genet 44(4):369-375. \[[PubMed ID: 22426310](http://www.ncbi.nlm.nih.gov/pubmed/22426310)\]

**Bivariate GREML method:**  
Lee SH, Yang J, Goddard ME, Visscher PM Wray NR (2012) Estimation of pleiotropy between complex diseases using SNP-derived genomic relationships and restricted maximum likelihood. Bioinformatics. 2012 Oct 28(19): 2540-2542. \[[PubMed ID: 22843982](http://www.ncbi.nlm.nih.gov/pubmed/22843982)\]

**Mixed linear model based association analysis:**  
Yang J, Zaitlen NA, Goddard ME, Visscher PM and Price AL (2013) Mixed model association methods: advantages and pitfalls. Nat Genet. 2014 Feb;46(2):100-6. \[[Pubmed ID: 24473328](http://www.ncbi.nlm.nih.gov/pubmed/24473328)\]

**GREML-LDMS method and LD-score calculation:**  
Yang et al. (2015) Genetic variance estimation with imputed variants finds negligible missing heritability for human height and body mass index. Nat Genet, doi: 10.1038/ng.3390.\[[PMID: 26323059](http://www.nature.com/ng/journal/vaop/ncurrent/full/ng.3390.html)\]

**Method to search for LD friends:**  
Yang et al. (2011) Genomic inflation factors under polygenic inheritance. Eur J Hum Genet. 19(7): 807-812. \[[Pubmed ID: 21407268](http://www.nature.com/ejhg/journal/v19/n7/full/ejhg201139a.html)\]

**fastBAT method:**  
Bakshi A., Zhu Z., Vinkhuyzen A.A.E., Hill W.D., McRae A.F., Visscher P.M., and Yang J. (2016). Fast set-based association analysis using summary data from GWAS identifies novel gene loci for human complex traits. Scientific Reports 6, 32894. \[[PMID: 27604177](https://www.nature.com/articles/srep32894)\]

<p style="color: rgb(204,0,0);font-weight:bold;">Last update: 8 Aug 2017</p>

