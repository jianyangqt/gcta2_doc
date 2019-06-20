
## Download
### Executable Files {: .notoc}

The executable files below only support a 64-bit operating system on the x86\_64 CPU platform. 

Linux [gcta\_1.92.2beta.zip](./gcta_1.92.2beta.zip)

Windows [gcta\_1.92.2beta\_win.zip](./gcta_1.92.2beta_win.zip)

Mac [gcta\_1.92.2beta\_mac.zip](./gcta_1.92.2beta_mac.zip)
 
The executable files are released under the MIT license. We recommend to use the Linux version because the Windows and Mac versions have not been fully tested.

> Note: GCTA 1.92.2beta is a beta version under testing. We have fixed a few bugs in the previous versions. If you find any bug in this version, please report it to Jian Yang at [jian.yang@uq.edu.au](mailto:jian.yang@uq.edu.au).

### Source code {: .notoc}

[gcta\_1.26.0\_src.zip](./gcta_1.26.0_src.zip)

The source code are released under GPL v3. The source code of the latest version will be released when it is stable.

### Update log {: .notoc}
#### Version 1.92.2beta (18 June 2019)
* Modified mtCOJO to accept LD score files with 4 columns. 
* Changed flag --gsmr-beta to --gsmr2-beta.
* Fixed a bug in GWAS simulation.
* Fixed a bug in fastGWA if a p-value is extremely small.
* Improved the performance of fastGWA.

#### Version 1.92.1beta6 (13 Apr 2019)
* Fixed a bug in --reml-bivar due to the update of Linux compiler.

#### Version 1.92.1beta5 (1 Apr 2019)
* Added a flag '[--reml-res-diag](#GREMLanalysis)' to specify the diagonal elements of the residual correlation matrix in REML.
* Added a new module [fastGWA](#fastGWA) (an extremely resource-efficient tool for mixed linear model association analysis of biobank-scale GWAS data).

#### Version 1.92.0beta3 (1 Feb 2019)
* Fixed a bug in COJO for some circumstances where the standard errors of SNP effects are extremely small.

#### Version 1.92.0beta (22 Jan 2019)
* Added a flag '--gsmr-beta' to use a testing version of the HEIDI-outlier method.
* Modified '--make-grm-alg' so that it can be used in combination with –make-grm-part. 

#### Version 1.91.7beta (8 Oct 2018)
* Added a flag (--mtcojo-bxy) in the mtCOJO analysis to read the effects of covariates on trait from a user-specified file.
* Added a multi-SNP-based HEIDI-outlier test in the HEIDI-outlier analysis.
* Added a function in mtCOJO to compute the effects of covariates on trait from a genetic correlation analysis if there are not enough SNPs to perform the GSMR analysis.
* Modified mtCOJO and GSMR to read  summary data from compressed text files.
* Modified HEIDI-outlier to save the removed pleiotropic SNPs in text file.
* Modified GSMR to sort SNPs by chi-squared values in the clumping analysis.
* Fixed a bug in the GSMR effect plot.
* Fixed a bug in --reml-bivar.
* Changed the flag --mlma-no-adj-covar to --mlma-no-preadj-covar.

#### Version 1.91.6beta (17 Aug 2018)
* Changed the criterion of selecting the top associated SNP by p-value in COJO to that by chi-squared value to avoid the problem of having extremely small p-values (e.g. those = 1e-300). 
* Fixed a bug in COJO in the Windows version.
* Fixed an issue related to allele frequency in COJO when the first allele differs between the GWAS summary data and the LD reference sample.
* Added a function to check the consistency of allele frequency between the GWAS summary data and the reference sample in COJO (the --diff-freq flag).
* Fixed a bug when manipulating the GRM in PCA.
* Added a new flag --unify-grm to unify the order of the IDs in multiple GRM files.

#### Version 1.91.5beta (7 Jul 2018)
* Fixed a bug in GSMR when there is a very small number of SNPs used to run an HEIDI-outlier analysis.
* Added a flag --diff-freq to check difference in allele frequency between data sets in the GSMR and mtCOJO analyses.
* Removed flags --clump-p1 and --heidi-snp from the GSMR and mtCOJO analyses
* The flag --gsmr-snp has been superseded by --gsmr-snp-min.
* Improved compatibility with old Linux version.
* Added a flag --mbfile to merge multiple BED files (e.g. genotype data of each chromsome saved in a separate BED file) into a single BED file.
* Fixed a memory issue with the flag --make-grm-x.
* Fixed a build stack issue in the Windows version.
* Fixed a rare thread freezing with --make-grm.

#### Version 1.91.4beta (17 Apr 2018)
* Fixed a bug in GSMR when there are multiple outcome variables.
* Fixed a bug in COJO when the standard error is extremely small.
* Improved the speed and memory usage of --make-grm-xchr, and added an option --make-grm-xchr-part to reduce the memory usage further.
* Added --mbfile in GRM functions to proceed genotypes stored in multiple PLINK files.
* Updated the options --update-sex, --update-ref-allele and --update-freq to be compatible with the new GRM functions.
* Fixed a bug of reporting "Illegal instruction" error for old CPUs (earlier than 2009). 
* Added an additional option --threads to specify the number of threads (the same as --thread-num). The number of threads will be obtained from standard OpenMP environment variable OMP\_NUM\_THREADS if --thread-num or --threads is not specified. 

#### Version 1.91.3beta (14 Mar 2018)
* Speeded up [dominance GRM](#GREMLfordominancevariance) and added a flag --make-grm-d-part to partition the computation.
* Fixed a bug in REML, REML bivar, MLMA and LD when the number of threads (specified by --thread-num) is larger than 1.
* Redirected the log output to both screen and .log file.
* Fixed a bug in COJO for the X chromosome when there is no gender information in the .fam file.
* Fixed a bug in mtCOJO.
* Added a flag (--effect-plot) in [GSMR](#GSMR) for visualization.

#### Version 1.91.2beta (2 Feb 2018)
* Added a new module [GSMR](#GSMR).
* Added a flag (--mbfile) to read multiple PLINK binary files for [GSMR](#GSMR) and [mtCOJO](#mtCOJO).
* Fixed a bug in SBLUP, and improved the speed by 40%.
* Fixed a bug in MLMA when dealing with the individuals' ID.
* Fixed unreadable characters in the output of some computer clusters.

#### Version 1.91.1beta (25 Nov 2017)
* Fixed a bug in --mtcojo.
* Fixed a memory issue in REML analysis and improved the speed by 3 times in the Linux version.
* Changed to use the shared library glibc avoid segmentation fault in higher versions of Linux kernel.

#### Version 1.91.0beta (21 Oct 2017)
* Added a new module [mtCOJO](#mtCOJO)
* Fixed an issue of file path in the Windows version

#### Version 1.90.2beta (24 Sep 2017)
* Fixed a bug in --mlma-loco with the --mlma-no-adj-covar option.
* Fixed a bug in --make-grm-part when the sample size of one partition is larger than 69K.
* Fixed the performance issue in reading the PLINK .fam file.
* Fixed an issue with --autosome-num.
* Removed the VC++ runtime dependency in the Windows version. 

#### Version 1.90.1beta (13 Sep 2017)
* Fixed a bug in estimating allele frequency in some occasions.
* Fixed a bug in computing a GRM occasionally in small sample.
* Fixed an issue in computing a GRM including rare variants.
* Fixed an issue to run Linux binary in the Linux subsystem on Windows 10.
* Fixed a memory issue in the Windows version.
* Removed --grm-no-relative and added --grm-singleton to get singleton subjects from a sample.
* Fixed a memory issue in --make-bK.

#### Version 1.90.0beta (8 Aug 2017)
* Improved the speed and memory usage of GRM computation by orders of magnitude.
* Added a new option --make-grm-part to partition the GRM computation into a large of parts to facilitate the analysis in large data set such as the UK Biobank.
* Improved the memory usage of the --grm-cutoff option.
* Added the --grm-no-relative option to extract the GRM of a subset of individuals who do not have any close relative in the sample.
* Improved the speed and memory usage of --freq by orders of magnitude.
* Improved the approximation accuracy of the COJO analysis.
* Added an option --cojo-sblup to perform a summary-data-based BLUP prediction analysis.
* Added the Haseman-Elston regression analysis to estimate the SNP-based heritability for a trait and genetic correlation between traits.
* Improved the speed of the bivariate GREML analysis (5X faster than original version).
* Added the Mac and Windows versions.
* Update the package dependencies to the latest, such as Intel MKL and Eigen. This improved the performance by ~40%.
* Fixed the memory issue when the sample size exceeds 500K in some functions (e.g. bivariate GREML and reading the GRM in gz format).

#### Version 1.26.0 (22 June 2016)
Download link: [gcta_1.26.0.zip](./gcta_1.26.0.zip)

* Fixed a bug in MLMA.
* Added a new module (GCTA-fastBAT) for a set- or gene-based association analysis using GWAS summary data.

#### Version 1.25.3 (27 April 2016)

Download link: [gcta_1.25.3.zip](./gcta_1.25.3.zip)

* Fixed a memory leaking issue in --mlma

#### Version 1.25.2 (22 Dec 2015)
Download link: [gcta_1.25.2.zip](./gcta_1.25.2.zip)

* A much more memory-efficient version of MLMA.
* Added a new option (--mlma-subtract-grm) for MLMA-LOCO with large data sets.
* Fst calculation has been changed to that based on a random model. The previous version was based on a fixed model. The difference is trivial for small Fst values but the random model has a good property that Fst is bounded at 1 for the most extreme allele frequency difference.
* Added a new option (--make-grm-inbred) to compute GRM for an inbred population (e.g. inbred mice or crops).
* Added a new option (--recode-std) to output standardised SNP genotypes.

#### Version 1.25.1 (8 Dec 2015)
Download link: [gcta_1.25.1.zip](./gcta_1.25.1.zip)

* Added an option <a href="#EstimatevarianceexplainedbyalltheSNPs" target="_blank">--reml-bendV</a>

#### Version 1.25.0 (30 Oct 2015)
Download link: [gcta_1.25.0.zip](./gcta_1.25.0.zip)

* Fixed a bug in --imp-rsq
* Added an option to calculate an unbiased estimate of LD score for LDSC regression analysis (see gcta.freeforums.net/thread/177/gcta-lds-calculating-score-snp); Added an option to calculate multi-component LD score following Finucane et al. (2015 Nat Genet).
* Added options to <a href="#Datamanagement" target="_blank">extract or exclude a region</a>.
* Add the --reml-bivar-no-constrain option to the <a href="#BivariateREMLanalysis" target="blank">bivariate GREML analysis</a>.
* Add an <a href="#GCTA-COJO:Conditional&jointassociationanalysis" target="_blank">option</a> to select a fixed number of top associated SNPs (taking LD into account) from GWAS.
* We have implemented the Zaitlen et al. method in GCTA which allows to estimate SNP-based h2 in family data without having to remove related individuals.

#### Version 1.24.7 (11 June 2015)
Download link: [gcta_1.24.7.zip](./gcta_1.24.7.zip)

* Mixed linear model association (MLMA) analysis with multiple GRMs
* Fst calculation
* Haseman-Elston regression
* LD score calculation
#### Version 1.24.4 (29 July 2014)
* changed the syntax for the conditional and joint analysis; fixed memory leak issues in mixed linear model based association analysis and bivariate GREML analysis with multiple GRMs; enabled the function converting dosage data to PLINK best guess.

#### Version 1.24.3 (5 Jun 2014)
* allows you to transform variance explained by all SNPs on the observed scale to that on the underlying scale in a bivariate analysis of a case-control study and a quantitative trait; pca
* only the top eigenvalues will be printed out.
#### GCTA-GREML Power Calculator (11 Apr 2014).
#### Version 1.24.2 (12 Mar 2014)
* fixed a bug in the conditonal and joint analysis (GCTA-COJO) when doing a backward model selection.
#### Version 1.24.1 (6 Mar 2014)
* a small change that allows you to use "Rsq" or "Rsq_hat" as the header for the last column of the *.mlinfo file from MACH imputation.
#### Version 1.24 (8 Jan 2014)
* fixed a bug in REML analysis as a result of a change made in v1.23 in transforming the estimate of genetic variance on the observed scale to that on the underlying scale; fixed a bug in GWAS simulation where the reported variance explained by a causal variant in the *.par file was incorrect.
#### Version 1.23 (18 Dec 2013)
* changed  --dosage-mach option and added a new option --dosage-mach-gz; fixed a bug in the --cojo-cond option when two SNPs are in very high LD and their allele frequencies are consistently higher in the reference sample than those in the discovery sample.
#### Version 1.22 (31 Oct 2013)
* fixed a bug in the --dosage-mach option when used in combined with the --imput-rsq option.
#### Version 1.21 (16 Oct 2013)
* fixed a bug in bivariate analysis including covariates; re-wrote the code for the option --dosage-mach; added a new option and changed syntax for the mixed linear model association analysis.
#### Version 1.20 (23 Aug 2013)
* added a new module mixed linear model association analysis; fixed a few bugs; made a few improvements.
#### Version 1.13 (19 Mar 2013)
* fixed a bug for the --make-grm-bin option.
#### Version 1.11 (14 Feb 2013)
* fixed a bug for the --mgrm-bin option and added the option to test for genetic correlation = 0 or 1 in a bivariate analysis.
#### Version 1.1 (10 Feb 2013)
* a much faster version which allows multi-thread computing (new option --thread-num); added new options --make-grm-bin and --grm-bin to more efficiently read and write the GRM files.
#### Version 1.04 (13 Sep 2012)
* added a new option to convert Minimac dosage data to PLINK binary PED format.
#### Version 1.03 (30 Aug 2012)
* fixed a few bugs and added a new option to convert MACH dosage data to PLINK binary PED format.
#### 29 July 2012
* fixed 2 bugs.
#### 16 July 2012
* fixed a few bugs.
#### 14 May 2012
* version 1.0 released!
#### 30 Nov 2011
* latest version (version 0.93.9) of source codes released.
#### Version 0.93.9 (18 Nov 2011)
* modified the --dosage-mach option to be compatiable with the latest MACH version; fixed a bug with the option --ld.
#### Version 0.93.8 (30 Sep 2011)
* fixed a bug for the option --grm-adj when the genotype data of some individuals are completely missing.
#### Version 0.93.7 (10 Sep 2011)
* fixed a bug when the option --ibc is used in combined with the option --keep or --remove, which causes wrong IDs in the *.ibc fie; fixed a bug in --gxe option when there are missing values for the environmental factor; and modified the function for converting Illumina raw genotype data to that in PLINK format.
#### Version 0.93.6 (28 Aug 2011)
* fixed a bug in the new option --reml-lrt which caused memory leak.
#### Version 0.93.5 (26 Aug 2011)
* added an option to turn off the LRT and fixed a bug in the case that the IDs of multiple GRM files are not in the same order.
#### Version 0.93.4 (15 Aug 2011)
* added a function to calculate the LRT for the REML analysis.
#### Version 0.93.2 (18 Jul 2011)
* fixed a bug in the matrix bending subroutine.
#### Version 0.93.1 (12 Jul 2011)
* improved the efficiency of reading PLINK binary data. 
#### Version 0.93.0 (8 Jul 2011)
* added a subroutine to deal with the issue when the variance-covariance matrix V is negative-definite; changed the default number of maximum REML iterations from 30 to 100; changed the method of calculating the diagonal elements of GRM to be the same as that for the off-diagonal elements; modified REML procedure to allow some elements of the GRM to be missing (printing a warning on the screen in stead of an error message).   
#### 8 Apr, 2011
* fixed a bug in GWAS simulation.
#### 2 Apr, 2011
* fixed a bug in a REML analysis, i.e. the estimate may be stuck at zero if the true parameter is very small.
#### 24 Mar, 2011
* modified the output of LD estimation and the input format of GWAS simulation
#### 10 Feb, 2011
* fixed a few bugs.
#### 24 Dec, 2010
* added a few new functions, e.g. convert the raw genotype data into PLINK binary format.
#### 23 Nov, 2010
* source codes released.
#### 14 Oct, 2010
* fixed a bug in reading the PLINK FAM file.
#### 13 Oct, 2010
* MacOS version released.
#### 11 Oct, 2010
* fixed a bug in transforming the estimate of variance explained by the SNPs on the observed scale to that on the underlying scale for a case-control study.
#### 17 Sep, 2010
* fixed a bug in the estimation of LD and compiled the program statically (more compatible 
#### 30 Aug, 2010
* first release.

