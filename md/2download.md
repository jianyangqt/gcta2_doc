
## Download
### Executable Files {: .notoc}

[gcta_1.26.0.zip](./gcta_1.26.0.zip)

The executable files (binary code) are release under MIT lincense.
Unfortunately, we have stopped updating the Windows and Mac versions. An earlier version for Windows (gcta.exe) and Mac OS (gcta\_mac) can be found at [gcta\_1.02.zip](./gcta_1.02.zip).

### Source code {: .notoc}

[gcta\_1.26.0\_src.zip](./gcta_1.26.0_src.zip)

The source code are released under GPL v3. 

### Update log {: .notoc}
#### Version 1.26.0 (22 June 2016)

Download link: [gcta_1.26.0.zip](./gcta_1.26.0.zip)

* Fixed a bug in MLMA.
* Added a new module (GCTA-fastBAT) for a set- or gene-based association analysis using GWAS summary data.
* Released the latest version of the GCTA source code. Download link: [gcta\_1.26.0\_src.zip](./gcta_1.26.0_src.zip)

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
