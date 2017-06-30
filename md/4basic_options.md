
## Basic options

### Input and output

--bfile test  
Input PLINK binary PED files, e.g. test.fam, test.bim and test.bed (see PLINK user manual for details).

--dosage-mach test.mldose test.mlinfo  
Input files in MACH output format (uncompressed), e.g. test.mldose and test.mlinfo (see MACH user manual for details).

--dosage-mach-gz test.mldose.gz test.mlinfo.gz  
Input files in MACH output format (compressed), e.g. test.mldose.gz and test.mlinfo.gz.

> Formats of the input files
> test.mldose
```nohighlight
001->0011 ML_DOSE 2.000 0.000 0.000 0.000 2.000 0.001 0.028 0.017 1.992 0.027  
002->0021 ML_DOSE 2.000 1.000 1.000 1.000 1.999 1.001 1.280 1.010 1.985 1.028  
003->0031-000 ML_DOSE 1.036 1.132 1.000 2.000 1.003 1.999 0.986 1.013 1.030 1.984  
...  
```

> test.mlinfo
```nohighlight
SNP	Al1	Al2	Freq1	MAF	Quality	Rsq  
rs1	G	T	0.8633	0.1367	0.9595	0.8697  
rs2	C	T	0.4654	0.4654	0.9702	0.9543  
rs3	G	T	0.4459	0.4459	0.9997	0.9995  
...  
```

**Note**: the --dosage-mach option was designed to read output files from an early version of MACH, which might not be compatible with output files from the latest version of MACH or Minimac.

--out test  
Specify output root filename

### Data management
--keep test.indi.list  
Specify a list of individuals to be included in the analysis.

--remove test.indi.list  
Specify a list of individuals to be excluded from the analysis.

--chr 1  
Include SNPs on a specific chromosome in the analysis, e.g. chromosome 1.

--autosome-num 22  
Specify the number of autosomes for a species other than human. For example, if you specify the number of autosomes to be 19, then chromosomes 1 to 19 will be recognized as autosomes and chromosome 20 will be recognized as the X chromosome. The default number is 22 if this option not specified. 

--autosome  
Include SNPs on all of the autosomes in the analysis. Note: this option will be overided by the --chr chr_num option, if you want to include all autosomes, please remove the --chr option. 

--extract test.snplist  
Specify a list of SNPs to be included in the analysis.  
> Input file format  
> test.snplist  
```nohighlight
rs103645  
rs175292  
......  
```

--exclude test.snplist  
Specify a list of SNPs to be excluded from the analysis.

--extract-snp rs123678  
Specify a SNP to be included in the analysis.

--exclude-snp rs123678  
Specify a single SNP to be excluded from the analysis.

--extract-region-snp rs123678 1000  
Extract a region centred around a specified SNP, e.g. +-1000Kb region centred around rs123678. 

--exclude-region-snp rs123678 1000  
Exclude a region centred around a specified SNP, e.g. +-1000Kb region centred around rs123678. 

--extract-region-snp 1 120000 1000  
Extract a region centred around a specified bp, e.g. +-1000Kb region centred around 120,000bp of chr 1. 

--exclude-region-snp 1 120000 1000  
Exclude a region centred around a specified bp, e.g. +-1000Kb region centred around 120,000bp of chr 1. This option is particularly useful for a analysis excluding the MHC region.

--maf 0.01  
Exclude SNPs with minor allele frequency (MAF) less than a specified value, e.g. 0.01.

--max-maf 0.1  
Include SNPs with MAF less than a specified value, e.g. 0.1.

--update-sex test.indi.sex.list  
Update sex information of the individuals from a file.  
> Input file format  
> test.indi.sex.list (no header line; columns are family ID, individual ID and sex). Sex coding: "1" or "M" for male and "2" or "F" for female.  
```nohighlight
011 0101 1  
012 0102 2  
013 0103 1  
......  
```

--update-ref-allele test\_reference\_allele.txt  
Assign a list of alleles to be the reference alleles for the SNPs included in the analysis. By default, the first allele listed in the *.bim file (the 5th coloumn) or *.mlinfo.gz file (the 2nd conlumn) is assigned to be the reference allele. NOTE: This option is invalid for the imputed dosage data only.  
> Input file format  
> test\_reference\_allele.txt (no header line; columns are SNP ID and reference allele)  
```nohighlight
rs103645 A  
rs175292 G  
......  
```

--imput-rsq 0.3  
Include SNPs with imputation R<sup>2</sup> (squared correlation between imputed and true genotypes) larger than a specified value, e.g. 0.3.

--update-imput-rsq test.imput.rsq  
Update imputation R<sup>2</sup> from a file. For the imputed dosage data, you do not have to use this option because GCTA can read the imputation R<sup>2</sup> from the *.mlinfo.gz file unless you want to write them. For the best guess data (usually in PLINK format), if you want to use a R<sup>2</sup> cut-off to filter SNPs, you need to use this option to read the imputation R<sup>2</sup> values from the specified file.  
> Input file format  
> test.imput.rsq (no header line; columns are SNP ID and imputation R<sup>2</sup>)  
```nohighlight
rs103645 0.976  
rs175292 1.000  
......  
```

--freq  
Output allele frequencies of the SNPs included in the analysis (in plain text format), e.g.  
> Output file format  
> test.freq (no header line; columns are SNP ID, reference allele and its frequency)  
```nohighlight
rs103645 A 0.312  
rs175292 G 0.602  
......  
```

--update-freq test.freq  
Update allele frequencies of the SNPs from a file rather than calculating from the data. The format of the input file is the same as the output format for the option --freq.

--recode   
Output SNP genotypes based on additive model (i.e. x coded as 0, 1 or 2) in compressed text format, e.g. test.xmat.gz.

--recode-nomiss  
Output SNP genotypes based on additive model without missing data. Missing genotypes are replaced by their expected values i.e. 2p where p is the frequency of the coded allele (also called the reference allele) of a SNP.  

--recode-std  
Output standardised SNP genotypes without missing data. The standardised genotype is *w = (x - 2p) / sqrt[2p(1-p)]*. Missing genotypes are replaced by zero.
> Output file format  
> test.xmat.gz (The first line contains family ID, individual ID and SNP ID. The second line contains two nonsense words "Reference Allele" and the reference alleles of the SNPs. Missing genotype is represented by "NA").  
```nohighlight
FID IID rs103645 rs175292  
Reference Allele A G  
011 0101 1 0  
012 0102 2 NA  
013 0103 0 1  
......
```

--make-bed  
Save the genotype data in PLINK binary PED files (*.fam, *.bim and *.bed).  
Example
```bash
# Convert MACH dosage data to PLINK binary PED format
gcta64  --dosage-mach  test.mldose.gz  test.mlinfo.gz  --make-bed --out test
```
**Note**: the --dosage-mach option was designed to read output files from an early version of MACH, which might not be compatible with output files from the latest version of MACH or Minimac.

### Multi-thread computing

We have made most of the analyses in GCTA being able to run on multiple CPUs.
 
--thread-num   10  
Specify the number of threads on which the program will be running.
 
> Examples
```bash
gcta64 --bfile test --make-grm --out test --thread-num 10
gcta64 --reml --grm test --pheno test.pheno --out test --thread-num 10
```
