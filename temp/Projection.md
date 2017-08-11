### PC loading and projection

--pc-loading ref_pca  
Generate the SNP loading from the GCTA-PCA results of the reference genotype

--project-loading ref\_snp\_loading 20  
Project the target genotype into the reference SNP loading (ref\_snp_loading), and generate 20 principal component

> Note: 
> * The mismatch SNPs in the reference and target genotype will be obsoleted in the PC loadings, so it might be biased if lots of SNPs missing in the target genotype.
> * The default mode uses the MAF in reference genotype to calculate the projection, which coincides with EIGENSOFT projection.

**Example**  
REF: reference genotype, the reference genotype to generate PC loading;  
TAR: target genotype to project the PC to REF;
```bash
# make GRM
gcta64 --bfile REF --maf 0.01 --autosome --make-grm --out REF
# PCA analysis
gcta64 --grm REF --pca 20 --out REF_pca20

# use the pca generated above to produce the SNP loading
gcta64 --bfile REF --pc-loading REF_pca20 --out REF_snp_loading

# project the TAR to the SNP loading. The number 20 is the PC want to project.
# It is the same if do chromosome by chromosome calculation add them up after each calculation 
gcta64 --bfile TAR --project-loading REF_snp_loading 20 --out TAR_pca20
```
