# PheWAS (Phenome-Wide Association Study)
This pipeline was used in the following publication: 

Moysés-Oliveira et al 2025. Pleiotropic effects of APOE variants on a sleep-based adult epidemiological cohort. Sleep Medicine.
https://doi.org/10.1016/j.sleep.2025.106490


# Usage

This Github contains the workflow pipeline that was used to run PheWAS. The workflow is organized as a [snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline. 

To run the pipeline with snakemake:
```console
cd  PheWAS/
mamba activate snakemake
snakemake --cores 7 
# number of cores can be modified accordingly
```

# Input data
## Plink Files
Path to .bed, .bim, .fam files from plink

## Codebook: 
.txt file with the following essential columns (other columns may be present):

|NewVariableName|VariableType|covariates|
| -------- | --- | --- | 
|BMI|CONTINUOUS|PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,age,age2,sex,cce|

VariableType should be CONTINUOUS, BINARY, ORDERED, UNORDERED

## Traits:
.txt file with IDs and one phenotype per column. Continuous traits should be normalized when needed.

## Traits_raw
.txt file with IDs and one phenotype per column. Continuous traits before normalization. This file is used for non-parametric testing.

## SNV file 
.snps file with SNV ID to be extracted from plink files.

# Config file
*config_PheWAS.yaml* with run name, input data and other parameters:
- chunks: number of traits that will be run at once. Full data is divided in chunks in script TraitByType.R
- Number of chunks per trait type: binary, ordinal, continuous and unordered. Example: if 180 continuous traits, n_cont will be 4 (if chunks = 50).
- snp_model: can be A (additive) or D (dominant). If A, plink will generate input file (.raw) and Non_parametric_A.R will be applied. If D, input file needs to be done in preprocessing and Non_parametric_D.R will be applied.
- n_snps: number of SNVs in the PheWAS
  
# Covariates script
*scripts/covariates/{wildcards.phewas}_covariates.R*

- Each run needs a new script to create a .txt file with ID and all covariates to be used.

# Tools
- [Plink.v1.9](https://www.cog-genomics.org/plink/):  *../envs/plink.yaml*
- R 
- qpdf v.11.9.1 *../envs/qpdf.yaml*
