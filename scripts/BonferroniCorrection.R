#### Bonferroni p.value correction
####### Packages #####################
if(!require(data.table)) install.packages("data.table")
if(!require(dplyr)) install.packages("dplyr")

library(data.table)
library(dplyr)


#################################################
#### Inputs and output ###################
#################################################

# Inputs
codebook <- snakemake@input[["codebook"]] 
mult <- snakemake@input[["mult"]] 
cont <- snakemake@input[["cont"]] 
binary <- snakemake@input[["binary"]] 
ordinal <- snakemake@input[["ordinal"]] 

# Outputs
mult_sig <- snakemake@output[["mult_sig"]]
cont_sig <- snakemake@output[["cont_sig"]]
binary_sig <- snakemake@output[["binary_sig"]]
ordinal_sig <- snakemake@output[["ordinal_sig"]]
summary <- snakemake@output[["summary"]]

#################################################
#### Load inputs ###################
#################################################

codebook <- as.data.frame(fread(codebook))
mult <- as.data.frame(fread(mult))
cont <- as.data.frame(fread(cont))
binary <- as.data.frame(fread(binary))
ordinal <- as.data.frame(fread(ordinal))



#####################################################
### Calculate  number of tests and Bonferroni pvalue
######################################################
# Number of traits
traits <- codebook[codebook$NewVariableName != "ID", ]
n_traits <- nrow(traits)

# Number of snps
n_snps <- snakemake@params[["n_snps"]]

# p-value after correction
pvalue <- 0.05/(n_traits * n_snps)


#####################################################
### Filter significant results
######################################################
mult_sig_df <- mult[mult$p.value < pvalue, ]
cont_sig_df <- cont[cont$p.value < pvalue, ]
binary_sig_df <- binary[binary$p.value < pvalue, ]
ordinal_sig_df <- ordinal[ordinal$p.value < pvalue, ]

#####################################################
### Check number of traits and SNPs that were in sig results
######################################################
mult_trait_snp <- mult_sig_df[, c("term", "Pheno")]
cont_trait_snp <- cont_sig_df[, c("term", "Pheno")]
binary_trait_snp <- binary_sig_df[, c("term", "Pheno")]
ordinal_trait_snp <- ordinal_sig_df[, c("term", "Pheno")]

# Combine data
list_df <- list(c(mult_trait_snp, cont_trait_snp, binary_trait_snp, ordinal_trait_snp))
trait_snp <- do.call(rbind, list_df)

# Select SNPs out of covariates + SNPs
cov <- unique(codebook[, "covariates"])
cov <- unlist(strsplit(cov, ","))

rsids <- grep("^rs", trait_snp$term, value = TRUE)
trait_snp <- trait_snp[trait_snp$term %in% rsids, ]

# Number of unique traits and snps
n.sig <- nrow(trait_snp)

sig_snps <- length(unique(trait_snp$term))
sig_traits <- length(unique(trait_snp$Pheno))


#####################################################
### Select traits and SNPs from significant results
######################################################

#############
## Save descriptive in text format
#################
text <- paste0("Total number of traits: ", n_traits, "\n", 
               "Total number of SNPs: ", n_snps, "\n",
               "p.value after Bonferroni correction 0.05/(N_traits * N_snps): ", pvalue, "\n",
               "Number of significant findings (after Bonferroni): ", n.sig, "\n",
               "Number of SNPs in significant findings: ", sig_snps, "\n",
               "Number of traits in significant findings: ", sig_traits)

writeLines(text, summary)

#####################################################
### Save significant results
######################################################
fwrite(mult_sig_df, mult_sig, quote = FALSE, row.names = FALSE, col.names = TRUE,
       sep = "\t", na = NA)
fwrite(cont_sig_df, cont_sig, quote = FALSE, row.names = FALSE, col.names = TRUE,
       sep = "\t", na = NA)
fwrite(binary_sig_df, binary_sig, quote = FALSE, row.names = FALSE, col.names = TRUE,
       sep = "\t", na = NA)
fwrite(ordinal_sig_df, ordinal_sig, quote = FALSE, row.names = FALSE, col.names = TRUE,
       sep = "\t", na = NA)
