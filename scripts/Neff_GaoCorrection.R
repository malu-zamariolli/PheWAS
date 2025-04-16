############ Filter sig results by Gao et all ##########################

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
mult <- snakemake@input[["mult"]] 
cont <- snakemake@input[["cont"]] 
binary <- snakemake@input[["binary"]] 
ordinal <- snakemake@input[["ordinal"]] 

neff_bin <- snakemake@input[["neff_bin"]]
neff_cont <- snakemake@input[["neff_cont"]]
neff_cat <- snakemake@input[["neff_cat"]]
neff_snp <- snakemake@input[["neff_snp"]]

# Outputs
mult_sig <- snakemake@output[["mult_sig"]]
cont_sig <- snakemake@output[["cont_sig"]]
binary_sig <- snakemake@output[["binary_sig"]]
ordinal_sig <- snakemake@output[["ordinal_sig"]]
summary <- snakemake@output[["summary"]]

#################################################
#### Load inputs ###################
#################################################

mult <- as.data.frame(fread(mult))
cont <- as.data.frame(fread(cont))
binary <- as.data.frame(fread(binary))
ordinal <- as.data.frame(fread(ordinal))

neff_bin <- as.data.frame(fread(neff_bin))
neff_cont <- as.data.frame(fread(neff_cont))
neff_cat <- as.data.frame(fread(neff_cat))
neff_snp <- as.data.frame(fread(neff_snp))



#####################################################
### Calculate corrected p-value
######################################################
# Number of effective traits
n_bin <- neff_bin$Neff
n_cont <- neff_cont$Neff
n_cat <- neff_cat$Neff

n_traits <- n_bin + n_cont + n_cat

# Number of snps
n_snps <- neff_snp$Neff

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

# Select SNPs out of covariates and SNPs
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
               "Neff binary: ", n_bin, "\n",
               "Neff continuous: ", n_cont, "\n",
               "Neff categorical (ordinal + non-ordinal): ", n_cat, "\n",
               "Neff SNPs: ", n_snps, "\n",
               "p.value after correction 0.05/(N_traits * N_snps): ", pvalue, "\n",
               "Number of significant findings (after correction): ", n.sig, "\n",
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
