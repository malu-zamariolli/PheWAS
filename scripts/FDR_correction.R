######### False Discovery Rate: Benjamini Hochberg ##################

#### Packages ############
shh <- suppressPackageStartupMessages # https://stackoverflow.com/questions/18931006/how-to-suppress-warning-messages-when-loading-a-library

if(!require(data.table)) install.packages("data.table")
if(!require(dplyr)) install.packages("dplyr")
if(!require(broom)) install.packages("broom")


shh(library(data.table))
shh(library(dplyr))
shh(library(broom))

###########################
### Inputs and outpus ####
############################
# Inputs
mult <- snakemake@input[["mult"]] 
cont <- snakemake@input[["cont"]] 
bin <- snakemake@input[["binary"]] 
ord <- snakemake@input[["ordinal"]] 


# Outputs 
fdr <- snakemake@output[["fdr"]]
fdr_sig <- snakemake@output[["fdr_sig"]]

summary <- snakemake@output[["summary"]]

########################################
### Load and prepare input - association
#########################################

uno <- as.data.frame(fread(mult))
cont <- as.data.frame(fread(cont))
bin <- as.data.frame(fread(bin))
ord <- as.data.frame(fread(ord))

# Retain only results from models that converged (not applicable to cont)
uno_ok <- uno[uno$convergence %in% "YES", ]
ord_ok <- ord[ord$convergence %in% c("YES", 0), ]
bin_ok <- bin[bin$Convergence == TRUE, ]

# Combine all results in list
df_list <- list(uno_ok, ord_ok, bin_ok, cont)
names(df_list) <- c("Mult", "Ordi", "Log", "Linear")

#################################
### Select SNP associations
##################################
# Retain only associations with SNP, not covariates
rs_list <- list()
# loop to select SNPs 
for (i in 1:length(df_list)) {
  df <- as.data.frame(df_list[[i]])
  index <- grep("^rs", df[["term"]])
  df_rs <- df[index, ]
  df_rs <- df_rs %>%  select("term", "estimate", "std.error", "statistic", "p.value", 
                             "conf.low", "conf.high", "Pheno", "N", "method")
  rs_list[[i]] <- df_rs
}
names(rs_list) <- c("Mult", "Ordi", "Log", "Linear")

# Combine all dfs
df.pvalue <- do.call(rbind, rs_list)


##############################################
### Select multinomial representative pvalue
##############################################

# Filter multinomial regression results
multinomial_results <- df.pvalue %>%
  filter(method == "MultiLogReg")

# Group by phenotype and SNP and select a representative p-value (Minimum p-value)
representative_p_values <- multinomial_results %>%
  group_by(Pheno, term) %>%
  filter(p.value == min(p.value)) %>%
  ungroup()

# Combine with all p.values
p.value_sans_mult <- df.pvalue %>% filter(method != "MultiLogReg")
df.pvalue_representative <- rbind(p.value_sans_mult, representative_p_values)

##############################################
### Calculate FDR correction
##############################################
df.pvalue_representative$p.value.adj <- p.adjust(df.pvalue_representative$p.value, method = "BH")

# Significance by FDR
significant_fdr <- df.pvalue_representative[df.pvalue_representative$p.value.adj < 0.05, ]

#####################################################
### Create descriptive file
######################################################
n_sig <- nrow(significant_fdr)

#############
## Save descriptive in text format
#################
text <- paste0("Significant after FDR: ", n_sig)

writeLines(text, summary)

##############################################
### Save
##############################################
fwrite(df.pvalue_representative, fdr, quote = FALSE,
    row.names = FALSE, col.names = TRUE, na = NA, sep = "\t")

fwrite(significant_fdr, fdr_sig, quote = FALSE,
    row.names = FALSE, col.names = TRUE, na = NA, sep = "\t")

