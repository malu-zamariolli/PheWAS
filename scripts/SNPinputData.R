######################### Split SNP data per dataframes for parallelization ##################

if(!require(data.table)) install.packages("data.table")
if(!require(dplyr)) install.packages("dplyr")

library(data.table)
library(dplyr)

######################
#### Inputs and output
##############################

snps <- snakemake@input[["raw"]]


splitSNP <- snakemake@output[["splitSNP"]]
allsnps <- snakemake@output[["allSNPs"]]

########################
### Load input
#######################
df <- as.data.frame(fread(snps))

######################
### Select and split
#############################
# Remove columns that are not of interest
df <- df %>% select(-c(FID, PAT, MAT, SEX, PHENOTYPE))
names(df)[[1]] <- "ID"

# Save
fwrite(df, allsnps, quote = FALSE, col.names = TRUE, 
       row.names = FALSE, sep = "\t", na = NA)


# Loop to split: one dataframe per SNP and save
df.list <- list()
for (i in 2:ncol(df)) {
  df.split <- df[, c(1 ,i)]
  df.list[[i-1]] <- df.split
  fwrite(df.split, file = splitSNP[[i-1]], na = NA, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}







