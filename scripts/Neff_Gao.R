############ Number of effective tests ##########################

####### Packages #####################

if(!require(data.table)) install.packages("data.table")
if(!require(dplyr)) install.packages("dplyr")
if(!require(tidyr)) install.packages("tidyr")

library(data.table)
library(dplyr)
library(tidyr)

#################################################
#### Inputs and output ###################
#################################################

# Inputs
df_input <- snakemake@input[["df_input"]] # traits or snps 
codebook <- snakemake@input[["codebook"]]

# Params
trait_type <- snakemake@params[["trait_type"]]

# Outputs
Neff <- snakemake@output[["Neff"]] # data frame with Neff


#################################################
#### Load inputs ###################
#################################################

df_input <- as.data.frame(fread(df_input))
codebook <- as.data.frame(fread(codebook))

# Switch unordered for ordered
codebook[codebook$VariableType %in% "UNORDERED", "VariableType"] <- "ORDERED"


#################################################
#### Select variables ###################
#################################################

if (trait_type == "SNP") {
  df <- df_input 
} else {
  traits <- codebook[codebook$VariableType %in% trait_type, "NewVariableName"]
  df <- df_input %>% dplyr::select(all_of(c("ID", traits)))
}


#################################################
#### Prepare matrix ###################
#################################################

# Individuals into rownames and exclude ID column
rownames(df) <- df$ID
df <- as.data.frame(df[, -1])


# Make variables into numeric
df <- mutate_all(df, function(x) as.numeric(as.character(x)))

# Tranform in matrix
df_matrix <- as.matrix(df)

print(paste0("Final dimensions: ", nrow(df), " (samples; row) x ", ncol(df), " (variables; column)"))


#################################################
#### Calculate correlation ###################
#################################################

# Calculate the correlation matrix
df_matrix_cor <- cor(df_matrix, use = 'pairwise.complete.obs')

# Set missing values to 0 to allow svd()
df_matrix_cor[which(is.na(df_matrix_cor))] <- 0


#################################################
#### Calculate eigenvalues ###################
#################################################

# svd() computes the singular-value decomposition of a rectangular matrix, with $d being a vector containing the singular values of the decomposition
df_matrix_EV <- svd(df_matrix_cor)$d 


#################################################
#### Calculate Neff ###################
#################################################
# Neff is defined as the number of eigenvalues required to explain 99.5% of the variation from the matrix data 
sum_EV <- 0
count <- 0

while(sum_EV/sum(df_matrix_EV) < 0.995) {
  count <- count + 1
  sum_EV <- sum_EV + df_matrix_EV[count]
}

print(paste0("Neff for ", ncol(df), " variables: " , count))

#################################################
#### Store Neff in dataframe ###################
#################################################

output <- data.frame(Type = trait_type, 
                     Neff = count)

fwrite(output, Neff, quote = FALSE,
       row.names = FALSE, col.names = TRUE, 
       sep = "\t", na = NA)