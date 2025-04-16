##################### Split traits by types and chunks #######################################

if(!require(data.table)) install.packages("data.table")
if(!require(dplyr)) install.packages("dplyr")

#############################
## Inputs  #####
##############################
traits <- snakemake@input[["traits"]]
codebook <- snakemake@input[["codebook"]]

codebook <- fread(codebook)
traits <- fread(traits)

##################################
# Function to split columns into chunks
split_columns_into_chunks <- function(columns, chunk_size) {
  split(columns, ceiling(seq_along(columns) / chunk_size))
}

# Determine the types of traits from codebook
continuous_traits <- codebook[VariableType %in% "CONTINUOUS", NewVariableName]
binary_traits <- codebook[VariableType %in% "BINARY", NewVariableName]
unordered_traits <- codebook[VariableType %in% "UNORDERED", NewVariableName]
ordinal_traits <- codebook[VariableType %in% "ORDERED", NewVariableName]

# Filter data by trait types, retaining ID column
continuous_data <- traits[, c("ID", ..continuous_traits)]
binary_data <- traits[, c("ID", ..binary_traits)]
unordered_data <- traits[, c("ID", ..unordered_traits)]
ordinal_data <- traits[, c("ID", ..ordinal_traits)]

# Define chunk size
chunk_size <- snakemake@params[['chunks']]

# Split columns into chunks: -1 to skip ID
continuous_chunks <- split_columns_into_chunks(names(continuous_data)[-1], chunk_size)
binary_chunks <- split_columns_into_chunks(names(binary_data)[-1], chunk_size)
unordered_chunks <- split_columns_into_chunks(names(unordered_data)[-1], chunk_size)
ordinal_chunks <- split_columns_into_chunks(names(ordinal_data)[-1], chunk_size)

# Function: Write chunks to files, including the ID column
write_chunks <- function(data, chunks, output_paths) {
  for (i in seq_along(chunks)) {
    chunk_data <- data[, c("ID", chunks[[i]]), with = FALSE]
    fwrite(chunk_data, file = output_paths[i], na = NA, sep = "\t")
  }
}


# Write the chunks to separate files
write_chunks(continuous_data, continuous_chunks, snakemake@output[['cont']])
write_chunks(binary_data, binary_chunks, snakemake@output[['binary']])
write_chunks(unordered_data, unordered_chunks, snakemake@output[['unor']])
write_chunks(ordinal_data, ordinal_chunks, snakemake@output[['ordi']])