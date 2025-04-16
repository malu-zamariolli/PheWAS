# Script to combine VIF results

####################
### Packages
#####################
if (!require(data.table)) install.packages("data.table")
library(data.table)

########################
### Inputs and outputs
###########################
# input
input_mult <- snakemake@input[["mult"]]
input_cont <- snakemake@input[["cont"]]
input_binary <- snakemake@input[["binary"]]
input_ordinal <- snakemake@input[["ordinal"]]

# Combine all input files into a list
input_files <- c(input_mult, input_cont, input_binary, input_ordinal)

# output
output <- snakemake@output[[1]]  

# Print output to ensure it's correctly referenced
print(paste("Output file path is:", output))

########################
### Load inputs 
###########################
# Initialize an empty list to store data frames
data_frames <- list()

# Iterate over the input files and read them into data frames
for (file in input_files) {
  print(paste("Processing file:", file))  # Print file for debugging purposes
  # Check if file exists
  if (file.exists(file)) {
    # Read the file into a data frame
    df <- as.data.frame(fread(file))
    # Append the data frame to the list
    data_frames <- append(data_frames, list(df))
  } else {
    warning(paste("File does not exist:", file))
  }
}

########################
### Combine dataframes
###########################
if (length(data_frames) > 0) {
  df.all <- do.call(rbind, data_frames)
  
  ######################
  ### Save
  ###########################
  # Ensure output is a valid single file path
  if (is.character(output) && length(output) == 1 && !is.na(output)) {
    print("Saving combined data frame to output file")
    fwrite(df.all, output, quote = FALSE, col.names = TRUE, 
           row.names = FALSE, sep = "\t", na = "NA")
  } else {
    stop("Output is not a valid single file path. Please check your snakemake output configuration.")
  }
} else {
  stop("No data frames to combine. Please check your input files.")
}
