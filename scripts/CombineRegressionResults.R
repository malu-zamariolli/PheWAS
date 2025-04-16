# Script to combine multinomial results

####################
### Packages
#####################
if(!require(data.table)) install.packages("data.table")
library(data.table)

########################
### Inputs and outputs
###########################
# input
input_files <- snakemake@input

# output
res_all <- snakemake@output[["res_all"]]
res_nom <- snakemake@output[["res_nom"]]

########################
### Load inputs 
###########################

# Initialize an empty list to store data frames
data_frames <- list()

# Iterate over the input files and read them into data frames
for (file in input_files) {
  # Read the file into a data frame
  df <- as.data.frame(fread((file)))
  # Append the data frame to the list
  data_frames <- append(data_frames, list(df))
}


########################
### Combine dataframes
###########################
df.all <- do.call(rbind, data_frames) 

######################
### Filter by p-value
###########################
df.nom <- df.all[df.all$p.value < 0.05, ]

######################
### Save
###########################
# All
fwrite(df.all, res_all, quote = FALSE, col.names = TRUE, 
       row.names = FALSE, sep = "\t", na = NA)

# Nominal
fwrite(df.nom, res_nom, quote = FALSE, col.names = TRUE, 
       row.names = FALSE, sep = "\t", na = NA)

