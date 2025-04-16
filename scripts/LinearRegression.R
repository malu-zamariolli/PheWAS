
# Linear regression

####### Packages #####################

if(!require(data.table)) install.packages("data.table")
if(!require(dplyr)) install.packages("dplyr")
if(!require(broom)) install.packages("broom")
if(!require(car)) install.packages("car")
if(!require(ggplot2)) install.packages("ggplot2")
if(!require(ggfortify)) install.packages("ggfortify")
if(!require(gridExtra)) install.packages("gridExtra")

library(data.table)
library(dplyr)
library(broom)
library(car)
library(ggplot2)
library(ggfortify)
library(gridExtra)

#################################################
#### Inputs and output ###################
#################################################

# Inputs
snp <- snakemake@input[["snp"]] 
trait <- snakemake@input[["trait"]]
cov <- snakemake@input[["cov"]]
codebook <- snakemake@input[["codebook"]]

# Outputs
res_all <- snakemake@output[["res_all"]]
vif.output <- snakemake@output[["vif"]]
plot <- snakemake@output[["plot"]]

#################################################
#### Load inputs ###################
#################################################

snp <- as.data.frame(fread(snp))
trait <- as.data.frame(fread(trait))
covariates <- as.data.frame(fread(cov))
codebook <- as.data.frame(fread(codebook))

##############################################
#### Functions
#############################################
# Create function for linear regression and plot
linear <- function(df, pheno, cov, geno) {
  # Organize covariates for formula
  cov_formula <- paste(cov, collapse = " + ")
  # Construct formulas
  formula_basic <- reformulate(cov, response = pheno)
  formula_full <- reformulate(c(cov, geno), response = pheno)
  # Model with covariates only
  fit_basic <-lm(formula_basic, data = df, model = TRUE)  # model TRUE is required for ANOVA below
  # Model with genotype
  fit_full <- lm(formula_full, data = df, model = TRUE)
  table.full <- broom::tidy(fit_full, conf.int = TRUE)
  # Add R2 of models
  table.full$AdjustedR2_null <- summary(fit_basic)$adj.r.squared
  table.full$AdjustedR2_full  <- summary(fit_full)$adj.r.squared
  # Test if model with genotype is significant in comparison with basic model
  model.aov <-  anova(fit_basic, fit_full) # Test if models differ
  # Add it to table
  table.full$anova.p <- model.aov$`Pr(>F)`[2]
  table.full$method <- "LinearRegression"
  # Plot normality and homecedasticity of residuals
  p <- autoplot(fit_full, which = c(1, 2)) + theme_minimal() + ggtitle(paste0(pheno, "-", geno))
  # return
  return(list(table.full, p))
}

# Create function for vif
vif_function <- function(df, pheno, cov, geno) {
  formula_full <- reformulate(c(cov, geno), response = pheno)
  df[, pheno] <- as.numeric(df[, pheno])
  m <- lm(formula_full, data = df, model = TRUE)
  vif_values <- as.data.frame(car::vif(m))
  vif_values <- data.frame(predictor = rownames(vif_values), vif = vif_values)
  vif_values$Pheno <- pheno 
  return(vif_values)
  
}

############################################
### Prepare dataframe
###########################################
# Convert all traits to continuous
trait <- trait %>%
  mutate(across(-ID, as.numeric))

# Covariates: convert characters to factor
covariates <- covariates %>%
  mutate_if(is.character, as.factor)

# Merge pheno, cov and geno dataframe
df <- inner_join(trait, snp)
df <- inner_join(df, covariates)

##############################################
#### Loop through phenotypes
#############################################

linear.res <- list()
vif.results <- list()
plot.list <- list()

for (i in 2:ncol(trait)) {
  pheno <- names(trait)[i]
  geno <- names(snp)[2] # always the same
  cov <- codebook[codebook$NewVariableName %in% pheno, "covariates"]
  cov <- unlist(strsplit(cov, ","))
  # Prepare df for model (only columns of interest) and exclude NAs
  df.model <- df %>% select(all_of(c("ID", pheno, geno, cov)))
  df.model <- na.omit(df.model)
  N <- nrow(df.model)
  # Run linear model with tryCatch
  linear.list <- tryCatch({
    linear(df.model, pheno, cov, geno)
  }, error = function(e) {
    # Return a list with a dataframe containing "error" in each column and an empty ggmultiplot object if an error occurs
    list(
      data.frame(
        term = "error",
        estimate = NA,
        std.error = NA,
        statistic = NA,
        p.value = NA,
        conf.low = NA,
        conf.high = NA,
        AdjustedR2_null = NA,
        AdjustedR2_full = NA,
        anova.p = NA,
        method = "error",
        stringsAsFactors = FALSE
      ),
      ggplot() + theme_void()
    )
  })
  # Extract results table
  table.full <- linear.list[[1]]
  table.full$Pheno <- pheno # add phenotype to y.level column
  table.full$N <- N
  # Store it in list 
  linear.res[[i-1]] <- table.full
  # Store plot in list
  p <- linear.list[[2]]
  plot.list[[i-1]] <- p
  # Run VIF with tryCatch
  vif <- tryCatch({
    vif_function(df.model, pheno, cov, geno)
  }, error = function(e) {
    # Return a dataframe with "error" in each column if an error occurs
    data.frame(
      predictor = "error",
      car..vif.m. = NA,
      Pheno = pheno,
      stringsAsFactors = FALSE
    )
  })
  
  vif.results[[i-1]] <- vif
}


######################################
### Unlist df and combine
#####################################
linear.combined <- do.call(rbind, linear.res)

vif.results.combined <- do.call(rbind, vif.results)

######################################
### Save
#####################################
# Linear results
fwrite(linear.combined, res_all, quote = FALSE, col.names = TRUE, 
       row.names = FALSE, sep = "\t", na = NA)
# VIF
fwrite(vif.results.combined, vif.output, quote = FALSE, col.names = TRUE, 
       row.names = FALSE, sep = "\t", na = NA)

# Plots
# Open PDF device
pdf(plot, width = 10, height = 10)

# Loop through each plot and save it to PDF
for (i in seq_along(plot.list)) {
  plot_obj <- plot.list[[i]]
  print(plot_obj)  # Print the plot to the current PDF device
}

# Close PDF device
dev.off()



