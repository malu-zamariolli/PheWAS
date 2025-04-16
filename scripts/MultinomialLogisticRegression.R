
# Multinomial logistic regression

####### Packages #####################

if(!require(nnet)) install.packages("nnet")
if(!require(data.table)) install.packages("data.table")
if(!require(dplyr)) install.packages("dplyr")
if(!require(broom)) install.packages("broom")
if(!require(car)) install.packages("car")

library(nnet)
library(data.table)
library(dplyr)
library(broom)
library(car)

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

# Create function for multinominal
multinomial <- function(df, pheno, cov, geno) {
  # Organize covariates for formula
  cov_formula <- paste(cov, collapse = " + ")
  # Construct formulas
  formula_basic <- reformulate(cov, response = pheno)
  formula_full <- reformulate(c(cov, geno), response = pheno)
  # Model with covariates only
  fit_basic <- multinom(formula_basic, data = df, model = TRUE)  # model TRUE is required for ANOVA below
  # Model with genotype
  fit_full <- multinom(formula_full, data = df, model = TRUE)
  table.full <- broom::tidy(fit_full, conf.int = TRUE)
  # Add OR to model (exponenciate logOR)
  table.full$OR <- exp(table.full$estimate)
  table.full$OR.conf.low <- exp(table.full$conf.low)
  table.full$OR.conf.high <- exp(table.full$conf.high)
  # Test if model with genotype is significant in comparison with basic model
  model.aov <-  anova(fit_basic, fit_full) # Test if models differ
  # Add it to table
  table.full$method <- "MultiLogReg"
  table.full$anova.p <- model.aov$`Pr(Chi)`[2]
  # Create column with both models' convergence status
  table.full$convergence <- ifelse(fit_basic$convergence == 0 & fit_full$convergence == 0, "YES", "NO")
  return(table.full)
  # Obs: nnet:: should not be specify with multinom() or PseudoR2 won't work
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
# Convert all traits to factors
trait <- trait %>%
  mutate(across(-ID, as.factor))

# Covariates: convert characters to factor
covariates <- covariates %>%
  mutate_if(is.character, as.factor)

# Merge pheno, cov and geno dataframe
df <- inner_join(trait, snp)
df <- inner_join(df, covariates)

############################################
### Loop throught phenotypes
###########################################

multinomial.res <- list()
vif.results <- list()

for (i in 2:ncol(trait)) {
  pheno <- names(trait)[i]
  geno <- names(snp)[2] # always the same
  cov <- codebook[codebook$NewVariableName %in% pheno, "covariates"]
  cov <- unlist(strsplit(cov, ","))
  # Prepare df for model (only columns of interest) and exclude NAs
  df.model <- df %>% select(all_of(c("ID", pheno, geno, cov)))
  df.model <- na.omit(df.model)
  N <- nrow(df.model)
  # Run multinomial with tryCatch
  table.full <- tryCatch({
    multinomial(df.model, pheno, cov, geno)
  }, error = function(e) {
    # Return a dataframe with "error" in each column if an error occurs
    data.frame(
      y.level = "error",
      term = "error",
      estimate = NA,
      std.error = NA,
      statistic = NA,
      p.value = NA,
      conf.low = NA,
      conf.high = NA,
      OR = NA,
      OR.conf.low = NA,
      OR.conf.high = NA,
      method = "error",
      anova.p = NA,
      convergence = "error",
      stringsAsFactors = FALSE
    )
  })
  table.full$Pheno <- pheno # add phenotype to y.level column
  table.full$N <- N
  multinomial.res[[i-1]] <- table.full
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
multinomial.combined <- do.call(rbind, multinomial.res)

vif.results.combined <- do.call(rbind, vif.results)

######################################
### Save
#####################################

fwrite(multinomial.combined, res_all, quote = FALSE, col.names = TRUE, 
       row.names = FALSE, sep = "\t", na = NA)

fwrite(vif.results.combined, vif.output, quote = FALSE, col.names = TRUE, 
       row.names = FALSE, sep = "\t", na = NA)
