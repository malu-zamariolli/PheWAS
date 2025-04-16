
# Logistic regression with firth regression

####### Packages #####################

if(!require(data.table)) install.packages("data.table")
if(!require(dplyr)) install.packages("dplyr")
if(!require(broom)) install.packages("broom")
if(!require(car)) install.packages("car")
if(!require(logistf)) install.packages("logistf")

library(data.table)
library(dplyr)
library(car)
library(broom)
library(logistf)


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

# Logistic regression
perform_logistic_regression <- function(df, pheno, cov, geno) {
  # Organize covariates for formula
  cov_formula <- paste(cov, collapse = " + ")
  # Construct formulas
  formula_basic <- reformulate(cov, response = pheno)
  formula_full <- reformulate(c(cov, geno), response = pheno)
  # Model with covariates only
  fit_basic <- glm(formula_basic, family = binomial(link = "logit"), data = df, model = TRUE)  # model TRUE is required for ANOVA below
  # Model with genotype
  fit_full <- glm(formula_full, family = binomial(link = "logit"), data = df, model = TRUE)
  table.full <- broom::tidy(fit_full, conf.int = TRUE)
  # Add OR to model (exponenciate logOR)
  table.full$OR <- exp(table.full$estimate)
  table.full$OR.conf.low <- exp(table.full$conf.low)
  table.full$OR.conf.high <- exp(table.full$conf.high)
  # Test if model with genotype is significant in comparison with basic model
  model.aov <- anova(fit_basic, fit_full, test="Chisq") # Test if models differ
  # Add it to table
  table.full$method <- "LogisticRegression"
  table.full$anova.p <- model.aov$`Pr(>Chi)`[2]
  # Convergence status
  table.full$Convergence <- fit_full$converged
  return(table.full)
}

# Firth Regression
perform_firth_regression <- function(df, pheno, cov, geno) {
  # Organize covariates for formula
  cov_formula <- paste(cov, collapse = " + ")
  # Construct formulas
  formula_basic <- reformulate(cov, response = pheno)
  formula_full <- reformulate(c(cov, geno), response = pheno)
  # Model with covariates only
  fit_basic <- logistf(formula_basic, data = df, model = TRUE)  # model TRUE is required for ANOVA below
  # Model with genotype
  fit_full <- tryCatch({
    logistf(formula_full, data = df, model = TRUE)
  }, warning = function(w) {
    warning.message <<- w$message  # Store warning message in global variable
    return(logistf(formula_full, data = df, model = TRUE)) 
  })
  #fit_full <- logistf(formula_full, data = df, model = TRUE)
  # Create table with results
  table.full <- data.frame(term = fit_full$terms,
                           estimate = fit_full$coefficients,
                           std.error = sqrt(diag(fit_full$var)),
                           statistic = qchisq(1-fit_full$prob,df=1),
                           p.value = fit_full$prob,
                           conf.low = fit_full$ci.lower,
                           conf.high = fit_full$ci.upper)
  # Add OR to model (exponenciate logOR)
  table.full$OR <- exp(table.full$estimate)
  table.full$OR.conf.low <- exp(table.full$conf.low)
  table.full$OR.conf.high <- exp(table.full$conf.high)
  # Test if model with genotype is significant in comparison with basic model
  model.aov <- anova(fit_basic, fit_full) # Test if models differ
  # Add it to table
  table.full$method <- "FirthLogisticRegression"
  table.full$anova.p <- model.aov$pval
  # Convergence status
  table.full$Convergence <- paste(warning.message, collapse = "; ")
  return(table.full)
}

# Calculate vif
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

logistic.res <- list()
vif.results <- list()

# Loop through each phenotype
for (i in 2:ncol(trait)) {
  pheno <- names(trait)[i]
  geno <- names(snp)[2] 
  cov <- codebook[codebook$NewVariableName %in% pheno, "covariates"]
  cov <- unlist(strsplit(cov, ","))
  # Prepare dataframe for model (select columns of interest and remove NAs)
  df.model <- df %>% select(all_of(c("ID", pheno, geno, cov)))
  df.model <- na.omit(df.model)
  N <- nrow(df.model)
  # Attempt logistic regression
  table.full <- tryCatch({
    perform_logistic_regression(df.model, pheno, cov, geno)
  }, error = function(e) {
    return(NULL)
  })
  if (!is.null(table.full)) {
    # Logistic regression succeeded
    print("Logistic done")
  } else {
    # Attempt Firth's regression if logistic regression fails
    table.full <- tryCatch({
      perform_firth_regression(df.model, pheno, cov, geno)
    }, error = function(e) {
      return(NULL)
    })
    
    if (!is.null(table.full)) {
      # Firth's regression succeeded
      print("Firth Regression done")
    } else {
      # Both regressions failed
      print("both regressions failed")
      table.full <- data.frame(
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
        Convergence = NA, 
        stringsAsFactors = FALSE
      )
    }
  }
  table.full$Pheno <- pheno
  table.full$N <- N
  logistic.res[[i-1]] <- table.full
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
logistic.combined <- do.call(rbind, logistic.res)

vif.results.combined <- do.call(rbind, vif.results)

######################################
### Save
#####################################

fwrite(logistic.combined, res_all, quote = FALSE, col.names = TRUE, 
       row.names = FALSE, sep = "\t", na = NA)

fwrite(vif.results.combined, vif.output, quote = FALSE, col.names = TRUE, 
       row.names = FALSE, sep = "\t", na = NA)
