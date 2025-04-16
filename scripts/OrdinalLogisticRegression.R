####### Ordinal logistic regression X Multinomial Regression (Proportional Odds) ###################

########## Packages ########################################
if(!require(data.table)) install.packages("data.table")
if(!require(dplyr)) install.packages("dplyr")
if(!require(MASS)) install.packages("MASS")
if(!require(car)) install.packages("car")
if(!require(broom)) install.packages("broom")
if(!require(nnet)) install.packages("nnet")

library(data.table)
library(dplyr)
library(MASS)
library(car)
library(broom)
library(nnet)


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
prop_odds <- snakemake@output[["prop_odds"]]


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
# Ordinal Logistic Regression
proportional_odds <- function(df, pheno, cov, geno) {
  # Organize covariates for formula
  cov_formula <- paste(cov, collapse = " + ")
  # Construct formulas
  formula_full <- reformulate(c(cov, geno), response = pheno)
  # Model with genotype
  fit_full <- MASS::polr(formula_full, data = df, Hess = T)
  # Assess proportional odds assumption - Brant Test
  prop.odds <- car::poTest(fit_full) 
  # Output - overall
  chi_square_value <- prop.odds$chisq
  df_overall <- prop.odds$df
  p_value <- 1 - pchisq(chi_square_value, df_overall)
  overall <- data.frame(Term = c("Overall"), 
                        p.value = p_value, 
                        Chisquare = chi_square_value, 
                        df = df_overall)
  # Predictors p-value
  chi_pred <- prop.odds$chisq.p
  df_pred <- prop.odds$df.p
  p_predictors <- 1 - pchisq(chi_pred, df_pred)
  predictors <- data.frame(Term = prop.odds$coef.names, 
                           p.value = p_predictors, 
                           Chisquare = chi_pred,
                           df = df_pred)
  # Combine data
  combine.df <- rbind(overall, predictors)
  # Filter overall and SNP
  prop.odds.df <- combine.df[combine.df$Term %in% c("Overall", geno), ]
  # Add trait being tested and snp 
  prop.odds.df$Pheno <- pheno
  prop.odds.df$SNP <- geno
  return(prop.odds.df)
}
# p-value can't be extracted from package: calculation https://www.r-bloggers.com/2022/05/calculate-the-p-value-from-chi-square-statistic-in-r/ 

# Ordinal Logistic Regression
ordinal_regression <- function(df, pheno, cov, geno) {
  # Organize covariates for formula
  cov_formula <- paste(cov, collapse = " + ")
  # Construct formulas
  formula_basic <- reformulate(cov, response = pheno)
  formula_full <- reformulate(c(cov, geno), response = pheno)
  # Model with covariates only
  fit_basic <- MASS::polr(formula_basic, data = df, Hess = T)
  # Model with genotype
  fit_full <- MASS::polr(formula_full, data = df, Hess = T)
  # Build table.full
  sum.model <- summary(fit_full)
  coef <- as.data.frame(sum.model$coefficients)
  table.full <- data.frame(y.level = "ordinal",cbind(term = row.names(coef), coef)) # y.level is too match multinomial results
  # p-value
  table.full$p.value <- pnorm(abs(table.full$t.value), lower.tail = FALSE) * 2
  # Add confidence intervals
  CI.table <- as.data.frame(confint(fit_full, level = 0.95))
  CI.table <- data.frame(cbind(term = row.names(CI.table), CI.table))
  colnames(CI.table)[c(2,3)] <- c("conf.low", "conf.high")
  table.full <- left_join(table.full, CI.table)
  # Add OR
  table.full$OR <- exp(table.full$Value)
  table.full$OR.conf.low <- exp(table.full$conf.low)
  table.full$OR.conf.high <- exp(table.full$conf.high)
  # Test if model with genotype is significant in comparison with basic model
  model.aov <-  anova(fit_basic, fit_full) # Test if models differ
  # Add it to table
  table.full$method <- "OrdinalLogisticRegression"
  table.full$anova.p <- model.aov$`Pr(Chi)`[2]
  # Convergence status
  table.full$convergence <- fit_full$convergence
  # Change colnames to match multinomial
  colnames(table.full) <- c("y.level","term", "estimate",	"std.error", 	"statistic",
                            "p.value", "conf.low",	"conf.high",	"OR",	
                            "OR.conf.low",	"OR.conf.high",	"method",	
                            "anova.p",	"convergence")
 # Return
  return(table.full)
}

# p-value reference https://stats.oarc.ucla.edu/r/dae/ordinal-logistic-regression/ 
#https://www.youtube.com/watch?v=wZLlL173uVc&t=1s



# Multinomial regression
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

# VIF
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
#######################################################
# Convert all traits to factors
trait <- trait %>%
  mutate(across(-ID, ~ factor(.x, ordered = TRUE)))

# Covariates: convert characters to factor
covariates <- covariates %>%
  mutate_if(is.character, as.factor)

# Merge pheno, cov and geno dataframe
df <- inner_join(trait, snp)
df <- inner_join(df, covariates)


###################################################
### Loop throught phenotypes - prop odds and VIF
################################################
# Start empty lists
prop.odds.res <- list()
vif.results <- list()


for (i in 2:ncol(trait)) {
  pheno <- names(trait)[i]
  geno <- names(snp)[2] # always the same
  cov <- codebook[codebook$NewVariableName %in% pheno, "covariates"]
  cov <- unlist(strsplit(cov, ","))
  # Prepare df for model (only columns of interest) and exclude NAs
  df.model <- df %>% dplyr::select(all_of(c("ID", pheno, geno, cov)))
  df.model <- na.omit(df.model)
  N <- nrow(df.model)
  # Run multinomial with tryCatch
  prop.odds.df <- tryCatch({
    proportional_odds(df.model, pheno, cov, geno)
  }, error = function(e) {
    # Return a list with df and printed message
    data.frame(Term = "error", 
               p.value = "NA",
               Chisquare = NA,
               df = NA,
               Pheno = pheno,
               SNP = geno,
               stringsAsFactors = FALSE)
  })
  prop.odds.res[[i-1]] <- prop.odds.df
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


# Rbind list dfs
prop.odds.combined <- do.call(rbind, prop.odds.res)

vif.results.combined <- do.call(rbind, vif.results)


###################################################
### Split dfs for Multinomial and Ordinal
################################################
# If prop odds <= 0.05 (overall or SNP), do multinomial (assumption was broken)
mult_leads <- prop.odds.combined[prop.odds.combined$p.value < 0.05, ]
mult_leads <- mult_leads %>% dplyr::select(Pheno, SNP)
mult_leads <- unique(mult_leads)
  

# If prop odds > 0.05, do multinomial
ord_leads <- prop.odds.combined %>% dplyr::select(Pheno, SNP)
ord_leads <- setdiff(ord_leads, mult_leads)

###################################################
### Loop throught phenotypes - ordinal regression
################################################
ordinal.res <- list()

if (nrow(ord_leads) > 0) {
  for (i in 1:nrow(ord_leads)) {
    pheno <- ord_leads[i, "Pheno"]
    geno <- names(snp)[2] # always the same
    cov <- codebook[codebook$NewVariableName %in% pheno, "covariates"]
    cov <- unlist(strsplit(cov, ","))
    # Prepare df for model (only columns of interest) and exclude NAs
    df.model <- df %>% dplyr::select(all_of(c("ID", pheno, geno, cov)))
    df.model <- na.omit(df.model)
    N <- nrow(df.model)
    # Run multinomial with tryCatch
    table.full <- tryCatch({
      ordinal_regression(df.model, pheno, cov, geno)
    }, error = function(e) {
      # Return a list with df and printed message
      data.frame(
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
        y.level = "error",
        stringsAsFactors = FALSE)
    })
    # Add pheno and N
    table.full$Pheno <- pheno 
    table.full$N <- N
    # Store df in list
    ordinal.res[[i]] <- table.full
  }
} else {
  print("No ordinal variables will be done with ordinal logistic - skipping loop")  
}

rm(geno, pheno, cov, df.model)


###################################################
### Loop throught phenotypes - multinomial regression
################################################
multinomial.res <- list()

if (nrow(mult_leads) > 0) {
  for (i in 1:nrow(mult_leads)) {
    pheno <- mult_leads[i, "Pheno"]
    geno <- names(snp)[2] # always the same
    cov <- codebook[codebook$NewVariableName %in% pheno, "covariates"]
    cov <- unlist(strsplit(cov, ","))
    # Prepare df for model (only columns of interest) and exclude NAs
    df.model <- df %>% dplyr::select(all_of(c("ID", pheno, geno, cov)))
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
        stringsAsFactors = FALSE)
    })
    table.full$Pheno <- pheno # add phenotype to y.level column
    table.full$N <- N
    multinomial.res[[i]] <- table.full
  }
} else {
    print("No ordinal variables will be done with multinomial - skipping loop")
}

###################################################
### Unlist dataframe and combine
################################################

ordinal.combined <- do.call(rbind, ordinal.res)
multinomial.combined <- do.call(rbind, multinomial.res)

if (is.null(ordinal.combined)) {
  all.results <- multinomial.combined
} else if (is.null(multinomial.combined)) {
  all.results <- ordinal.combined
} else {
  all.results <- rbind(ordinal.combined, multinomial.combined)
}

######################################
### Save
#####################################

fwrite(all.results, res_all, quote = FALSE, col.names = TRUE, 
       row.names = FALSE, sep = "\t", na = NA)

fwrite(vif.results.combined, vif.output, quote = FALSE, col.names = TRUE, 
       row.names = FALSE, sep = "\t", na = NA)

fwrite(prop.odds.combined, prop_odds, quote = FALSE, col.names = TRUE, 
       row.names = FALSE, sep = "\t", na = NA)
