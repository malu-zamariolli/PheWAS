######### Non-parametric testing for  SNP in additive model ##################

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
snps <- snakemake@input[["snps"]]
traits <- snakemake@input[["traits"]]

# Outputs 
all_res <- snakemake@output[["all_res"]]
descriptive <- snakemake@output[["descriptive"]]


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

########################################
### Load and prepare input - Genotype and traits
#########################################
# SNPs
snps <- as.data.frame(fread(snps))
# Traits
traits <- as.data.frame(fread(traits))

# Combine traits and SNPs
df.raw <- left_join(traits, snps)

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

# Filter nominal results just in case (input is already nominal)
df.nominal <- df.pvalue[df.pvalue$p.value < 0.05, ]

####################################################
##### Non-parametric testing - Linear Regression
##################################################
# For linear regression, we'll use Wilcoxon Rank Test
print("Starting linear regression - Wilcoxon Rank Test")

df.linear <- df.nominal[df.nominal$method %in% "LinearRegression", ]

wt.list <- list()
for (i in 1:nrow(df.linear)) {
  pheno <- df.linear[i, "Pheno"]
  snp <- df.linear[i, "term"]
  formula <- as.formula(paste(pheno, "~", "as.factor(", snp, ")"))
  wt.test <- tryCatch({
    broom::tidy(wilcox.test(formula, na.action = na.omit, data = df.raw))
    }, error = function(e) {
        # Return dataframe with error
        data.frame(
            statistic = NA,
            p.value = NA,
            method = "error",
            alternative = "error"
        )
    })
  wt.test$Pheno <- pheno
  wt.test$SNP <- snp
  wt.test$Group <- "full"
  wt.list[[i]] <- wt.test
}

df.linear.wt <- do.call(rbind, wt.list)

linear.confirmed <- df.linear.wt[df.linear.wt$p.value < 0.05, ]

rm(wt.list, pheno, snp, wt.test, formula, i, df.linear)

####################################################
##### Non-parametric testing - Logistic Regression
##################################################
# Cochran Armitage to consider ordinal nature of additive model
print("Starting logistic regression - Fisher's test")

df.logistic <- df.nominal[df.nominal$method %in% "LogisticRegression", ]

fi.list <- list()
for (i in 1:nrow(df.logistic)) {
  pheno <- df.logistic[i, "Pheno"]
  snp <- df.logistic[i, "term"]
  table.contigency <- table(df.raw[[pheno]], df.raw[[snp]])
  fi.test <- tryCatch({
    broom::tidy(fisher.test(table.contigency))
    }, error = function(e) {
        # Return dataframe with error
        data.frame(
            estimate = NA,
            p.value = NA,
            conf.low = NA,
            conf.high = NA,
            method = "error", 
            alternative = "error"
        )
    })
  fi.test$Pheno <- pheno
  fi.test$SNP <- snp
  fi.test$Group <- "full"
  fi.list[[i]] <- fi.test
}

df.logistic.fi <- do.call(rbind, fi.list)

logistic.confirmed <- df.logistic.fi[df.logistic.fi$p.value < 0.05, ]

rm(fi.list, pheno, snp, fi.test, i, df.logistic, table.contigency)

####################################################
##### Non-parametric testing - Multinomial Regression
##################################################
# Cochran Armitage to consider ordinal nature of additive model
print("Starting Multinomial regression - Fisher's Test")

df.mult <- df.nominal[df.nominal$method %in% "MultiLogReg", ]

# Select results from full, to obtain y.level
ord.mult <- ord[ord$term %in% df.mult$term & ord$Pheno %in% df.mult$Pheno &
                       ord$estimate %in% df.mult$estimate, ]

uno.mult <- uno[uno$term %in% df.mult$term & uno$Pheno %in% df.mult$Pheno &
                  uno$estimate %in% df.mult$estimate, ]

df.mult <- rbind(ord.mult, uno.mult)

fi.list <- list()
for (i in 1:nrow(df.mult)) {
  pheno <- df.mult[i, "Pheno"]
  snp <- df.mult[i, "term"]
  group <- df.mult[i, "y.level"]
  # Filter groups of interest
  df.test <- df.raw[, c(pheno, snp)]
  ref_group <- min(unique(df.test[[pheno]]), na.rm = TRUE)
  df.test <- df.test[df.test[[pheno]] %in% c(1, group), ]
  # Table and test
  table.contigency <- table(df.test[[pheno]], df.test[[snp]])
  fi.test <- tryCatch({
    broom::tidy(fisher.test(table.contigency))
  }, error = function(e) {
    data.frame(
        estimate = NA,
        p.value = NA,
        conf.low = NA,
        conf.high = NA,
        method = "error", 
        alternative = "error"
    )
  })
  fi.test$Pheno <- pheno
  fi.test$SNP <- snp
  fi.test$Group <- group
  fi.list[[i]] <- fi.test
}

df.mult.fi <- do.call(rbind, fi.list)

multinomial.confirmed <- df.mult.fi[df.mult.fi$p.value < 0.05, ]

rm(fi.list, pheno, snp, fi.test, i, df.mult, table.contigency, df.test, group)

####################################################
##### Non-parametric testing - Ordinal Regression
##################################################
# Jonckheere-Terpstra will be used to test trend with permutatios due to N > 100
print("Starting Ordinal regression - Wilcoxon Rank Test")

df.ord <- df.nominal[df.nominal$method %in% "OrdinalLogisticRegression", ]

wt.list <- list()
for (i in 1:nrow(df.ord)) {
  pheno <- df.ord[i, "Pheno"]
  snp <- df.ord[i, "term"]
  formula <- as.formula(paste(pheno, "~", "as.factor(", snp, ")"))
  wt.test <- tryCatch({
    broom::tidy(wilcox.test(formula, na.action = na.omit, data = df.raw))
    }, error = function(e) {
        # Return dataframe with error
        data.frame(
            statistic = NA,
            p.value = NA,
            method = "error",
            alternative = "error"
        )
    })
  wt.test$Pheno <- pheno
  wt.test$SNP <- snp
  wt.test$Group <- "full"
  wt.list[[i]] <- wt.test
}

df.ord.wt <- do.call(rbind, wt.list)

ordinal.confirmed <- df.ord.wt[df.ord.wt$p.value < 0.05, ]

rm(wt.list, pheno, i, snp, wt.test, df.ord)

####################################################
##### Combine results and save
##################################################
print("Combining results")
# Combine rows 
df.confirmed <- bind_rows(df.ord.wt, df.mult.fi, df.logistic.fi, df.linear.wt)

fwrite(df.confirmed, all_res, quote = FALSE, 
    row.names = FALSE, col.names = TRUE, sep = "\t", na = NA)

####################################################
##### Create and save descriptive
##################################################
n_nom <- nrow(df.nominal)

n_mult <- nrow(multinomial.confirmed)
n_bin <- nrow(logistic.confirmed)
n_ord <- nrow(ordinal.confirmed)
n_linear <- nrow(linear.confirmed)

#################
text <- paste0("Total number of nominal associations: ", n_nom, "\n", 
               "Confirmed in multinomial regression: ", n_mult, "\n",
               "Confirmed in logistic regression: ", n_bin, "\n",
               "Confirmed in ordinal regression: ", n_ord, "\n",
               "Confirmed in linear regression: ", n_linear, "\n")

writeLines(text, descriptive)