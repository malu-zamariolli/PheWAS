######### Non-parametric testing for  SNP in additive model ##################

#### Packages ############
shh <- suppressPackageStartupMessages # https://stackoverflow.com/questions/18931006/how-to-suppress-warning-messages-when-loading-a-library

if(!require(data.table)) install.packages("data.table")
if(!require(dplyr)) install.packages("dplyr")
if(!require(DescTools)) install.packages("DescTools")
if(!require(broom)) install.packages("broom")
if(!require(clinfun)) install.packages("clinfun")
if(!require(stringr)) install.packages("stringr")


shh(library(data.table))
shh(library(dplyr))
shh(library(DescTools))
shh(library(broom))
shh(library(clinfun))
shh(library(stringr))

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
# For linear regression, we'll use Kruskall Wallis
print("Starting linear regression - Kruskall Wallis")

df.linear <- df.nominal[df.nominal$method %in% "LinearRegression", ]

kw.list <- list()
for (i in 1:nrow(df.linear)) {
  pheno <- df.linear[i, "Pheno"]
  snp <- df.linear[i, "term"]
  formula <- as.formula(paste(pheno, "~", "as.factor(", snp, ")"))
  kw.test <- tryCatch({
    broom::tidy(kruskal.test(formula, data = df.raw))
    }, error = function(e) {
        # Return dataframe with error
        data.frame(
            statistic = NA,
            p.value = NA,
            parameter = NA,
            method = "error"
        )
    })
  kw.test$Pheno <- pheno
  kw.test$SNP <- snp
  kw.test$Group <- "full"
  kw.list[[i]] <- kw.test
}

df.linear.kw <- do.call(rbind, kw.list)

linear.confirmed <- df.linear.kw[df.linear.kw$p.value < 0.05, ]

rm(kw.list, pheno, snp, kw.test, formula, i, df.linear)

####################################################
##### Non-parametric testing - Logistic Regression
##################################################
# Cochran Armitage to consider ordinal nature of additive model
print("Starting logistic regression - Cochran Armitage")

df.logistic <- df.nominal[df.nominal$method %in% "LogisticRegression", ]

# List for Cochran
ca.list <- list()
# List for Fisher
fi.log.list <- list()

for (i in 1:nrow(df.logistic)) {
  pheno <- df.logistic[i, "Pheno"]
  snp <- df.logistic[i, "term"]
  table.contigency <- table(df.raw[[pheno]], df.raw[[snp]])
  ca.test <- tryCatch({
    broom::tidy(DescTools::CochranArmitageTest(table.contigency))
    }, error = function(e) {
        # Return dataframe with error
        data.frame(
            statistic = NA,
            p.value = NA,
            parameter = NA,
            method = "error", 
            alternative = "error"
        )
    })
  ca.test$Pheno <- pheno
  ca.test$SNP <- snp
  ca.test$Group <- "full"
  ca.list[[i]] <- ca.test
  # Fisher 
    fi.log <- tryCatch({
    broom::tidy(fisher.test(table.contigency))
    }, error = function(e) {
        # Return dataframe with error
        data.frame(
            p.value = NA,
            method = "error", 
            alternative = "error"
        )
    })
  fi.log <- fi.log[, c("p.value", "method", "alternative")]
  fi.log$Pheno <- pheno
  fi.log$SNP <- snp
  fi.log$Group <- "full"
  fi.log.list[[i]] <- fi.log
}


df.logistic.ca <- do.call(rbind, ca.list)

df.logistic.fi <- do.call(rbind, fi.log.list)

logistic.confirmed <- df.logistic.ca[df.logistic.ca$p.value < 0.05, ]

rm(ca.list, pheno, snp, ca.test,i, df.logistic, table.contigency, fi.log.list)

####################################################
##### Non-parametric testing - Multinomial Regression
##################################################
# Cochran Armitage to consider ordinal nature of additive model
print("Starting Multinomial regression - Cochran Armitage")

df.mult <- df.nominal[df.nominal$method %in% "MultiLogReg", ]

# Select results from full, to obtain y.level
ord.mult <- ord[ord$term %in% df.mult$term & ord$Pheno %in% df.mult$Pheno &
                       ord$estimate %in% df.mult$estimate, ]

uno.mult <- uno[uno$term %in% df.mult$term & uno$Pheno %in% df.mult$Pheno &
                  uno$estimate %in% df.mult$estimate, ]

df.mult <- rbind(ord.mult, uno.mult)

# Cochran list
ca.list <- list()
# Fisher list
fi.mult.list <- list()

for (i in 1:nrow(df.mult)) {
  pheno <- df.mult[i, "Pheno"]
  snp <- df.mult[i, "term"]
  group <- df.mult[i, "y.level"]
  # Filter groups of interest
  df.test <- df.raw[, c(pheno, snp)]
  ref_group <- min(unique(df.test[[pheno]]), na.rm = TRUE)
  df.test <- df.test[df.test[[pheno]] %in% c(ref_group, group), ]
  # Table and test
  table.contigency <- table(df.test[[pheno]], df.test[[snp]])
  ca.test <- tryCatch({
    broom::tidy(DescTools::CochranArmitageTest(table.contigency))
  }, error = function(e) {
    data.frame(
        statistic = NA,
        p.value = NA, 
        parameter = NA,
        method = "error",
        alternative = "error"
    )
  })
  ca.test$Pheno <- pheno
  ca.test$SNP <- snp
  ca.test$Group <- group
  ca.list[[i]] <- ca.test
  # Fisher
  fi.mult <- tryCatch({
    broom::tidy(fisher.test(table.contigency))
  }, error = function(e) {
    data.frame(
        p.value = NA,
        method = "error", 
        alternative = "error"
    )
  })
  fi.mult <- fi.mult[, c("p.value", "method", "alternative")]
  fi.mult$Pheno <- pheno
  fi.mult$SNP <- snp
  fi.mult$Group <- group
  fi.mult.list[[i]] <- fi.mult
}

df.mult.ca <- do.call(rbind, ca.list)
df.mult.fi <- do.call(rbind, fi.mult.list)

multinomial.confirmed <- df.mult.ca[df.mult.ca$p.value < 0.05, ]

rm(ca.list, pheno, snp, ca.test,i, df.mult, table.contigency, df.test, group, fi.mult.list)

####################################################
##### Non-parametric testing - Ordinal Regression
##################################################
# Jonckheere-Terpstra will be used to test trend with permutatios due to N > 100
print("Starting Ordinal regression - Jonckheere-Terpstra")

df.ord <- df.nominal[df.nominal$method %in% "OrdinalLogisticRegression", ]

# Jonckheere-Terpstra list
jt.list <- list()
# Fisher list 
fi.ord.list <- list()

for (i in 1:nrow(df.ord)) {
  pheno <- df.ord[i, "Pheno"]
  snp <- df.ord[i, "term"]
  jt.test <- tryCatch({
    broom::tidy(jonckheere.test(df.raw[[pheno]], df.raw[[snp]], nperm = 1000))
  }, error = function(e) {
    data.frame(
        statistic = NA,
        p.value = NA,
        method = "error", 
        alternative = "error"
    )
  })
  jt.test$Pheno <- pheno
  jt.test$SNP <- snp
  jt.test$Group <- "full"
  jt.list[[i]] <- jt.test
  # Fisher
  table.contigency <- table(df.raw[[pheno]], df.raw[[snp]])
  fi.ord <- tryCatch({
    broom::tidy(fisher.test(table.contigency, simulate.p.value = TRUE, B = 10000))
  }, error = function(e) {
    data.frame(
        p.value = NA,
        method = "error", 
        alternative = "error"
    )
  })
  fi.ord <- fi.ord[, c("p.value", "method", "alternative")]
  fi.ord$Pheno <- pheno
  fi.ord$SNP <- snp
  fi.ord$Group <- "full"
  fi.ord.list[[i]] <- fi.ord
}

df.ord.jt <- do.call(rbind, jt.list)
df.ord.fi <- do.call(rbind, fi.ord.list)

ordinal.confirmed <- df.ord.jt[df.ord.jt$p.value < 0.05, ]

rm(jt.list, pheno, i, snp, jt.test, df.ord)

####################################################
##### Combine results and save
##################################################
print("Combining results")
# Combine rows 
df.confirmed <- bind_rows(df.ord.jt, df.ord.fi, df.mult.ca, df.mult.fi, df.logistic.ca, df.logistic.fi, df.linear.kw)

# Remove spaces, \n and \t in method (that create weird outputs)
df.confirmed$method <- stringr::str_remove_all(df.confirmed$method, " ")
df.confirmed$method <- stringr::str_remove_all(df.confirmed$method, "\n")
df.confirmed$method <- stringr::str_remove_all(df.confirmed$method, "\t")


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
               "Initial assesment of significance - fisher needs to be checked",  "\n",
               "Confirmed in multinomial regression: ", n_mult, "\n",
               "Confirmed in logistic regression: ", n_bin, "\n",
               "Confirmed in ordinal regression: ", n_ord, "\n",
               "Confirmed in linear regression: ", n_linear, "\n")

writeLines(text, descriptive)