
# A. Config file
configfile: "config_PheWAS.yaml"

# B. Define wildcard
wildcard_constraints:
     phewas = config["phewas"]
     
# C. non-wildcard definitions
codebook = config["CodeBook"]
traits = config["Traits"]
traits_raw = config["Traits_raw"]
plinkfiles = config["PlinkFiles"]
model = config["snp_model"]


#----------------------- Rule all --------------------------------------------------#
rule all:
    input: 
        expand("{phewas}/config/config_PheWAS.yaml", phewas = config["phewas"]),
        expand("{phewas}/results/multinomial_{phewas}_bonferroniSig.txt", phewas = config["phewas"]),
        expand("{phewas}/results/continuous_{phewas}_bonferroniSig.txt", phewas = config["phewas"]),
        expand("{phewas}/results/binary_{phewas}_bonferroniSig.txt", phewas = config["phewas"]),
        expand("{phewas}/results/ordinal_{phewas}_bonferroniSig.txt", phewas = config["phewas"]),
        expand("{phewas}/results/DescriptiveAfterBonferroni.txt", phewas = config["phewas"]),
        expand("{phewas}/SNPs/SNPdata_allSNPs_{phewas}.txt", phewas = config["phewas"]),
        expand("{phewas}/results/continuous_{phewas}_NormalityHeteroscedasticity.pdf", phewas = config["phewas"]),
        expand("{phewas}/results/{phewas}_vif.txt", phewas = config["phewas"]),
        expand("{phewas}/results/NeffGao/multinomial_{phewas}_NeffSig.txt", phewas = config["phewas"]),
        expand("{phewas}/results/NeffGao/continuous_{phewas}_NeffSig.txt", phewas = config["phewas"]),
        expand("{phewas}/results/NeffGao/binary_{phewas}_NeffSig.txt", phewas = config["phewas"]),
        expand("{phewas}/results/NeffGao/ordinal_{phewas}_NeffSig.txt", phewas = config["phewas"]),
        expand("{phewas}/results/NeffGao/DescriptiveAfterNeff.txt", phewas = config["phewas"]),
        expand("{phewas}/results/non_param/{phewas}_all_results.txt", phewas = config["phewas"]),
        expand("{phewas}/results/FDR/{phewas}_FDRadjustment_all.txt", phewas = config["phewas"])
  

#---------------------- Rules ---------------------------------------------------------#   
# Extract SNPs
rule extractSNPs:
    input:
        multiext(f'{plinkfiles}', ".bed", ".bim", ".fam"),
        snplist = "{phewas}/SNPlist/{phewas}.snps"
    output:
        raw = temp(multiext("{phewas}/SNPlist/{phewas}_recoded", ".raw", ".log"))
    conda: 
        "../envs/plink.yaml"
    params:
        binaries = f'{plinkfiles}',
        raw = "{phewas}/SNPlist/{phewas}_recoded",
        model = config["snp_model"]
    shell:
        """
        plink --bfile {params.binaries} --extract {input.snplist} --recode {params.model} --out {params.raw}
        """
# recode A gives encoding of additive model: 0, 1, 2

# Split SNPs in individual dfs
rule splitSNP:  
    input:
        raw = expand("{phewas}/SNPlist/{phewas}_recoded.raw", phewas=config["phewas"])
    output:
        splitSNP = temp(expand("{phewas}/SNPs/SNPdata_{snp}_{phewas}.txt", phewas=config["phewas"], snp=range(1, int(config["n_snps"]) + 1))),
        allSNPs = expand("{phewas}/SNPs/SNPdata_allSNPs_{phewas}.txt", phewas=config["phewas"])
    script:
        "scripts/SNPinputData.R"


rule covariates:
    input:
        ids = f'{plinkfiles}' + ".fam"
    output:
        cov = "{phewas}/covariates/covariates_{phewas}.txt"
    script:
        "scripts/covariates/{wildcards.phewas}_covariates.R"


rule traits:
    input:
        codebook = f'{codebook}',
        traits = f'{traits}'
    output:
        cont = temp(expand("{phewas}/traits/continuous_{chunk}_{phewas}.txt", phewas = config["phewas"], chunk = range(1, int(config["n_cont"]) + 1))),
        binary = temp(expand("{phewas}/traits/binary_{chunk}_{phewas}.txt", phewas = config["phewas"], chunk = range(1, int(config["n_bin"]) + 1))),
        unor = temp(expand("{phewas}/traits/unordered_{chunk}_{phewas}.txt", phewas = config["phewas"], chunk = range(1, int(config["n_uno"]) + 1))),
        ordi = temp(expand("{phewas}/traits/ordinal_{chunk}_{phewas}.txt", phewas = config["phewas"], chunk = range(1, int(config["n_ord"]) + 1)))
    params:
        chunks = config["chunks"]
    script:
        "scripts/TraitByType.R"


# Multinomial: categorical unordered
rule multinomial:
    input:
        snp = "{phewas}/SNPs/SNPdata_{snp}_{phewas}.txt",
        trait = "{phewas}/traits/unordered_{chunk}_{phewas}.txt",
        cov = "{phewas}/covariates/covariates_{phewas}.txt",
        codebook = f'{codebook}'
    output:
        res_all = temp("{phewas}/results/temp/multinomial_{snp}_{chunk}_{phewas}.txt"),
        vif = temp("{phewas}/results/temp/multinomial_{snp}_{chunk}_{phewas}_vif.txt")
    benchmark:
        "{phewas}/benchmarks/{snp}.{chunk}.multinomial.benchmark.txt"
    script:
        "scripts/MultinomialLogisticRegression.R"

rule combine_multinomial:
    input:
        expand("{phewas}/results/temp/multinomial_{snp}_{chunk}_{phewas}.txt", phewas = config["phewas"], chunk = range(1, int(config["n_uno"]) + 1), snp = range(1, int(config["n_snps"]) + 1)) 
    output:
        res_all = "{phewas}/results/multinomial_{phewas}_combined.txt",
        res_nom = "{phewas}/results/multinomial_{phewas}_nominal.txt"
    script:
        "scripts/CombineRegressionResults.R"


# Linear: continuous
rule linear:
    input:
        snp = "{phewas}/SNPs/SNPdata_{snp}_{phewas}.txt",
        trait = "{phewas}/traits/continuous_{chunk}_{phewas}.txt",
        cov = "{phewas}/covariates/covariates_{phewas}.txt",
        codebook = f'{codebook}'
    output:
        res_all = temp("{phewas}/results/temp/continuous_{snp}_{chunk}_{phewas}.txt"),
        vif = temp("{phewas}/results/temp/continuous_{snp}_{chunk}_{phewas}_vif.txt"),
        plot = temp("{phewas}/results/plots/linearRegression/continuous_{snp}_{chunk}_{phewas}.pdf")
    benchmark:
        "{phewas}/benchmarks/{snp}.{chunk}.linear.benchmark.txt"
    script:
        "scripts/LinearRegression.R"

rule combine_linear:
    input:
        expand("{phewas}/results/temp/continuous_{snp}_{chunk}_{phewas}.txt", phewas = config["phewas"], chunk = range(1, config["n_cont"] + 1), snp=range(1, int(config["n_snps"]) + 1)) 
    output:
        res_all = "{phewas}/results/continuous_{phewas}_combined.txt",
        res_nom = "{phewas}/results/continuous_{phewas}_nominal.txt"
    script:
        "scripts/CombineRegressionResults.R"

# Combine linear regression pdf files
rule combine_linear_plots:
    input:
        expand("{phewas}/results/plots/linearRegression/continuous_{snp}_{chunk}_{phewas}.pdf", phewas = config["phewas"], chunk = range(1, config["n_cont"] + 1), snp=range(1, int(config["n_snps"]) + 1)) 
    output:
        plots_all = "{phewas}/results/continuous_{phewas}_NormalityHeteroscedasticity.pdf"
    conda: 
        "../envs/qpdf.yaml"
    shell:
        """
        qpdf --empty --pages {input} -- {output.plots_all}
        """


# Logistic : binary
rule logistic:
    input:
        snp = "{phewas}/SNPs/SNPdata_{snp}_{phewas}.txt",
        trait = "{phewas}/traits/binary_{chunk}_{phewas}.txt",
        cov = "{phewas}/covariates/covariates_{phewas}.txt",
        codebook = f'{codebook}'
    output:
        res_all = temp("{phewas}/results/temp/binary_{snp}_{chunk}_{phewas}.txt"),
        vif = temp("{phewas}/results/binary_{snp}_{chunk}_{phewas}_vif.txt")
    benchmark:
        "{phewas}/benchmarks/{snp}.{chunk}.logistic.benchmark.txt"
    script:
        "scripts/LogisticRegression.R"

rule combine_logistic:
    input:
        expand("{phewas}/results/temp/binary_{snp}_{chunk}_{phewas}.txt", phewas = config["phewas"], chunk = range(1, config["n_bin"] + 1), snp=range(1, int(config["n_snps"]) + 1)) 
    output:
        res_all = "{phewas}/results/binary_{phewas}_combined.txt",
        res_nom = "{phewas}/results/binary_{phewas}_nominal.txt"
    script:
        "scripts/CombineRegressionResults.R"

# Ordinal Logistic Regression
rule ordinal:
    input:
        snp = "{phewas}/SNPs/SNPdata_{snp}_{phewas}.txt",
        trait = "{phewas}/traits/ordinal_{chunk}_{phewas}.txt",
        cov = "{phewas}/covariates/covariates_{phewas}.txt",
        codebook = f'{codebook}'
    output:
        res_all = temp("{phewas}/results/temp/ordinal_{snp}_{chunk}_{phewas}.txt"),
        vif = temp("{phewas}/results/temp/ordinal_{snp}_{chunk}_{phewas}_vif.txt"),
        prop_odds = "{phewas}/results/ordinal_{snp}_{chunk}_{phewas}_OverallPropOdds.txt"
    benchmark:
        "{phewas}/benchmarks/{snp}.{chunk}.ordinal.benchmark.txt"
    script:
        "scripts/OrdinalLogisticRegression.R"

rule combine_ordinal:
    input:
        expand("{phewas}/results/temp/ordinal_{snp}_{chunk}_{phewas}.txt", phewas = config["phewas"], chunk = range(1, config["n_ord"] + 1), snp=range(1, int(config["n_snps"]) + 1)) 
    output:
        res_all = "{phewas}/results/ordinal_{phewas}_combined.txt",
        res_nom = "{phewas}/results/ordinal_{phewas}_nominal.txt"
    script:
        "scripts/CombineRegressionResults.R"

# Bonferroni correction
rule bonferroni:
    input:
        codebook = f'{codebook}',
        mult = "{phewas}/results/multinomial_{phewas}_nominal.txt",
        cont = "{phewas}/results/continuous_{phewas}_nominal.txt",
        binary = "{phewas}/results/binary_{phewas}_nominal.txt",
        ordinal = "{phewas}/results/ordinal_{phewas}_nominal.txt"
    output:
        mult_sig = "{phewas}/results/multinomial_{phewas}_bonferroniSig.txt",
        cont_sig = "{phewas}/results/continuous_{phewas}_bonferroniSig.txt",
        binary_sig = "{phewas}/results/binary_{phewas}_bonferroniSig.txt",
        ordinal_sig = "{phewas}/results/ordinal_{phewas}_bonferroniSig.txt",
        summary = "{phewas}/results/DescriptiveAfterBonferroni.txt"
    params:
        n_snps = config["n_snps"]
    script:
        "scripts/BonferroniCorrection.R"

rule combine_vif:
    input:
        mult = expand("{phewas}/results/temp/multinomial_{snp}_{chunk}_{phewas}_vif.txt", phewas = config["phewas"], chunk = range(1, config["n_uno"] + 1), snp=range(1, int(config["n_snps"]) + 1)),
        cont = expand("{phewas}/results/temp/continuous_{snp}_{chunk}_{phewas}_vif.txt", phewas = config["phewas"], chunk = range(1, config["n_cont"] + 1), snp=range(1, int(config["n_snps"]) + 1)),
        binary = expand("{phewas}/results/binary_{snp}_{chunk}_{phewas}_vif.txt", phewas = config["phewas"], chunk = range(1, config["n_bin"] + 1), snp=range(1, int(config["n_snps"]) + 1)),
        ordinal = expand("{phewas}/results/temp/ordinal_{snp}_{chunk}_{phewas}_vif.txt", phewas = config["phewas"], chunk = range(1, config["n_ord"] + 1), snp=range(1, int(config["n_snps"]) + 1)) 
    output:
        "{phewas}/results/{phewas}_vif.txt"
    script:
        "scripts/CombineVIF.R"

#-------------------- Gao - Number of effective tests ----------------#
rule Neff_binary:
    input:
        df_input = f'{traits}',
        codebook = f'{codebook}'
    output:
        Neff = "{phewas}/results/NeffGao/Neff_binary_traits.txt"
    params:
        trait_type = "BINARY"
    benchmark:
        "{phewas}/benchmarks/Neff_binary_traits.benchmark.txt"
    script:
        "scripts/Neff_Gao.R"

rule Neff_continuous:
    input:
        df_input = f'{traits}',
        codebook = f'{codebook}'
    output:
        Neff = "{phewas}/results/NeffGao/Neff_continuous_traits.txt"
    params:
        trait_type = "CONTINUOUS"
    benchmark:
        "{phewas}/benchmarks/Neff_continuous_traits.benchmark.txt"
    script:
        "scripts/Neff_Gao.R"

# Categorical ordinal and non-ordinal will be considered together
rule Neff_categorical:
    input:
        df_input = f'{traits}',
        codebook = f'{codebook}'
    output:
        Neff = "{phewas}/results/NeffGao/Neff_categorical_traits.txt"
    params:
        trait_type = "ORDERED"
    benchmark:
        "{phewas}/benchmarks/Neff_categorical_traits.benchmark.txt"
    script:
        "scripts/Neff_Gao.R"

rule Neff_SNPs:
    input:
        df_input = "{phewas}/SNPs/SNPdata_allSNPs_{phewas}.txt",
        codebook = f'{codebook}'
    output:
        Neff = "{phewas}/results/NeffGao/Neff_SNPs.txt"
    params:
        trait_type = "SNP"
    benchmark:
        "{phewas}/benchmarks/Neff_SNPs.benchmark.txt"
    script:
        "scripts/Neff_Gao.R"

rule sig_Gao:
    input:
        mult = "{phewas}/results/multinomial_{phewas}_nominal.txt",
        cont = "{phewas}/results/continuous_{phewas}_nominal.txt",
        binary = "{phewas}/results/binary_{phewas}_nominal.txt",
        ordinal = "{phewas}/results/ordinal_{phewas}_nominal.txt",
        neff_bin = "{phewas}/results/NeffGao/Neff_binary_traits.txt",
        neff_cont = "{phewas}/results/NeffGao/Neff_continuous_traits.txt",
        neff_cat = "{phewas}/results/NeffGao/Neff_categorical_traits.txt",
        neff_snp = "{phewas}/results/NeffGao/Neff_SNPs.txt"
    output:
        mult_sig = "{phewas}/results/NeffGao/multinomial_{phewas}_NeffSig.txt",
        cont_sig = "{phewas}/results/NeffGao/continuous_{phewas}_NeffSig.txt",
        binary_sig = "{phewas}/results/NeffGao/binary_{phewas}_NeffSig.txt",
        ordinal_sig = "{phewas}/results/NeffGao/ordinal_{phewas}_NeffSig.txt",
        summary = "{phewas}/results/NeffGao/DescriptiveAfterNeff.txt"
    script:
        "scripts/Neff_GaoCorrection.R"

rule FDR:
    input:
        mult = "{phewas}/results/multinomial_{phewas}_combined.txt",
        cont = "{phewas}/results/continuous_{phewas}_combined.txt",
        binary = "{phewas}/results/binary_{phewas}_combined.txt",
        ordinal = "{phewas}/results/ordinal_{phewas}_combined.txt"
    output:
        fdr = "{phewas}/results/FDR/{phewas}_FDRadjustment_all.txt",
        fdr_sig = "{phewas}/results/FDR/{phewas}_FDRadjustment_sig.txt",
        summary = "{phewas}/results/FDR/DescriptiveAfterFDR.txt"
    script:
        "scripts/FDR_correction.R"

#----------------------------- Non-parametric testing ---------------------------------#
rule non_param:
    input:
        mult = "{phewas}/results/multinomial_{phewas}_nominal.txt",
        cont = "{phewas}/results/continuous_{phewas}_nominal.txt",
        binary = "{phewas}/results/binary_{phewas}_nominal.txt",
        ordinal = "{phewas}/results/ordinal_{phewas}_nominal.txt",
        snps = "{phewas}/SNPs/SNPdata_allSNPs_{phewas}.txt",
        traits = f'{traits_raw}'
    output:
        all_res = "{phewas}/results/non_param/{phewas}_all_results.txt",
        descriptive = "{phewas}/results/non_param/DescriptiveResults.txt"
    script:
        f"scripts/Non_parametric_{model}.R"

#----------------------------- Save config --------------------------------------------------#
#Copy config file to corresponding folder, so config informations are saved to each run
rule config:
    input:
        "config_PheWAS.yaml"
    output:
        "{phewas}/config/config_PheWAS.yaml"
    params:
        path = "{phewas}/config/"
    shell:
        """
        cp {input} {params.path}
        """

