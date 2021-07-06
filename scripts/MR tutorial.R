library(devtools)
install.packages("remotes")
remotes::install_github("MRCIEU/TwoSampleMR")
library(tidyverse)
library(TwoSampleMR)

#EXPOSURE DATA ---------------------------------------------
education.file <- "~/Desktop/Lee2018educ.chrall.CPRA_b37 2.tsv"
education_dat <- read_exposure_data(
    filename = education.file,
    sep = "\t",
    snp_col = "DBSNP_ID",
    beta_col = "BETA",
    se_col = "SE",
    effect_allele_col = "ALT",
    other_allele_col = "REF",
    eaf_col = "AF",
    pval_col = "P",
    samplesize_col = "N",
    chr_col = "CHROM",
    pos_col = "POS",
    phenotype_col = "TRAIT"
  )

#OUTCOME DATA ---------------------------------------------
AD.file <- "~/Desktop/Kunkle2019load_stage123.chrall.CPRA_b37.tsv"
AD_dat <- read_outcome_data(
    snps = education_dat$SNP,
    filename = AD.file,
    sep = "\t",
    snp_col = "DBSNP_ID",
    beta_col = "BETA",
    se_col = "SE",
    effect_allele_col = "ALT",
    other_allele_col = "REF",
    eaf_col = "AF",
    pval_col = "P",
    samplesize_col = "N",
    chr_col = "CHROM",
    pos_col = "POS"
  )

#HARMONISE DATA --------------------------------------------------------
dat <- harmonise_data(education_dat,AD_dat)

#PREFORM MR ------------------------------------------------------------
res.mr <- mr(dat)
p1 <- mr_scatter_plot(res,dat)
res.mr
mr_method_list()
mr(dat, method_list=c("mr_egger_regression,","mr_ivw"))

#MR-MoE-----------------------------------------------------------------
load("rf.rdata")
res <- mr_wrapper(dat)
res_moe <- mr_moe(res,rf)