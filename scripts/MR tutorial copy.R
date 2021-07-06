library(devtools)
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

#filter 
education_dat<- filter(education_dat, pval.exposure < 0.00000005) %>%
  distinct(SNP, .keep_all = TRUE)

#CLUMP DATA
education_dat <- clump_data(education_dat)

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
res.mr <- mr(dat) #each MR method for each combination of exposure-outcome traits
generate_odds_ratios(res.mr)

p1 <- mr_scatter_plot(res.mr,dat) #create a scatter plot
p1[[1]]
ggsave(p1[[1]], file="educationvsAD.pdf", width=7, height=7)

mr_method_list() #list of all the methods
mr(dat, method_list=c("mr_ivw")) #just ivw
mr_allmethods(dat,method=all)

res_single <- mr_singlesnp(dat) #create a forest plot
p2 <- mr_forest_plot(res_single)
p2[[1]]

res_single <- mr_singlesnp(dat, all_method=c("mr_ivw")) # forest plot with different methods
p2 <- mr_forest_plot(res_single)
p2[[1]]

mr_report(dat) #generate report

#MR-MoE-----------------------------------------------------------------
load("/Users/jahnvipatel/Downloads/rf.rdata")
res <- mr_wrapper(dat)
res_moe <- mr_moe(res,rf)
res_moe.estimates <- res_moe$w8C1MU.DWrQTt$estimates
res_moe.heterogeneity <-res_moe$w8C1MU.DWrQTt$heterogeneity
res_moe.egger.intercepts <-res_moe$w8C1MU.DWrQTt$directional_pleiotropy

mr.education <- mr(dat, method_list=c("mr_ivw_fe","mr_simple_median","mr_egger_regression","mr_simple_mode","mr_ivw_mre","mr_weighted_mode")) #has highest MOE AUROC score
oddsratios <- generate_odds_ratios(mr.education)