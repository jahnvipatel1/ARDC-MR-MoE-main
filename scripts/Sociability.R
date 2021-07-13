library(devtools)
library(tidyverse)
library(TwoSampleMR)

#EXPOSURE DATA ---------------------------------------------
sociability.file <- "~/Desktop/Bralten2021sociability.chrall.CPRA_b37 2.tsv"
sociability_dat <- read_exposure_data(
  filename = sociability.file,
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
sociability_dat<- filter(physicalactivity_dat, pval.exposure < 0.00000005) %>%
  distinct(SNP, .keep_all = TRUE)

#CLUMP DATA
sociability_dat <- clump_data(sociability_dat)

#OUTCOME DATA ---------------------------------------------
AD.file <- "~/Desktop/Kunkle2019load_stage123.chrall.CPRA_b37.tsv"
AD_dat <- read_outcome_data(
  snps = sociability_dat$SNP,
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
  pos_col = "POS",
  phenotype_col = "TRAIT"
)

#HARMONISE DATA --------------------------------------------------------
sociability.dat <- harmonise_data(sociability_dat,AD_dat)

#MR-MoE-----------------------------------------------------------------
load("/Users/jahnvipatel/Downloads/rf.rdata")
res.soc <- mr_wrapper(sociability.dat)
res_moe.soc <- mr_moe(res.soc,rf) 