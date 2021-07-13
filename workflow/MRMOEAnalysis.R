#!/usr/bin/Rscript
## ========================================================================== ##
## Run MR-MOE Analysis
## ========================================================================== ##

### ===== Command Line Arguments ===== ##
args = commandArgs(trailingOnly=TRUE)

exposure.summary = args[1] #exposure summary statistics
p.threshold = as.numeric(args[2])
outcome.summary = args[3] #outcome summary statistics
out.file = args[4] #specify output file

### ===== Load packages ===== ###
library(devtools)
library(tidyverse)
library(TwoSampleMR)
load("rf.rdata") #can find on Mendelian Randomization Dropbox

### ===== Read in Exposure Data ===== ###
message("\n READING IN EXPOSURE \n")

exposure_dat <- read_exposure_data(
  filename = exposure.summary,
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

### ===== Filter Exposure Data ===== ###
exposure_dat<- filter(exposure_dat, pval.exposure < p.threshold) %>%
  distinct(SNP, .keep_all = TRUE)

### ===== Clump Exposure ===== ###
message("\n CLUMPING EXPOSURE SNPS \n")

exposure_dat <- clump_data(exposure_dat)

### ===== Read in Outcome Data ===== ###
message("\n READING IN OUTCOME \n")

outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = outcome.summary,
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

### ===== Harmonization ===== ###
message("\n Begining Harmonization \n")

harmonized.dat <- harmonise_data(expsoure_dat,outcome_dat)

### ===== Run MR-MOE ===== ###
message("\n Running MR-MOE \n")
res <- mr_wrapper(harmonized.dat)
res_moe <- mr_moe(res,rf)

###Write out results from MR-MOE
message("\n Writing Out MR-MOE Results \n")
write_csv(res_moe, out.file)