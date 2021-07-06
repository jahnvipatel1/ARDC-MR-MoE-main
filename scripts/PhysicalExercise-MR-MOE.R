library(devtools)
library(tidyverse)
library(TwoSampleMR)

#EXPOSURE DATA ---------------------------------------------
physicalactivity.file <- "~/Desktop/Klimentidis2018mvpa.chrall.CPRA_b37.tsv"
physicalactivity_dat <- read_exposure_data(
  filename = physicalactivity.file,
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
physicalactivity_dat<- filter(physicalactivity_dat, pval.exposure < 0.00000005) %>%
  distinct(SNP, .keep_all = TRUE)

#CLUMP DATA
physicalactivity_dat <- clump_data(physicalactivity_dat)

#OUTCOME DATA ---------------------------------------------
AD.file <- "~/Desktop/Kunkle2019load_stage123.chrall.CPRA_b37.tsv"
AD_dat <- read_outcome_data(
  snps = physicalactivity_dat$SNP,
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
physicalactivity.dat <- harmonise_data(physicalactivity_dat,AD_dat)

#MR-MoE-----------------------------------------------------------------
load("/Users/jahnvipatel/Downloads/rf.rdata")
res.pa <- mr_wrapper(physicalactivity.dat)
res_moe.pa <- mr_moe(res.pa,rf) 

#PREFORM MR ------------------------------------------------------------
res.mr <- mr(physicalactivity.dat) #each MR method for each combination of exposure-outcome traits
generate_odds_ratios(res.mr)

#MAKE PLOTS ------------------------------------------------------------
#create a scatter plot
res.pa <- mr(physicalactivity.dat, 
          method_list=c(
            "mr_ivw",
            "mr_weighted_mode",
            "mr_weighted_median",
            "mr_simple_median",
            "mr_simple_mode",
            "mr_egger_regression")
          )
p1 <- mr_scatter_plot(res.pa, physicalactivity.dat)
p1[[1]]
ggsave(p1[[1]], file="physicalactivityscatterplot.png", width=7, height=7)

b<-c(0,0,0,0,-2.842025,0,0,0,0,0,0,0,-2.842025,0,0,0,0,0,0,0,-2.842025,0,0,0,0,0,0,0,-2.842025,0,0,-9.373453,-9.373453)
res_moe.pa$lXWUfu.VcdeTs$estimates$a=b

q <- mr_scatter_plot2(mrdat=physicalactivity.dat,res=res_moe.pa$lXWUfu.VcdeTs$estimates)
q + facet_wrap(vars(method))
#create a forest plot
res_single_pa <- mr_singlesnp(physicalactivity.dat, 
                           all_method=c(
                             "mr_ivw",
                             "mr_weighted_mode",
                             "mr_weighted_median",
                             "mr_simple_median",
                             "mr_simple_mode",
                             "mr_egger_regression")) 
p2 <- mr_forest_plot(res_single_pa)
p2[[1]]
ggsave(p2[[1]], file="physicalactivityforestplot.png", width=7, height=7)

#create a funnel plot
res_single <- mr_singlesnp(
  physicalactivity.dat,
  all_method=c("mr_ivw",
    "mr_weighted_mode",
    "mr_simple_median",
    "mr_egger_regression")
  )
p3 <- mr_funnel_plot(res_single)
p3[[1]]
ggsave(p3[[1]], file="physicalactivityfunnelplot.png", width=7, height=7)

#generate report
mr_report(physicalactivity.dat) 

#generate spreadsheet
library(writexl)
write_xlsx(res_moe.pa$`0b3TOb.Rin5Eq`,"\\Desktop\\Physical Exercise MR-MOE.xlsx")
