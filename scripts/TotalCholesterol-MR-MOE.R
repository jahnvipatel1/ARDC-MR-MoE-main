library(devtools)
library(tidyverse)
library(TwoSampleMR)

#EXPOSURE DATA ---------------------------------------------
tcholesterol.file <- "~/Desktop/Willer2013tc.chrall.CPRA_b37.tsv"
tcholesterol_dat <- read_exposure_data(
  filename = tcholesterol.file,
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
tcholesterol_dat<- filter(tcholesterol_dat, pval.exposure < 0.00000005) %>%
  distinct(SNP, .keep_all = TRUE)

#CLUMP DATA
tcholesterol_dat <- clump_data(tcholesterol_dat)

#OUTCOME DATA ---------------------------------------------
AD.file <- "~/Desktop/Kunkle2019load_stage123.chrall.CPRA_b37.tsv"
AD_dat <- read_outcome_data(
  snps = tcholesterol_dat$SNP,
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
tcholesterol.dat <- harmonise_data(tcholesterol_dat,AD_dat)

#MR-MoE-----------------------------------------------------------------
load("/Users/jahnvipatel/Downloads/rf.rdata")
res.tc <- mr_wrapper(tcholesterol.dat)
res_moe.tc <- mr_moe(res.tc,rf) 

#PREFORM MR ------------------------------------------------------------
res.mr <- mr(tcholesterol.dat) #each MR method for each combination of exposure-outcome traits
generate_odds_ratios(res.mr)

#MAKE PLOTS ------------------------------------------------------------
#create a scatter plot
res.tc <- mr(tcholesterol.dat, 
             method_list=c(
               "mr_ivw",
               "mr_weighted_mode",
               "mr_weighted_median",
               "mr_simple_median",
               "mr_simple_mode",
               "mr_egger_regression")
)
p1 <- mr_scatter_plot(res.tc, tcholesterol.dat)
p1[[1]]
ggsave(p1[[1]], file="totalcholesterolscatterplot.png", width=7, height=7)

d<-c(0,0,0,-0.9612088,-0.9612088,0,0,0,0,0,0,0,0,0,0.428420,0.428420,0,0,0,0,0,0,0,0,0,0,-1.2168229,-1.2168229,0,0,0,0,0,0,0,0,0,0.3203365,0.3203365,0,0,0,0,0)
res_moe.tc$hbMKxk.V9qffW$estimates$a=d

z <- mr_scatter_plot2(mrdat=tcholesterol.dat,res=res_moe.tc$hbMKxk.V9qffW$estimates)
z + facet_wrap(vars(method))

#create a forest plot
res_single_tc <- mr_singlesnp(tcholesterol.dat, 
                              all_method=c(
                                "mr_ivw",
                                "mr_weighted_mode",
                                "mr_weighted_median",
                                "mr_simple_median",
                                "mr_simple_mode",
                                "mr_egger_regression")) 
p2 <- mr_forest_plot(res_single_tc)
p2[[1]]
ggsave(p2[[1]], file="totalcholesterolforestplot.png", width=7, height=9)

#create a funnel plot
res_single <- mr_singlesnp(
  tcholesterol.dat,
  all_method=c("mr_ivw",
               "mr_simple_mode",
               "mr_simple_median",
               "mr_egger_regression")
)
p3 <- mr_funnel_plot(res_single)
p3[[1]]
ggsave(p3[[1]], file="totalcholesterolfunnelplot.png", width=7, height=10)

#generate report
mr_report(tcholesterol.dat) 

#generate spreadsheet
library(writexl)
write_xlsx(res_moe.tc$Izjty2.r44nmc,"\\Desktop\\Total Cholesterol MR-MOE.xlsx")