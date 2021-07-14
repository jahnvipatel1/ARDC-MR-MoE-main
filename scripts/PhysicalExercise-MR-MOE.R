library(devtools)
library(tidyverse)
library(TwoSampleMR)
library(ggplot2)
library(ggnewscale)
install.packages("openxlsx",dependencies = TRUE)
library(openxlsx)
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

res <- do.call(rbind, res_moe.pa)
class(res$estimates)
write_xlsx(res_moe.pa[[1]],"\\Desktop\\Physical Exercise MR-MOEahhhh.xlsx")
write.csv(res_moe,glue("~/{out.file}/{exposure.summary}_{outcome.summary}.csv"))
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
res_moe.pa$gEOyNc.6cf5D1$estimates$a=b

q <- mr_scatter_plot2(mrdat=physicalactivity.dat,res=res_moe.pa$gEOyNc.6cf5D1$estimates)
q[[1]]
q$plot_env$d$remove=c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE)
q$plot_env$d$action=c(2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,1)
z <- q$plot_env$d[-c(10,13)]
z <- z[-c(10,13),]
q$plot_env$d=z
q + facet_wrap(vars(selection)) + theme(legend.position = 'top') 

scatter.tophits.pa <- mr_scatter_plot(res.pa,res.filter2.pa)
scatter.tophits.pa.a <- scatter.tophits.pa[[1]] + theme_bw() + theme(legend.position = 'bottom')

scatter.steiger.pa <- mr_scatter_plot(res.pa,(slice(res.filter2.pa,-c(13,18))))
scatter.steigera.pa <- scatter.steiger.pa[[1]] + theme_bw() + theme(legend.position = 'bottom') #steiger

scatter.both.pa <- mr_scatter_plot(res.pa,(slice(res.filter2.pa,-c(13,18))))
scatter.botha.pa <- scatter.both.pa[[1]] + theme_bw() + theme(legend.position = 'bottom') #both

ggpubr::ggarrange(scatter.tophits.pa.a, scatter.steigera.pa, scatter.botha.pa,common.legend = T)
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
res_single.pa <- mr_singlesnp(
  physicalactivity.dat,
  all_method=c("mr_ivw",
    "mr_weighted_mode",
    "mr_simple_median",
    "mr_egger_regression")
  )
p3 <- mr_funnel_plot(res_single)
p3[[1]]
ggsave(p3[[1]], file="physicalactivityfunnelplot.png", width=7, height=7)

res.filter.pa <- res_single.pa %>% slice(1:19)
res_moe2.pa <- subset(res_moe.pa$gEOyNc.6cf5D1$estimates,select = -c(nsnp,ci_low,ci_upp,steiger_filtered,outlier_filtered,selection,method2,MOE))
res_methods.pa <- paste0("All - ",res_moe2.pa$method)
res_moe2.pa <- subset(res_moe2.pa,select=-c(method))
res_moe2.pa <- cbind(res_methods.pa,res_moe2.pa)
colnames(res_moe2.pa)[1] <- "SNP"

exposure.pa <- rep(c("MVPA"),times=c(33))
outcome.pa <- rep(c("AD"),times=c(33))
samplesize.pa <- rep(c(63926),times=c(33))
res_moe3.pa<- cbind(exposure.pa,outcome.pa,samplesize.pa)
colnames(res_moe3.pa)[1] <- "exposure"
colnames(res_moe3.pa)[2] <- "outcome"
colnames(res_moe3.pa)[3] <- "samplesize"

res_moe_filter.pa <- merge(res_moe2.pa,res_moe3.pa,by=0)
res_moe_filter.pa <- subset(res_moe_filter.pa,select = -c(Row.names))
colnames(res_moe_filter.pa)[6] <- "p"
res.filter2.pa <- rbind(res.filter.pa,res_moe_filter.pa)

tophits.pa <- mr_funnel_plot(res.filter2.pa)
tophits.pa.a <- tophits.pa[[1]] + theme_bw() + theme(legend.position = 'bottom')


steiger.pa <- mr_funnel_plot(slice(res.filter2.pa,-c(13,18)))
steigera.pa <- steiger.pa[[1]] + theme_bw() + theme(legend.position = 'bottom') #steiger

both.pa <- mr_funnel_plot(slice(res.filter2.pa,-c(13,18)))
botha.pa <- both.pa[[1]] + theme_bw() + theme(legend.position = 'bottom') #both

ggpubr::ggarrange(tophits.pa.a, steigera.pa, botha.pa,common.legend = T)

#generate report
mr_report(physicalactivity.dat) 

#generate spreadsheet
library(writexl)
write_xlsx(res_moe.pa$`0b3TOb.Rin5Eq`,"\\Desktop\\Physical Exercise MR-MOE.xlsx")
