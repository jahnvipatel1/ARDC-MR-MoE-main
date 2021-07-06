library(devtools)
library(tidyverse)
library(TwoSampleMR)
library(ggnewscale)
library(ggplot2)

#EXPOSURE DATA ---------------------------------------------
education.file <- "~/Desktop/Lee2018educ.chrall.CPRA_b37.tsv"
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
education.dat <- harmonise_data(education_dat,AD_dat)

#MR-MoE-----------------------------------------------------------------
load("/Users/jahnvipatel/Downloads/rf.rdata")
res <- mr_wrapper(education.dat)
res_moe <- mr_moe(res,rf)

#PREFORM MR ------------------------------------------------------------
res.mr <- mr(education.dat) #each MR method for each combination of exposure-outcome traits
generate_odds_ratios(res.mr)

#MAKE PLOTS ------------------------------------------------------------
#create a scatter plot
res.mr <- mr(education.dat, 
             method_list=c(
               "mr_ivw",
               "mr_weighted_mode",
               "mr_weighted_median",
               "mr_simple_median",
               "mr_simple_mode",
               "mr_egger_regression")
)
p1 <- mr_scatter_plot(res.mr,education.dat)
p1[[1]]
ggsave(p1[[1]], file="educationscatterplot.png", width=7, height=8)
view(education.dat)
view(p1$NCFz7e.oQE1lc$data)

mr_scatter_plot2 <- function(mrdat,res){
  
  message("Plotting Scatters: ", mrdat$exposure[1], " - ", mrdat$outcome[1])
  
  ## reoriante effect direction of negative exposure values
  
  d <- mrdat %>%
    
    mutate(beta.outcome = ifelse(beta.exposure < 0, beta.outcome * -1, beta.outcome),
           
           beta.exposure = ifelse(beta.exposure < 0, beta.exposure * -1, beta.exposure))
  
  ## Make scatter plot
  
  ggplot(data = d, aes(x = beta.exposure, y = beta.outcome)) +
    
    geom_errorbar(aes(ymin = beta.outcome - se.outcome, ymax = beta.outcome + se.outcome),
                  
                  colour = "grey", width = 0) +
    
    geom_errorbarh(aes(xmin = beta.exposure - se.exposure, xmax = beta.exposure + se.exposure),
                   
                   colour = "grey", height = 0) +
    
    geom_point(aes(colour = "!mrpresso_keep")) +
    
    scale_colour_manual(values = c("black", "#CD534CFF")) +
    
    new_scale_color() +
    
    geom_abline(data = res, aes(intercept = a, slope = b, color = method, linetype = method), show.legend = TRUE) +
    
    scale_color_brewer(palette = "Set3") +
    
    labs(x = paste("SNP effect on\n", d$exposure[1]),
         
         y = paste("SNP effect on\n", d$outcome[1])) +
    
    theme_bw() +
    
    theme(legend.position = "bottom",
          
          legend.direction = "horizontal",
          
          text = element_text(size=8)) +
    
    guides(linetype = guide_legend(nrow = 1),
           =
           colour_new = FALSE)
  
}
a<-c(0,0,-0.106,0.0842,0,0,0,0,0,0,0,0,0,0,0,0,0.0421,0.0842,0.0421,0,0.0853,0,0,0,0,0,0,0,0,0,0.0853,0,0,0,0,0,0,0,-0.0106,0,0,0,0,0)
res_moe$NCFz7e.oQE1lc$estimates$a=a

p <- mr_scatter_plot2(mrdat = education.dat, res = res_moe$NCFz7e.oQE1lc$estimates)
p + facet_wrap(vars(method))

#create a forest plot
res_single_educ <- mr_singlesnp(education.dat, 
                              all_method=c(
                                "mr_ivw",
                                "mr_weighted_mode",
                                "mr_weighted_median",
                                "mr_simple_median",
                                "mr_simple_mode",
                                "mr_egger_regression")) 
p2 <- mr_forest_plot(res_single_educ)
p2[[1]]
ggsave(p2[[1]], file="educationforestplot.png", width=7, height=20)

#create a funnel plot
res_single_educ <- mr_singlesnp(
  education.dat,
  all_method=c("mr_ivw",
               "mr_weighted_mode",
               "mr_weighted_median",
               "mr_simple_median",
               "mr_simple_mode",
               "mr_egger_regression")
)
p3 <- mr_funnel_plot(res_single_educ)
p3[[1]]
ggsave(p3[[1]], file="educationfunnelplot.png", width=7, height=7)

#generate report
mr_report(education.dat) 

#generate spreadsheet
library(writexl)
write_xlsx(res_moe$w8C1MU.DWrQTt,"\\Desktop\\Education MR-MOE.xlsx")