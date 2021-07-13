## Extract exposure instruments from outcome gwas

###Command Line Arguments 
args = commandArgs(trailingOnly = TRUE) # Set arguments from the command line
exposure.summary = args[1] # Exposure summary statistics
outcome.summary = args[2] # Outcome Summary statistics
out = args[3]

###Load Packages
library(TwoSampleMR)
load("rf.rdata") #can find on Mendelian Randomization Dropbox

###Read in SNPs
message("READING IN EXPOSURE \n")
exposure_dat <- read_tsv(exposure.summary)
