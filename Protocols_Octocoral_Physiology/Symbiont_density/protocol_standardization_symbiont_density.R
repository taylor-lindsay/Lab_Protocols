# Title: 
# Author:
# Date: 
# Input files:
#  
# Output files: 
# 
# Notes

# Libraries ---------------------------------------------------------------

library(tidyverse)

# TL_octo

sym_path <- "~/Desktop/GITHUB/PR_OCTO/Practice/sym/TL_octo_Sym_D_Counts.csv"                         
meta_path <- "~/Desktop/GITHUB/PR_OCTO/Practice/TL_octo_practice_raw.csv"
output_path <- '~/Desktop/GITHUB/PR_OCTO/Practice/Results/TL_octo_results_sym.csv'

# Data -------------------------------------------------------------------

#raw symbiont data 
sym <- read.csv(paste(sym_path)) %>%
  select(sample_id,average_per_square,CV) %>%
  filter(average_per_square !="#DIV/0!") 

# Load homogenate volume
meta <- read_csv(paste(meta_path)) %>%
  select(sample_id,sym_weight) 

# standardize -------------------------------------------------------------------
# Join DF with homogenate 
sym <- full_join(sym, meta)

# Multiply chlorophyll by the homogenate volume and divide by surface area
sym2 <- sym %>%
  mutate(sym.mg = as.numeric(average_per_square)/ as.numeric(sym_weight))

sym_small <- sym2 %>%
  #filter(as.numeric(CV) <= 15)  %>% 
  select(sample_id,sym.mg)

# write the file 
write.csv(sym_small, paste(output_path))    

