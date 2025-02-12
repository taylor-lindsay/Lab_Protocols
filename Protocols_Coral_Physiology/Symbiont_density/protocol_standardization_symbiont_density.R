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

# TLAP
sym_path <- "~/Desktop/GITHUB/TL_Astrangia/TLAP_Phys_CHL_sym/TLAP_Sym_D_Counts.csv"
meta_path <- "~/Desktop/GITHUB/TL_Astrangia/TLAP_Raw_Data/TLAP_ALL_Results.csv"
output_path <- "~/Desktop/GITHUB/TL_Astrangia/TLAP_Raw_Data/TLAP_Sym_Standardized.csv"

#TLPR21
sym_path <- '~/Desktop/GITHUB/TLPR21_2/Sym_density/TLPR21_Sym_D_Counts_Passed_Counts.csv'
meta_path <- '~/Desktop/GITHUB/TLPR21_2/TLPR21_Raw_Master.csv'
output_path <- '~/Desktop/GITHUB/TLPR21_2/Sym_density/TLPR21_Sym_Results.csv'

# TL_Trans
sym_path <- "~/Desktop/GITHUB/TL_Trans/Sym/TL_Trans_Sym_D_Counts.csv"
meta_path <- "~/Desktop/GITHUB/TL_Trans/TL_Trans_Raw_Master.csv"
output_path <- "~/Desktop/GITHUB/TL_Trans/Sym/TL_TRANS_Sym_Standardized.csv"

# Data -------------------------------------------------------------------

#raw symbiont data 
sym <- read.csv(paste(sym_path)) %>%
  select(sample_id,average_per_square) %>%
  filter(average_per_square !="#DIV/0!") 

# Load homogenate volume
meta <- read_csv(paste(meta_path)) %>%                                              #####
select(sample_id, airbrush_volume,surface_area) %>%
  filter(!is.na(airbrush_volume)) 

# standardize -------------------------------------------------------------------
# Join DF with homogenate 
sym <- full_join(sym, meta)

# Multiply chlorophyll by the homogenate volume and divide by surface area
sym2 <- sym %>%
  filter(!is.na(surface_area)) %>%
  mutate(sym.cm2 = (as.numeric(average_per_square)*10000) * as.numeric(airbrush_volume) / as.numeric(surface_area))

sym_small <- sym2 %>%
  #filter(as.numeric(CV) <= 15)  %>% 
  select(sample_id,sym.cm2)

# write the file 
write.csv(sym_small, paste(output_path))    

