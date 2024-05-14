## # Title: protocol_standardization_lipids.R
# Author: TL
# Date: 05.14.2024
# Input files: 
#  Raw lipids Data & Metadata 
# Output files: 
#  Standardized lipids data 
# Notes

# Libraries & Data  ---------------------------------------------------------------
library(tidyverse)
library(dplyr)
library(ggplot2)

theme_set(theme_bw())

# update the paths of the data as needed 
raw_path <- "~/Desktop/GITHUB/TL_Astrangia/Lipids/TLAP_Lipid_Extraction.csv"
meta_path <- "~/Desktop/GITHUB/TL_Astrangia/Raw_Data/AP23_ALL_Results.csv"
output_path <-"~/Desktop/GITHUB/TL_Astrangia/Raw_Data/TLAP_lipids_Standardized.csv"

# Load raw lipids data
raw <- read_csv(paste(raw_path)) %>%
  select(c(sample_id, AFDW, percent_lipids, total_lipids)) %>%
  filter(total_lipids >= 0)

# Load metadata (surface area & homogenate volume)
meta <- read_csv(paste(meta_path)) %>%
  select(c(sample_id, airbrush_volume, surface_area))

# standardize -------------------------------------------------------------------

# find the mean of sample replicates 
raw_means <- raw %>%
  group_by(sample_id) %>%
  summarise_at(vars(AFDW:total_lipids), mean)

# Join DF with homogenate 
lipids <- full_join(raw_means, meta)

# Multiply lipids & AFDW by the homogenate volume and divide by surface area
lipids_results <- lipids %>%
  filter(!is.na(surface_area)) %>%
  mutate(AFDW.cm2 = as.numeric(AFDW) * as.numeric(airbrush_volume) / as.numeric(surface_area)) %>%
  mutate(lipids.mg.cm2 = as.numeric(total_lipids*1000) * as.numeric(airbrush_volume) / as.numeric(surface_area)) %>%
  select(c(sample_id, percent_lipids, AFDW.cm2, lipids.mg.cm2)) %>%
  filter(!is.na(lipids.mg.cm2))

# write the file 
write.csv(lipids_results, paste(output_path))    




