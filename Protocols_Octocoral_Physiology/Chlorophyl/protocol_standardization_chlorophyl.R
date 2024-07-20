# Title: CHL analysis file OCTOCORALS
# Author: Taylor Lindsay (updated from Putnam Lab)
# Date: 06.29.2024
# Input files:
    # CHL platemaps       
    # CHL platereader outputs 
    # Homogenate volume 
# Output files:
    # Results file (chl a & chl c2)
# Notes 


# Packages ----------------------------------------------------------------

#install.packages("dplyr")

# load packages
library(tidyverse)
library(dplyr)

# Import Data -------------------------------------------------------------

chl_path <- "~/Desktop/GITHUB/PR_OCTO/CHL/CHL_platemap_results/"                         
meta_path <- "~/Desktop/GITHUB/PR_OCTO/widths_weights/TL_Octo_Phys_prep.csv"
output_path <- '~/Desktop/GITHUB/PR_OCTO/CHL/TL_octo_results_chl.csv'

#chl_path <- "~/Desktop/GITHUB/PR_OCTO/Practice/chl/"                         
#meta_path <- "~/Desktop/GITHUB/PR_OCTO/Practice/TL_octo_practice_raw.csv"
#output_path <- '~/Desktop/GITHUB/PR_OCTO/Practice/Results/TL_octo_results_chl.csv'

# Load sample weight
meta <- read_csv(paste(meta_path)) %>%                                           
  select(sample_id, chl_weight)

# Define function to read in chl data
read_chl <- function(file) {
  chl_data <- read_csv(file, skip = 39, n_max = 39) %>%
    #select(-1) %>%
    magrittr::set_colnames(c("row", 1:12, "wavelength")) %>%
    fill(row) %>%
    gather("col", "absorbance", -wavelength, -row) %>%
    unite("well", c(row, col), sep = "")
}

# List chlorophyll data files
all_chl_files <- list.files(path = chl_path, pattern = "*.csv")          # List all files in directory
chl_platemaps <- list.files(path = chl_path, pattern = "platemap")       # List platemap files
chl_data_files <- setdiff(all_chl_files, chl_platemaps)                  # List absorbance data files

# Read in all files into tibble
df <- tibble(file = chl_data_files) %>%
  mutate(platemap = map(file, ~ read_csv(paste0(chl_path, tools::file_path_sans_ext(.), "_platemap.csv"))),
         chl_data = map(file, ~ read_chl(paste0(chl_path, .))))

# Merge platemap and data for each plate
df <- df %>%
  mutate(merged = map2(platemap, chl_data, ~ right_join(.x, .y))) %>% 
  unnest(merged) %>%
  filter(!is.na(sample_id)) %>%
  select(c(file, well, sample_id,wavelength,absorbance)) %>%
  ungroup()

# Calculate Chl Concentrations --------------------------------------------

# average all technical replicates for each plate/sample/wavelength, including all acetone blanks together (per plate)
df2 <- df %>%                    
  group_by(sample_id,wavelength,file) %>%
  summarize(mean_abs = mean(absorbance)) %>%
  spread(wavelength, mean_abs)

# get the acetone blank 750 absorbace for each file (i.e., plate), and subtract from 630 and 663 values for each sample
df2 <- df2 %>%
  group_by(file) %>%
  mutate(blank750 = `Chl:750`[sample_id == "BK"]) %>%
  ungroup() %>%
  mutate(adj630 = `Chl:630` - blank750,
         adj663 = `Chl:663` - blank750)

# calculate chla and chlc2 values based on equations from Jeffrey and Humphrey 1975
# units µg/mg
#path length adjustment = 0.6, which is based on the 200ul path length 

df2 <- df2 %>%
  mutate(chla.ug.ml = (11.43 * adj663)/0.6 - (0.64 * adj630)/0.6,
         chlc2.ug.ml = (27.09 * adj630)/0.6 - (3.63 * adj663)/0.6)

# Normalize to Surface Area -----------------------------------------------

# Join DF with homogenate 
chl <- full_join(df2, meta)

# Multiply chlorophyll by the homogenate volume and divide by surface area
chl <- chl %>%
  mutate(chla.ug.mg = chla.ug.ml / as.numeric(chl_weight),
         chlc2.ug.mg = chlc2.ug.ml / as.numeric(chl_weight))

# remove blanks and NAs
chl <- filter(chl, !sample_id %in% c("NA", "BK"))

# write chlorophyll data to file
final <- chl %>%
  select(sample_id, chla.ug.mg, chlc2.ug.mg) %>%
  filter(!is.na(chla.ug.mg))%>%
  filter(!is.na(chlc2.ug.mg)) 

final %>%
  write.csv(., paste(output_path))                                
