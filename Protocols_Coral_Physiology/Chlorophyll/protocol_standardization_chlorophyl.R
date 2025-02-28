# Title: CHL analysis file
# Author: Taylor Lindsay (updated from Putnam Lab)
# Date: 01.16.2022
# Input files:
    # CHL platemaps       
    # CHL platereader outputs 
    # Homogenate volume 
    # Surface area 
# Output files:
    # Results file (chl a & chl c2)
# Notes 


# Packages ----------------------------------------------------------------

## install packages if you don't already have them
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("plotrix")) install.packages("plotrix")

# load packages
library(plotrix)
library(tidyverse)

# Import Data -------------------------------------------------------------

# TLPR21 Path 
chl_path <- "~/Desktop/GITHUB/TLPR21_2/CHL/TLPR21_CHL_raw_files/"                         
surface_path <- "~/Desktop/GITHUB/TLPR21_2/Surface_Area/TLPR21_Surface_Area.csv"
vol_path <- "~/Desktop/GITHUB/TLPR21_2/TLPR21_Raw_Master.csv"
output_path <- '~/Desktop/GITHUB/TLPR21_2/CHL/TLPR21_Results_CHL.csv'

# TLPR21 Path 
chl_path <- "~/Desktop/Matias_samples/Chlorophyll/raw/"                         
surface_path <- "~/Desktop/Matias_samples/MGPR21_Raw.csv"
vol_path <- "~/Desktop/Matias_samples/MGPR21_Raw.csv"
output_path <- '~/Desktop/Matias_samples/Chlorophyll/MGPR21_Results_CHL.csv'

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
  mutate(merged = map2(platemap, chl_data, ~ right_join(.x, .y)))

# Load homogenate volume
vol <- read_csv(paste(vol_path)) %>%                                           
  select(sample_id, airbrush_volume) %>%
  filter(!is.na(airbrush_volume))

# Load surface area
sa <- read_csv(paste(surface_path)) %>%                              
  select(sample_id, surface_area) %>%
  filter(!is.na(surface_area))

# Calculate Chl Concentrations --------------------------------------------


# average all technical replicates for each plate/sample/wavelength, including all acetone blanks together (per plate)
df <- df %>%
  unnest(merged) %>%
  filter(!is.na(sample_id)) %>%                         # remove empty wells (colony_id is NA)
  group_by(file, sample_id, wavelength) %>%
  summarise(n = n(), mean_abs = mean(absorbance)) %>%
  spread(wavelength, mean_abs)

# get the acetone blank 750 absorbace for each file (i.e., plate), and subtract from 630 and 663 values for each sample
df <- df %>%
  group_by(file) %>%
  mutate(blank750 = `Chl:750`[sample_id == "BK"]) %>%
  ungroup() %>%
  mutate(adj630 = `Chl:630` - blank750,
         adj663 = `Chl:663` - blank750)

# calculate chla and chlc2 values based on equations from Jeffrey and Humphrey 1975
# units µg/ml
#path length adjustment = 0.6 

df <- df %>%
  mutate(chla.ug.ml = (11.43 * adj663)/0.6 - (0.64 * adj630)/0.6,
         chlc2.ug.ml = (27.09 * adj630)/0.6 - (3.63 * adj663)/0.6)

#previous, with no pathlength adjustment
#df <- df %>%
#mutate(chla.ug.ml = (11.43 * adj663) - (0.64 * adj630),
#chlc2.ug.ml = (27.09 * adj630) - (3.63 * adj663))



# Normalize to Surface Area -----------------------------------------------

# Join DF with homogenate 
chl <- full_join(df, vol)

# Join df with surface area 
chl <- full_join(chl, sa)

# Multiply chlorophyll by the homogenate volume and divide by surface area
chl <- chl %>%
  filter(!is.na(surface_area)) %>%
  mutate(chla.ug.cm2 = chla.ug.ml * as.numeric(airbrush_volume) / as.numeric(surface_area),
         chlc2.ug.cm2 = chlc2.ug.ml * as.numeric(airbrush_volume) / as.numeric(surface_area))

# remove blanks and NAs
chl <- filter(chl, !sample_id %in% c("NA", "BK"))

# write chlorophyll data to file
chl %>%
  select(sample_id, chla.ug.cm2, chlc2.ug.cm2) %>%
  #mutate(timepoint="FEB")%>%
  filter(!is.na(chla.ug.cm2))%>%
  filter(!is.na(chlc2.ug.cm2))%>%
  write.csv(., paste(output_path))                                
