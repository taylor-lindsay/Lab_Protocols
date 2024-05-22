# Title: TL_Trans_AFDW.R
# Author: Taylor Lindsay
# Date: 04.09.2023
# Input files:
#  TL_Trans_AFDW.csv
# Output files: 
# 
# Notes 

# Packages & Data Import --------------------------------------------------

# Install Packages
library(tidyverse)
library(tidyr)
library(ggplot2)
library(dplyr)

AFDW_path <- "~/Desktop/GITHUB/TLPR21_2/AFDW/TLPR21_AFDW.csv"
meta_path <- "~/Desktop/GITHUB/TLPR21_2/TLPR21_Raw_Master.csv"
surface_path <-"~/Desktop/GITHUB/TLPR21_2/Surface_Area/TLPR21_Surface_Area.csv"
output_path <- "~/Desktop/GITHUB/TLPR21_2/AFDW/TLPR21_AFDW_Results.csv"

#AFDW data
AFDW <- read.csv(paste(AFDW_path)) %>%
  select(sample_id,sym_host,AFDW.g.ml)

# Load homogenate volume
vol <- read_csv(paste(meta_path)) %>%                                              #####
select(sample_id, airbrush_volume) %>%
  filter(!is.na(airbrush_volume))

# Load surface area
sa <- read_csv(paste(surface_path)) %>%                                #####
select(sample_id, surface_area) %>%
  filter(!is.na(surface_area))

# Merge & Edit Data -------------------------------------------------------

# filter negatives and remove duplicates 
AFDW2 <- AFDW %>%
  filter(., AFDW.g.ml > 0)%>%
  group_by(sample_id, sym_host) %>%
  summarise(AFDW_mean = mean(AFDW.g.ml))

#separate sym and host data 
Sym <- AFDW2 %>%
  filter(.,sym_host=="SYM") %>%
  setNames(c("sample_id","sym_host","AFDW_sym")) %>%
  select(c(sample_id, AFDW_sym))

Host <- AFDW2 %>%
  filter(.,sym_host=="HOST") %>%
  setNames(c("sample_id","sym_host","AFDW_host")) %>%
  select(c(sample_id, AFDW_host))

#merge two together
merged1 <- full_join(Sym,Host)
#merge with metadata 
AFDW_Merge <- full_join(sa,merged1, by="sample_id")
AFDW_Merge <- full_join(vol,AFDW_Merge, by="sample_id")

# Original data was g/ml

AFDW_fin <- AFDW_Merge %>%
  mutate(Sym_AFDW_mg.cm2 = ((as.numeric(AFDW_sym) * as.numeric(airbrush_volume) / as.numeric(surface_area))*1000)) %>%
  mutate(Host_AFDW_mg.cm2 = ((as.numeric(AFDW_host) * as.numeric(airbrush_volume) / as.numeric(surface_area))*1000))

# just the columns I want 
AFDW_small <- AFDW_fin %>%
  select(sample_id, Host_AFDW_mg.cm2, Sym_AFDW_mg.cm2)

ggplot(AFDW_small, aes(x=sample_id,y=Sym_AFDW_mg.cm2)) +
  geom_point()

# write the file 
write.csv(AFDW_small, paste(output_path))

