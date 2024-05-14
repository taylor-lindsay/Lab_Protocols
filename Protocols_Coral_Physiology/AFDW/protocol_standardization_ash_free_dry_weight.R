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

AFDW_path <- "~/Desktop/GITHUB/TLPR21/AFDW/TLPR21_AFDW.csv"
meta_path <- "~/Desktop/GITHUB/TLPR21/TLPR21_Raw_Master.csv"
surface_path <-"~/Desktop/GITHUB/TLPR21/Surface_Area/TLPR21_Surface_Area.csv"
output_path <- "~/Desktop/GITHUB/TLPR21/AFDW/TLPR21_AFDW_Results.csv"

#AFDW data
AFDW <- read.csv(paste(AFDW_path)) %>%
  select(colony_id,sym_host,AFDW.g.ml)

# Load homogenate volume
vol <- read_csv(paste(meta_path)) %>%                                              #####
select(colony_id, airbrush_volume) %>%
  filter(!is.na(airbrush_volume))

# Load surface area
sa <- read_csv(paste(surface_path)) %>%                                #####
select(colony_id, surface_area) %>%
  filter(!is.na(surface_area))

# Merge & Edit Data -------------------------------------------------------

#separate sym and host data 


Sym <- AFDW %>%
  filter(.,sym_host=="SYM") %>%
  filter(.,AFDW.g.ml>0)%>%
  setNames(c("colony_id","sym_host","AFDW_sym")) %>%
  .[,c(1,3)]

Host <- AFDW %>%
  filter(.,sym_host=="HOST") %>%
  filter(.,AFDW.g.ml>0)%>%
  setNames(c("colony_id","sym_host","AFDW_host")) %>%
  .[,c(1,3)]

#merge two together
merged1 <- full_join(Sym,Host)
#merge with metadata 
AFDW_Merge <- full_join(sa,merged1, by="colony_id")
AFDW_Merge <- full_join(vol,AFDW_Merge, by="colony_id")

# Original data was g/ml

AFDW_fin <- AFDW_Merge %>%
  mutate(Sym_AFDW_mg.cm2 = ((as.numeric(AFDW_sym) * as.numeric(airbrush_volume) / as.numeric(surface_area))*1000)) %>%
  mutate(Host_AFDW_mg.cm2 = ((as.numeric(AFDW_host) * as.numeric(airbrush_volume) / as.numeric(surface_area))*1000))

# just the columns I want 
AFDW_small <- AFDW_fin %>%
  select(colony_id, Host_AFDW_mg.cm2, Sym_AFDW_mg.cm2)

ggplot(AFDW_small, aes(x=colony_id,y=Sym_AFDW_mg.cm2)) +
  geom_point()

# write the file 
write.csv(AFDW_small, paste(output_path))

