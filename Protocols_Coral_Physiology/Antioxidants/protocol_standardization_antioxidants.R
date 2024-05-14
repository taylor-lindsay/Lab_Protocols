#title: "Total Antioxidant Capacity (TAC) analysis"


# load packages
library(tidyverse)
library(broom)
library(ggpubr) 
library(ggplot2)

# write function to read files 
read_tac <- function(file) {
  tac_data <- read_csv(file, skip = 30, n_max = 24) %>%
    #select(-1) %>%
    magrittr::set_colnames(c("row", 1:12, "wavelength")) %>%
    #fill(row) %>%
    gather("col", "absorbance", -wavelength, -row) %>%
    unite("well", c(row, col), sep = "")
}

# List data files (taken from CHL sheet, just ignore that part)
chl_path <- "~/Desktop/GITHUB/TL_Astrangia/Antioxidants/AP23_Antioxidants/"                        # Path to chlorophyll data directory     #####
#chl_path <- "~/Desktop/GITHUB/TL_Astrangia/Antioxidants/Antioxidant_practice/data/"                        # Path to chlorophyll data directory     #####
all_chl_files <- list.files(path = chl_path, pattern = "*.csv")          # List all files in directory
chl_platemaps <- list.files(path = chl_path, pattern = "platemap")       # List platemap files
chl_data_files <- setdiff(all_chl_files, chl_platemaps)                  # List absorbance data files

# Read in all files into tibble
df <- tibble(file = chl_data_files) %>%
  mutate(platemap = map(file, ~ read_csv(paste0(chl_path, tools::file_path_sans_ext(.), "_platemap.csv"))), 
         tac_data = map(file, ~ read_tac(paste0(chl_path, .))))
df2 <- df %>%
  mutate(merged = map2(platemap, tac_data, ~ right_join(.x, .y)))

# Plot standard curve
# Create standards following kit instructions
standards <- tribble(
  ~std, ~conc_mM,
  1,           1,
  2,         0.5,
  3,        0.25,
  4,       0.125,
  5,      0.0625,
  6,     0.03125,
  7,     0.01560,
  8,     0.00780,
  9,     0.00390,
  10,           0
)

std_curve <- df2 %>%
  unnest(merged) %>%
  separate(file, into = c("initials", "protocol", "plate", "type"), remove = FALSE) %>%
  filter(grepl("standard", sample_id)) %>%
  select(type, well, sample_id, wavelength, absorbance) %>%
  rename(std = sample_id) %>% 
  mutate(std = as.numeric(str_sub(std, -1))) %>% 
  group_by(std, type, well) %>%
  summarise(mean_absorbance = mean(absorbance)) %>%
  spread(type, mean_absorbance) %>%
  group_by(std) %>%
  summarise(mean_Buffer = mean(buffer, na.rm = TRUE),
            mean_Cu = mean(Cu, na.rm = TRUE)) %>%
  mutate(net_abs = mean_Cu - mean_Buffer)


std_curve$std <- std_curve$std %>% 
  as.numeric(gsub("0", "10", .)) 

std_curve <- std_curve %>% 
  left_join(standards)

## Visualizing standard curves
std_curve %>%
  ggplot(aes(x = net_abs, y = conc_mM)) +
  geom_point() + theme_bw() +
  labs(title = "Standard curve") 

### check that data looks linear or use non linear 

# Fit linear model
lmod <- lm(conc_mM ~ net_abs, data = std_curve)
# Get fitted values
fitted <- lmod %>% 
  broom::augment(newdata = tibble(net_abs = seq(0, max(std_curve$net_abs), 0.005)))
std_curve_plot <- std_curve %>%
  ggplot(aes(x = net_abs, y = conc_mM)) +
  geom_point(color = "red", size = 3)
std_curve_plot + geom_line(data = fitted, aes(x = net_abs, y = .fitted)) +
  labs(title = "Standard curve")

# Calculate total antioxidant capacity for all samples based on standard curve
tac <- df2 %>%
  unnest(merged) %>%
  separate(file, into = c("initials", "protocol", "plate", "type"), remove = FALSE) %>%
  filter(!grepl("standard", sample_id)) %>%
  filter(!is.na(sample_id)) %>%
  select(type, well, sample_id, wavelength, absorbance) %>% 
  spread(type, absorbance) %>%
  group_by(sample_id) %>%
  summarise(mean_Buffer = mean(buffer), mean_Cu=mean(Cu)) %>%
  # Use standard curve to predict concentrations based on sample absorbance values
  mutate(net_abs = mean_Cu - mean_Buffer,
         uae.mmol.L = map_dbl(net_abs, ~ predict(lmod, newdata = data.frame(net_abs = .))))
std_curve_plot +
  geom_point(data = tac, aes(x = net_abs, y = uae.mmol.L), pch = "X", cex = 3, alpha = 0.3) +
  labs(title = "All samples projected on standard curve")

# Tissue homogenate volume data
metadata <- read_csv("~/Desktop/GITHUB/TL_Astrangia/Raw_data/AP23_ALL_Results.csv") %>% 
  select(sample_id, cage, treatment, full_treatment, Apo_Sym, airbrush_volume, surface_area, prot_ug)
         
# Protein data
# Join homogenate volumes with sample metadata
#metadata <- full_join(metadata, homog_vols) %>%
  #full_join(prot)%>%
  #full_join(sa)%>%
  #full_join(biomass)
# Join TAC data with metadata
tac_results <- left_join(tac, metadata)%>%
  mutate(cre.umol.L = uae.mmol.L * 2189) %>%   # Convert to CRE (see product manual) per unit sample volume
  mutate(cre.umol = cre.umol.L * (airbrush_volume / 1000)) %>% # Convert to CRE per coral by multiplying by homog. vol.
  mutate(cre.umol.mgprot = cre.umol / (prot_ug / 1000)) %>% # Convert to CRE per mg protein by dividing by total protein
  filter(cre.umol >= 0) %>%
  filter(!is.na(Apo_Sym))

#summary table 
summary <- tac_results %>% 
  group_by(full_treatment) %>% 
  summarise(
    count = n(),
    AO_mean = mean(cre.umol),
    AO_sd = sd(cre.umol),
    AO_prot_mean = mean(cre.umol.mgprot),
    AO_prot_sd = sd(cre.umol.mgprot))

# Write data to output file
fin <- tac_results %>% 
  select(sample_id, cre.umol, cre.umol.mgprot) %>%
  write_csv(., file = "~/Desktop/GITHUB/TL_Astrangia/Raw_Data/AP23_Results_Antioxidants_full.csv")

# set comparisons 
treatment_comparisons <- list(c("apo_deep","sym_deep"), c("apo_control","sym_control"), c("apo_shade","sym_shade"))

# Summarize results
ggplot(tac_results) +
  geom_boxplot(aes(full_treatment, cre.umol, group=full_treatment)) +
  labs(x = "", y = "Copper Reducing Equivalents (µmol/mg protein)",
       title = "Total antixodidant capacity") + 
  stat_compare_means(comparisons = treatment_comparisons, method = "wilcox.test", 
                     symnum.args = list(cutpoints = c(0, 0.001 ,0.01, 0.05, Inf), symbols = c("***", "**", "*", "")))

ggplot(tac_results, aes(treatment, cre.umol, color=Apo_Sym)) +
  geom_boxplot()

tac_results %>% 
  group_by(treatment) %>%
  summarise(across(cre.umol, list(mean = mean, min = min, max = max), na.rm = TRUE))

# Plot results by species and site
# Plot all data points with mean ± se
tacplot <- tac_results %>%
  filter(!is.na(cre.umol.mgprot)) %>%
  ggplot(aes(x = full_treatment, y = cre.umol.mgprot, color = Apo_Sym)) +
  #coord_cartesian(ylim = c(0, 1))+
  labs(x = "", y = "Copper Reducing Equivalents (µmol/mg protein)",
       title = "Total antixodidant capacity") +
  geom_jitter(width = 0.1) +                                            # Plot all points
  stat_summary(fun.data = mean_cl_normal, fun.args = list(mult = 1),    # Plot standard error
               geom = "errorbar", color = "black", width = 0.5) +
  stat_summary(fun.y = mean, geom = "point", color = "black")           # Plot mean
tacplot

