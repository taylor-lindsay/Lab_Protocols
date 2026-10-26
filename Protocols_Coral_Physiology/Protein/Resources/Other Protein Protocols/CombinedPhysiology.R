# Written by Hannah Aichelman (hannahaichelman@gmail.com)
# Physiology Analysis
# Each individual physiology parameter is named as a level in the table of contents for ease of navigating the code.
# After you run the code in the 'Read in and format data' and 'Surface Area Measurements' sections, you can jump around to different physiology parameters.

library(ggplot2)
library(lme4)
library(plotly)
library(ggridges)
library(tidyverse)
library(arsenal) #easily compare two data frames
library(Rmisc)
library(lmerTest)
library(emmeans)
library(reshape)
library(readxl)
library(wesanderson)
library(ggpubr)
library(car)
library(sjPlot)
library(effects)
library(glmmTMB)
library(performance)
library(patchwork)
library(magrittr)
library(SciViews)

#for PCAs
library(ggpubr)
library(ggfortify)
library(ggplot2)
library(cluster)
library(FactoMineR) # lots of options for pca visuals and summary stats
library(factoextra)
library(corrplot)
library(dplyr)
library(cowplot)
library(vegan)

##### Read in and format data #####
# read in the data
post_phys <- read.csv('Physiology_Data/data_files/dtvmaster.csv') # physiology data taken at the end of the experiment
init_phys <- read.csv('Physiology_Data/data_files/initial-phys-mod.csv') # physiology data taken at the start of the experiment

# did it load correctly?
head(post_phys)
str(post_phys)

# not reading in calcification or fv/fm data here because we don't have that information for the initial physiology.
# so i am only reading in the data that matches across the initial and final physiology data frames so we can combine them effectively
post_phys2 = post_phys %>%
  dplyr::select(frag, survivedtoend, treat, t3sastan1, t3sastan2, t3sastan3, t3sarec1, t3sarec2, t3sarec3, blastvol, blaster,
                hprotplate, hprot1, hprot2, hprot3, chl630_1, chl630_2, chl630_3, chl663_1, chl663_2, chl663_3, symcount1, symcount2, symcount3) %>%
  dplyr::rename(sastan1 = t3sastan1, sastan2 = t3sastan2, sastan3 = t3sastan3, sarec1 = t3sarec1, sarec2 = t3sarec2, sarec3 = t3sarec3) %>%
  mutate(treat = as.factor(treat), sarec3 = as.numeric(sarec3), blastvol = as.numeric(blastvol)) %>%
  dplyr::filter(frag!="I2A12", frag!="I2A7", frag!="O4D8", frag!="I3E10", frag!="I3I10") #these frags are being removed because they were duplicate genotypes within treatment and had the least complete information of the two

init_phys2 = init_phys %>%
  dplyr::select(frag, survivedtoend, treat, pasastan1, pasastan2, pasastan3, pasarec1, pasarec2, pasarec3, blastvol, blaster,
                hprotplate, hprot1, hprot2, hprot3, chl630_1, chl630_2, chl630_3, chl663_1, chl663_2, chl663_3, symcount1, symcount2, symcount3) %>%
  dplyr::rename(sastan1 = pasastan1, sastan2 = pasastan2, sastan3 = pasastan3, sarec1 = pasarec1, sarec2 = pasarec2, sarec3 = pasarec3) %>%
  mutate(treat = as.factor(treat), survivedtoend = as.character(survivedtoend))

dim(post_phys2)
dim(init_phys2)

combined_phys <- dplyr::bind_rows(post_phys2,init_phys2)
dim(combined_phys)
head(combined_phys)

# add in the carbohydrate data we are confident in that Olivia Nieves worked on
carbs <- read.csv('Physiology_Data/data_files/host_sym_carbs_ON.csv')
dim(carbs)

# merge data frames to keep all of the individuals in the combined physiology spreadsheet and add NA's for any that we don't have carb info for
phys <- left_join(combined_phys, carbs, by = "frag")
dim(phys)
head(phys)

# extract the original two letter site code (inshore/offshore #1-4) from the nubbin ID
phys$origsitecode <- substr(phys$frag, 1, 2)
phys$origsitecode <- as.factor(phys$origsitecode)

# add in site name with inshore/offshore indicator
phys$sitename <- ifelse(phys$origsitecode == 'I2', 'SP',
                        ifelse(phys$origsitecode == 'I3', 'CI',
                               ifelse(phys$origsitecode == 'I4', 'PD',
                                      ifelse(phys$origsitecode == 'O2', 'BS',
                                             ifelse(phys$origsitecode == 'O3', 'CA',
                                                    'BN')))))
phys$sitename <- as.factor(phys$sitename)

# make new nubbin IDs based on the new, more informative site codes
# extract the genotype and frag number
phys$fragid <- substr(phys$frag,3,5)

# add inshore/offshore designation
phys$reef <- substr(phys$frag,1,1)
phys$reef <- ifelse(phys$reef == 'O', 'Outer Reef', 'Inner Reef')
phys$reef <- as.factor(phys$reef)

phys$genet <- substr(phys$fragid,1,1)

#create a new column of combined genotype and site for stats later
phys = phys %>%
  unite(gen_site, c(origsitecode,genet), sep = "", remove = FALSE) %>%
  mutate(gen_site = as.factor(gen_site))

# re-level and re-name treatment
phys$treat <- factor(phys$treat, levels = c("init", "1", "2", "3","4","5"))
levels(phys$treat) <- c("Initial","Control","Low Var","Mod Var","High Var","Control 2")
phys$treat <- as.factor(phys$treat)

str(phys)

## merge physiology data with lineage
# this includes data frame without duplicated preps and with no clones (I4G removed because it was smaller bam file)
lineages = read.csv("Physiology_Data/data_files/tve_lineages_noclones.csv")
head(lineages)

phys <- left_join(phys, lineages, by = "gen_site")
phys$lineage = as.factor(phys$lineage)
head(phys)

## combine physiology data with dominant symbiont type data
head(phys)
# need to separate phys data frame, add its2 types, then re-combine
# first add its2 data to T0 phys dataframe
phys_t0 = phys %>%
  dplyr::filter(treat == "Initial")

its2_types_t0 = read.csv("Physiology_Data/data_files/ITS2.dominanttype.T0.csv")
head(its2_types_t0)

its2_types_t0 = its2_types_t0 %>%
  select(-Sample_or_Control, -treat, -sitename, -reef, -lineage)

phys_its_t0 = left_join(phys_t0, its2_types_t0, by = "gen_site")
head(phys_its_t0)

phys_its_t0 = phys_its_t0 %>%
  dplyr::select(-frag.y) %>%
  dplyr::rename(frag = frag.x)
head(phys_its_t0)

# now do the same for the prestress data frame
phys_ps = phys %>%
  dplyr::filter(treat != "Initial")

its2_types_ps = read.csv("Physiology_Data/data_files/ITS2.dominanttype.prestress.csv")
head(its2_types_ps)

its2_types_ps = its2_types_ps %>%
  select(-Sample_or_Control, -treat, -gen_site, -sitename, -reef, -lineage)

phys_its_ps = left_join(phys_ps, its2_types_ps, by = "frag")
head(phys_its_ps)

dim(phys_its_ps)
dim(phys_its_t0)
colnames(phys_its_ps) == colnames(phys_its_t0)
# same dimensions and column names - ready to recombine

phys = rbind(phys_its_t0, phys_its_ps)
str(phys)
phys$dominant_type = as.factor(phys$dominant_type)
phys$minor_type = as.factor(phys$minor_type)
head(phys)


## make metadata to use throughout the code

phys_metadata = phys %>%
  select(frag, treat, gen_site, sitename, reef, dominant_type, lineage)

head(phys_metadata)
#write.csv(phys_metadata, "Physiology_Data/data_files/phys_metadata.csv",row.names=FALSE)

# set color palettes used throughout
cols_site <- c("CI" = "#543005", "PD"= "#bf812d",  "SP"= "#dfc27d",  "BN" = "#003c30", "BS"= "#35978f", "CA"= "#80cdc1")
#cols_treat <- c("darkgrey", "#FF9966","#CC3300","#7f0000")
cols_treat <- c("darkgrey","#CC3300")
cols_lineage <- c("L1" = "#3f007d", "L2" = "#807dba", "L3" = "#bcbddc")
its2_cols_greens = c("A4" = "#ffeda0", "C1" = "#edf8e9", "C3af" = "#238b45","C3" = "#a1d99b","D1" = "#00441b")

#### Test for Lineage Distribution ####
# Want to test whether distribution of lineages is significantly different across inshore and offshore sites.

str(phys_metadata) # includes all 3 lineages

# first try logit regression: https://stats.oarc.ucla.edu/r/dae/logit-regression/
library(aod)

mylogit = glm(formula = lineage ~ reef, data = phys_metadata, family = "binomial")
summary(mylogit)

# now try chi-square test:
chisq_dat = phys_metadata %>%
  distinct(gen_site, .keep_all = TRUE)

table(chisq_dat$lineage, chisq_dat$reef) # make contingency table

test = chisq.test(table(chisq_dat$lineage, chisq_dat$reef))
test
# Warning message:
#   In chisq.test(table(chisq_dat$lineage, chisq_dat$reef)) :
#   Chi-squared approximation may be incorrect

# Because chisq.test requires an assumption of a minimum of 5 expected counts in the contingency table, we will use a Fisher's exact test.

fisher.test = fisher.test(table(chisq_dat$lineage, chisq_dat$reef))
fisher.test

# Fisher's Exact Test for Count Data
#
# data:  table(chisq_dat$lineage, chisq_dat$reef)
# p-value = 8.707e-08
# alternative hypothesis: two.sided

# Or can remove low count L3 from dataset and still use Chi-Square test - THIS IS REPORTED IN MANUSCRIPT
chisq_dat_noL3 = phys_metadata %>%
  distinct(gen_site, .keep_all = TRUE) %>%
  dplyr::filter(lineage != "L3") %>%
  droplevels()

test = chisq.test(table(chisq_dat_noL3$lineage, chisq_dat_noL3$reef))
test

# Pearson's Chi-squared test with Yates' continuity correction
#
# data:  table(chisq_dat_noL3$lineage, chisq_dat_noL3$reef)
# X-squared = 23.577, df = 1, p-value = 1.2e-06

#### Surface Area Measurements ####
# get surface area measurements in cm^2
# pa (pre-acclimation) measurements for surface area are suitable only for initial phys nubbins. physiology for the experimental nubbins must use T3 surface area
# I used the appropriate surface area measurement for each time point, but combined into one column for simplicity
# look above for more details - read in initial phys and post-variability phys separately, each data frame with its own SA measurements included

phys$SAcm2 <- (4*((phys$sarec1+phys$sarec2+phys$sarec3)/3))/((phys$sastan1+phys$sastan2+phys$sastan3)/3)

# check to see if surface areas look normal
# histogram of nubbin sizes
ggplot(phys, aes(x=SAcm2, color = sitename, fill = sitename))+
  theme_bw()+
  geom_histogram(binwidth = 1)+
  scale_fill_manual(values = cols_site)+
  scale_color_manual(values = cols_site)+
  labs(x=expression(paste("Surface Area (cm"^2*')')))+
  facet_wrap(~sitename, 6,1)

# boxplot of nubbin sizes
plot <- ggplot(phys, aes(x=sitename, y=SAcm2, fill = sitename, text = paste("Nubbin:",frag)))+
  theme_bw()+
  geom_boxplot()+
  geom_jitter(aes(color = sitename),width = 0.3)+
  scale_fill_manual(values = cols_site)+
  scale_color_manual(values = cols_site)+
  theme(axis.title.x = element_blank())+
  theme(axis.text.x=element_text(angle = 45, vjust = 0.5))
#labs(y=expression(paste("Surface Area (cm"^2*')')))
ggplotly(subplot(list(plot),nrows=1,titleY=F) %>% layout(showlegend=T))

##### Tissue Thickness (T0) #####

# tissue thickness is measured in mm and was only taken from nubbins from the start of the experiment (T0)

# read in initial physiology data file again
init_phys_tiss <- read.csv('Physiology_Data/data_files/initial-phys-mod.csv')

tiss_phys = init_phys_tiss %>%
  select(frag,treat,tisthick1,tisthick2,tisthick3,tisthick4,tisthick5,tisthick6,
         origsitecode,sitename,fragid,nubbin,reef,SAcm2) %>%
  mutate(treat = as.factor(treat), reef = as.factor(reef), sitename = as.factor(sitename))

# re-name treatment
levels(tiss_phys$treat) <- c("Initial")

#re-name site names
tiss_phys$sitename <- factor(tiss_phys$sitename, levels = c("I-Cristobal","I-Punta Donato","I-STRI Point",
                                                            "O-Bastimentos N","O-Bastimentos S","O-Cayo de Agua"))
levels(tiss_phys$sitename) <- c("CI","PD","SP","BN","BS","CA")

# add identifying data
# make new nubbin IDs based on the new, more informative site codes
tiss_phys$genet <- substr(tiss_phys$fragid,1,1)

#create a new column of combined genotype and site for stats later
# calculate average tissue thickness for each nubbin
tiss_phys = tiss_phys %>%
  unite(gen_site, c(origsitecode,genet), sep = "", remove = FALSE) %>%
  mutate(gen_site = as.factor(gen_site)) %>%
  mutate(avgtiss=rowMeans(.[ , c("tisthick1","tisthick2","tisthick3","tisthick4","tisthick5","tisthick6")], na.rm=TRUE)) %>%
  dplyr::filter(frag!= "I3E6_v2") %>% # redundant info, same as other I3E sample
  dplyr::filter(gen_site != "I4G") # clone pair, removed this sample from 2brad so removing from phys too


# merge with lineage info for later plotting
# this includes data frame without duplicated preps, and with one of the clone pair (I4G) removed
lineages = read.csv("Physiology_Data/data_files/tve_lineages_noclones.csv")

tiss_phys_all_lin <- left_join(tiss_phys, lineages, by = "gen_site")

str(tiss_phys_all_lin)
tiss_phys_all_lin$gen_site = as.factor(tiss_phys_all_lin$gen_site)
tiss_phys_all_lin$lineage = as.factor(tiss_phys_all_lin$lineage)

tiss_phys_2_lin = tiss_phys_all_lin %>%
  dplyr::filter(is.na(lineage) | lineage!="L3") # want to keep NA values for lineage here since they still have other info, will remove na's for lineage specific plots

# can't combine with dominant sym type here because tissue thickness is from initial phys and sym types were from end of variability

# STATS
m1 <- lm(avgtiss ~ lineage, data = tiss_phys_all_lin)
summary(m1)
# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)
# (Intercept)   5.1923     0.1526  34.031  < 2e-16 ***
# lineageL2    -0.6256     0.2583  -2.422 0.019739 *
# lineageL3    -1.7534     0.4905  -3.575 0.000881 ***
#
Anova(m1)
# Anova Table (Type II tests)
#
# Response: avgtiss
# Sum Sq Df F value   Pr(>F)
# lineage   10.432  2  8.0019 0.001111 **
# Residuals 28.028 43

lsmeans(m1, pairwise~lineage, adjust="tukey")
# $contrasts
# contrast estimate    SE df t.ratio p.value
# L1 - L2     0.626 0.258 43   2.422  0.0507
# L1 - L3     1.753 0.490 43   3.575  0.0025
# L2 - L3     1.128 0.511 43   2.209  0.0811

# PLOTS

#SummarySE to format data for plotting
tiss_means_site_all_lin <- summarySE(tiss_phys_all_lin, measurevar="avgtiss", groupvars=c("treat","sitename"))
tiss_means_site_2_lin <- summarySE(tiss_phys_2_lin, measurevar="avgtiss", groupvars=c("treat","sitename"))

# plot, treatment x axis colored by site data figure
tiss_plot_site <- ggplot(tiss_means_site_2_lin,aes(x = sitename, y = avgtiss, fill = sitename, color = sitename))+
  theme_bw()+
  geom_errorbar(aes(x = sitename, ymax = avgtiss+se, ymin = avgtiss-se, color = sitename), width = .2, position = position_dodge(width=0.3)) +
  geom_point(size = 3, pch = 21, position = position_dodge(width=0.3), color = "black")+
  scale_fill_manual(name = "Site",
                     values = cols_site)+
  scale_color_manual(name = "Site",
                     values = cols_site)+
  xlab("Site Name")+
  ylab("Tissue Thickness (mm)")+
  ylim(3.75,6) +
  #geom_vline(xintercept = 1.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  theme(legend.position = "none")
tiss_plot_site

ggsave(tiss_plot_site, filename = "/Users/hannahaichelman/Documents/BU/TVE/TissueThickness/tissuethickness_site_2_lin.pdf", width=4, height=4, units=c("in"), useDingbats=FALSE)


#SummarySE to format data for plotting with lineage
tiss_phys_all_lin_nona = tiss_phys_all_lin %>%
  drop_na(lineage)

tiss_phys_2_lin_nona = tiss_phys_2_lin %>%
  drop_na(lineage)

tiss_means_all_lin <- summarySE(tiss_phys_all_lin_nona, measurevar="avgtiss", groupvars=c("lineage"))
tiss_means_2_lin <- summarySE(tiss_phys_2_lin_nona, measurevar="avgtiss", groupvars=c("lineage"))

# plot, lineage on x axis
tiss_plot_lineage <- ggplot(tiss_means_all_lin,aes(x = lineage, y = avgtiss, fill = lineage, color = lineage))+
  theme_bw()+
  geom_errorbar(aes(x = lineage, ymax = avgtiss+se, ymin = avgtiss-se, color = lineage), width = .2, position = position_dodge(width=0.3)) +
  geom_point(size = 3, position = position_dodge(width=0.3), shape = 21, color = "black")+
  scale_fill_manual(name = "Lineage",
                    values = cols_lineage)+
  scale_color_manual(name = "Lineage",
                     values = cols_lineage)+
  xlab("Lineage")+
  ylab("Tissue Thickness (mm)")+
  #ylim(3.75,6) +
  #geom_vline(xintercept = 1.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  theme(legend.position = "none")
tiss_plot_lineage

ggsave(tiss_plot_lineage, filename = "/Users/hannahaichelman/Documents/BU/TVE/TissueThickness/tissuethickness_lineage_alllin.pdf", width=4, height=4, units=c("in"), useDingbats=FALSE)

# plot just CI corals to see if lineage difference holds
tiss_phys_lineage_CI = tiss_phys_2_lin_nona %>%
  subset(sitename == "CI")

tiss_means_lineage_CI <- summarySE(tiss_phys_lineage_CI, measurevar="avgtiss", groupvars=c("lineage"))

m1 <- lm(avgtiss ~ lineage, data = tiss_phys_lineage_CI)
summary(m1)
anova(m1)
# Response: avgtiss
# Df Sum Sq Mean Sq F value Pr(>F)
# lineage    1 0.0764 0.07640  0.2372 0.6469
# Residuals  5 1.6107 0.32215

lsmeans(m1, pairwise~lineage, adjust="tukey")

# plot, lineage on x axis
tiss_plot_lineage_CI <- ggplot(tiss_means_lineage_CI,aes(x = lineage, y = avgtiss, color = lineage))+
  theme_bw()+
  geom_point(size = 3, position = position_dodge(width=0.3))+
  geom_errorbar(aes(x = lineage, ymax = avgtiss+se, ymin = avgtiss-se), width = .2, position = position_dodge(width=0.3)) +
  scale_color_manual(name = "Lineage",
                     labels = c("L1","L2"),
                     values = cols_lineage)+
  xlab("Lineage")+
  ylab("Tissue Thickness (mm)")+
  #ylim(3.75,6) +
  #geom_vline(xintercept = 1.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  theme(legend.position = "none")
tiss_plot_lineage_CI

ggsave(tiss_plot_lineage_CI, filename = "/Users/hannahaichelman/Documents/BU/TVE/TissueThickness/tissuethickness_lineage_CIonly.pdf", width=4, height=4, units=c("in"), useDingbats=FALSE)


# plot just SP corals to see if lineage difference holds
tiss_phys_lineage_SP = tiss_phys_all_lin_nona %>%
  subset(sitename == "SP")

tiss_means_lineage_SP <- summarySE(tiss_phys_lineage_SP, measurevar="avgtiss", groupvars=c("lineage"))

# plot, reef zone x axis
tiss_plot_lineage_SP <- ggplot(tiss_means_lineage_SP,aes(x = lineage, y = avgtiss, color = lineage))+
  theme_bw()+
  geom_point(size = 3, position = position_dodge(width=0.3))+
  geom_errorbar(aes(x = lineage, ymax = avgtiss+se, ymin = avgtiss-se), width = .2, position = position_dodge(width=0.3)) +
  scale_color_manual(name = "Lineage",
                     labels = c("L2","L3"),
                     values = cols_lineage)+
  xlab("Lineage")+
  ylab("Tissue Thickness (mm)")+
  #ylim(3.75,6) +
  #geom_vline(xintercept = 1.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  theme(legend.position = "none")
tiss_plot_lineage_SP

ggsave(tiss_plot_lineage_SP, filename = "/Users/hannahaichelman/Documents/BU/TVE/TissueThickness/tissuethickness_lineage_SPonly.pdf", width=4, height=4, units=c("in"), useDingbats=FALSE)

##### Protein Concentrations #####
# first, subset physiology data add in standard curve data needed to correct for separate protein plates
hprot_phys = phys %>%
  select(frag,treat,hprotplate,hprot1,hprot2,hprot3,
         survivedtoend,blastvol,gen_site,origsitecode,sitename,fragid,
         reef,genet,SAcm2, dominant_type, lineage)

head(hprot_phys)

# we need to adjust the mean absorbance values by the intercept value of each plate
hprot_phys$protx2coef <- as.numeric(ifelse(hprot_phys$hprotplate == 'D1', 0.472799,
                                           ifelse(hprot_phys$hprotplate == 'D2', 0.385059,
                                                  ifelse(hprot_phys$hprotplate == 'D3', 0.618482,
                                                         ifelse(hprot_phys$hprotplate == 'D4', 0.539057,
                                                                ifelse(hprot_phys$hprotplate == 'D5', 0.637505,
                                                                       ifelse(hprot_phys$hprotplate == 'D6', 0.4492627,
                                                                              ifelse(hprot_phys$hprotplate == 'D7', 0.456354,
                                                                                     ifelse(hprot_phys$hprotplate == 'D8', 0.526895,
                                                                                            ifelse(hprot_phys$hprotplate == 'O1', 0.0650601,
                                                                                                   ifelse(hprot_phys$hprotplate == 'O2', 0.753409,
                                                                                                          ifelse(hprot_phys$hprotplate == 'O3', 0.4325502,
                                                                                                                 ifelse(hprot_phys$hprotplate == 'O4', 0.554534,
                                                                                                                        ifelse(hprot_phys$hprotplate == 'O5', 0.680815,
                                                                                                                               ifelse(hprot_phys$hprotplate == 'O6', 0.772513,
                                                                                                                                      'NA')))))))))))))))

hprot_phys$protxcoef <- as.numeric(ifelse(hprot_phys$hprotplate == 'D1', -0.061659,
                                          ifelse(hprot_phys$hprotplate == 'D2', 0.038351,
                                                 ifelse(hprot_phys$hprotplate == 'D3', -0.073079,
                                                        ifelse(hprot_phys$hprotplate == 'D4', -0.048588,
                                                               ifelse(hprot_phys$hprotplate == 'D5', -0.029648,
                                                                      ifelse(hprot_phys$hprotplate == 'D6', 0.0452952,
                                                                             ifelse(hprot_phys$hprotplate == 'D7', -0.014998,
                                                                                    ifelse(hprot_phys$hprotplate == 'D8', -0.031275,
                                                                                           ifelse(hprot_phys$hprotplate == 'O1', 0.0777724,
                                                                                                  ifelse(hprot_phys$hprotplate == 'O2', -0.021612,
                                                                                                         ifelse(hprot_phys$hprotplate == 'O3', 0.0287098,
                                                                                                                ifelse(hprot_phys$hprotplate == 'O4', -0.029577,
                                                                                                                       ifelse(hprot_phys$hprotplate == 'O5', -0.027620,
                                                                                                                              ifelse(hprot_phys$hprotplate == 'O6', -0.048054,
                                                                                                                                     'NA')))))))))))))))

hprot_phys$protint <- as.numeric(ifelse(hprot_phys$hprotplate == 'D1', 0.004528,
                                        ifelse(hprot_phys$hprotplate == 'D2',0.002527,
                                               ifelse(hprot_phys$hprotplate == 'D3', 0.004711,
                                                      ifelse(hprot_phys$hprotplate == 'D4', 0.003552,
                                                             ifelse(hprot_phys$hprotplate == 'D5', 0.003320,
                                                                    ifelse(hprot_phys$hprotplate == 'D6', 0.0018613,
                                                                           ifelse(hprot_phys$hprotplate == 'D7', 0.002955,
                                                                                  ifelse(hprot_phys$hprotplate == 'D8', 0.003667,
                                                                                         ifelse(hprot_phys$hprotplate == 'O1', 0.0005224,
                                                                                                ifelse(hprot_phys$hprotplate == 'O2', 0.002474,
                                                                                                       ifelse(hprot_phys$hprotplate == 'O3', 0.0007926,
                                                                                                              ifelse(hprot_phys$hprotplate == 'O4', 0.002101,
                                                                                                                     ifelse(hprot_phys$hprotplate == 'O5', 0.003050,
                                                                                                                            ifelse(hprot_phys$hprotplate == 'O6', 0.002921,
                                                                                                                                   'NA')))))))))))))))

head(hprot_phys)

# now we can use the quadratic equation to solve for absorbance
# y = ax^2 + bx + c where
# y = protein concentration
# x = raw absorbance from Bradford Assay

# these equations get the protein values in ug/ul
hprot_phys$calcprot1 <- ((hprot_phys$protx2coef*(hprot_phys$hprot1^2)) + (hprot_phys$protxcoef*hprot_phys$hprot1) + hprot_phys$protint)
hprot_phys$calcprot2 <- ((hprot_phys$protx2coef*(hprot_phys$hprot2^2)) + (hprot_phys$protxcoef*hprot_phys$hprot2) + hprot_phys$protint)
hprot_phys$calcprot3 <- ((hprot_phys$protx2coef*(hprot_phys$hprot3^2)) + (hprot_phys$protxcoef*hprot_phys$hprot3) + hprot_phys$protint)

# calculate average protein concentration
hprot_phys = hprot_phys %>%
  mutate(avgprot=rowMeans(.[ , c("calcprot1","calcprot2","calcprot3")], na.rm=TRUE))

# [prot concentration in ug/ul] = ([sample protein value] * [1000 ul in 1 ml / 80 ul of sample volume (where sample volume includes dilution with seawater)] * number of mls in slurry / [surface area of fragment])
# 80/6.4 is the dilution factor - 80 uL total volume with 6.4 uL of sample
hprot_phys$prot_ugcm2 <- (hprot_phys$avgprot*(80/6.4)*1000*hprot_phys$blastvol)/hprot_phys$SAcm2
hprot_phys$prot_mgcm2 <- hprot_phys$prot_ugcm2/1000

# How many protein NA's?
sum(is.na(hprot_phys$prot_mgcm2))
# 110

# But how many of these are missing and how many are just nubbins that died?
missing_prot_inds = hprot_phys %>%
  filter(survivedtoend=="yes") %>%
  filter(is.na(hprot1))
# 5 have no protein record from the absorbance stage: I2C3, I3A7, O2C8, I3B4, I4E7, but all of these have no blaster or blastvol info

missing_blastinfo = hprot_phys %>%
  filter(survivedtoend=="yes") %>%
  filter(is.na(blastvol))
# "O3F2" "O3H2" "I2C3" "I3A7" "O2C8" "I4G3" "O2G7" "I3B4" "I3D8" "I3F4" "I4E7" "O4F4" "I4G4" "I3E5" "I3F5" "I3G5"

# make data subsets for stats and plotting - removing low and high variability treatments this time through.
hprot_phys_all_lin = hprot_phys %>%
  drop_na(prot_mgcm2) %>%
  dplyr::filter(treat!="Control 2") %>%
  dplyr::filter(treat!="Low Var") %>%
  dplyr::filter(treat!="High Var") %>%
  dplyr::filter(frag != "I3E6_v2") %>% # no consensus on this I3E6 individual phenotype for this measure - so removing both
  dplyr::filter(frag != "I3E6") %>%
  dplyr::filter(gen_site != "I4G") %>% # clone with I4F, remove from dataset
  select(frag, gen_site, treat, sitename, dominant_type, lineage, prot_mgcm2)

hprot_phys_2_lin = hprot_phys_all_lin %>%
  dplyr::filter(is.na(lineage) | lineage!="L3") %>%
  drop_na(sitename, lineage, treat, dominant_type) # doing this for plotting

initial_hprot_phys_all_lin = hprot_phys_all_lin %>%
  dplyr::filter(treat=="Initial")

initial_hprot_phys_2_lin = initial_hprot_phys_all_lin %>%
  dplyr::filter(is.na(lineage) | lineage!="L3") # want to keep NA values for lineage here since they still have other info, will remove na's for lineage specific plots

post_hprot_phys_all_lin = hprot_phys_all_lin %>%
  dplyr::filter(treat!="Initial")

post_hprot_phys_2_lin = post_hprot_phys_all_lin %>%
  dplyr::filter(is.na(lineage) | lineage!="L3")

# the values we have here seem to be similar to protein values from Wright et al. GCB, good sanity check.

# STATS
## NEED TO RUN STATS SEPARATELY FOR INITIAL AND FINAL TIME POINTS

model.prot.initial = lm(prot_mgcm2 ~ lineage, data = initial_hprot_phys_2_lin)
summary(model.prot.initial)
anova(model.prot.initial)
# Response: prot_mgcm2
#         Df   Sum Sq  Mean Sq F value   Pr(>F)
# lineage    1 0.078953 0.078953  12.002 0.001282 **
# Residuals 40 0.263142 0.006579

check_model(model.prot.initial)


model.prot.final <- lmer(prot_mgcm2 ~ treat+lineage + (1|gen_site), data = post_hprot_phys_2_lin)
summary(model.prot.final)
anova(model.prot.final)
# Type III Analysis of Variance Table with Satterthwaite's method
#           Sum Sq  Mean Sq NumDF  DenDF F value   Pr(>F)
# treat   0.003558 0.003558     1 21.929  0.3924 0.537489
# lineage 0.107332 0.107332     1 26.358 11.8368 0.001947 **

check_model(model.prot.final)


#specify model (because we are interested in pairwise, have to include the interaction)
m.emm<- lmer(prot_mgcm2 ~ treat*lineage + (1|gen_site), data = post_hprot_phys_2_lin, REML=FALSE)

emms<-emmeans(m.emm, ~lineage|treat) #, adjust="Bonferoni"
pairs(emms, interaction = "pairwise") %>% rbind(adjust="fdr")
# treat   lineage_pairwise estimate     SE   df t.ratio p.value
# Control L1 - L2             0.128 0.0420 79.9   3.045  0.0063
# Mod Var L1 - L2             0.102 0.0394 77.1   2.579  0.0118


# PLOTS
#SummarySE to format data for plotting
prot_means_site_2_lin <- summarySE(hprot_phys_2_lin, measurevar="prot_mgcm2", groupvars=c("treat","sitename"))

# plot, treatment x axis colored by site data figure
prot_plot_site <- ggplot(prot_means_site_2_lin,aes(x = treat, y = prot_mgcm2, fill = sitename))+
  theme_bw()+
  geom_errorbar(aes(x = treat, ymax = prot_mgcm2+se, ymin = prot_mgcm2-se), width = .2, position = position_dodge(width=0.3)) +
  geom_point(size = 3, pch = 21, position = position_dodge(width=0.3))+
  scale_fill_manual(name = "Site",
                     values = cols_site)+
  xlab("Treatment")+
  ylab(bquote("Total Protein (mg" ~cm^-2~')'))+
  ylim(0,0.4) +
  geom_vline(xintercept = 1.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        legend.position = c(.95, .99),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(5, 5, 5, 5)
        #legend.key = element_rect(fill = "none")
)
prot_plot_site

ggsave(prot_plot_site, filename = "/Users/hannahaichelman/Dropbox/BU/TVE/Protein/plots/protein_site_2_lin.pdf", width=5, height=4, units=c("in"), useDingbats=FALSE)


# plot, treatment x axis and lineage color
#SummarySE to format data for plotting
prot_means_lineage_2_lin <- summarySE(hprot_phys_2_lin, measurevar="prot_mgcm2", groupvars=c("treat","lineage"))

# plot, treatment x axis colored by site data figure
prot_plot_lineage <- ggplot(prot_means_lineage_2_lin,aes(x = treat, y = prot_mgcm2, fill = lineage))+
  theme_bw()+
  geom_errorbar(aes(x = treat, ymax = prot_mgcm2+se, ymin = prot_mgcm2-se), width = .2, position = position_dodge(width=0.3)) +
  geom_point(size = 3, pch = 21, position = position_dodge(width=0.3))+
  scale_fill_manual(name = "Lineage",
                     values = cols_lineage)+
  xlab("Treatment")+
  ylab(bquote("Total Protein (mg" ~cm^-2~')'))+
  ylim(0,0.4) +
  geom_vline(xintercept = 1.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        legend.position = c(.95, .99),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(5, 5, 5, 5)
        #legend.key = element_rect(fill = "none")
)
prot_plot_lineage

ggsave(prot_plot_lineage, filename = "/Users/hannahaichelman/Dropbox/BU/TVE/Protein/plots/protein_lineage_2_lin.pdf", width=5, height=4, units=c("in"), useDingbats=FALSE)


# plot, treatment x axis and symtype color
#SummarySE to format data for plotting
prot_means_sym <- summarySE(hprot_phys_2_lin, measurevar="prot_mgcm2", groupvars=c("treat","dominant_type"))

# plot, treatment x axis colored by site data figure
prot_plot_sym <- ggplot(prot_means_sym,aes(x = treat, y = prot_mgcm2, fill = dominant_type))+
  theme_bw()+
  geom_errorbar(aes(x = treat, ymax = prot_mgcm2+se, ymin = prot_mgcm2-se), width = .2, position = position_dodge(width=0.3)) +
  geom_point(size = 3, pch = 21, position = position_dodge(width=0.3))+
  scale_fill_manual(name = "Dominant Symbiont",
                     values = its2_cols_greens)+
  xlab("Treatment")+
  ylab(bquote("Total Protein (mg" ~cm^-2~')'))+
  ylim(0,0.4) +
  geom_vline(xintercept = 1.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        legend.position = c(.95, .99),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(5, 5, 5, 5)
        #legend.key = element_rect(fill = "none")
)
prot_plot_sym

ggsave(prot_plot_sym, filename = "/Users/hannahaichelman/Dropbox/BU/TVE/Protein/plots/protein_lineage_sym.pdf", width=6, height=4, units=c("in"), useDingbats=FALSE)



# Combine these protein plots
prot_plots_all = ggarrange(prot_plot_site, prot_plot_lineage, prot_plot_sym,
                           ncol = 3, nrow = 1)
ggsave(prot_plots_all, filename = "/Users/hannahaichelman/Dropbox/BU/TVE/Protein/plots/protein_all_plots.pdf", width=12, height=4, units=c("in"), useDingbats=FALSE)



# plot just CI corals to see if lineage difference holds
hprot_phys_2_lin_CI = hprot_phys_2_lin %>%
  subset(sitename == "CI")

hprot_means_lineage_CI <- summarySE(hprot_phys_2_lin_CI, measurevar="prot_mgcm2", groupvars=c("treat","lineage"))

prot_plot_lineage_CI <- ggplot(hprot_means_lineage_CI,aes(x = treat, y = prot_mgcm2, color = lineage, pch = lineage))+
  theme_bw()+
  geom_point(size = 3, position = position_dodge(width=0.3))+
  geom_errorbar(aes(x = treat, ymax = prot_mgcm2+se, ymin = prot_mgcm2-se), width = .2, position = position_dodge(width=0.3)) +
  scale_color_manual(name = "Lineage",
                     breaks = c("L1","L2"),
                     values = cols_lineage)+
  scale_shape_manual(name = "Lineage",
                     breaks = c("L1","L2"),
                     values=c(19,17))+
  xlab("Treatment")+
  ylab(bquote("Total Protein (mg" ~cm^-2~')'))+
  #ylim(0,0.5) +
  geom_vline(xintercept = 1.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))
prot_plot_lineage_CI

ggsave(prot_plot_lineage_CI, filename = "/Users/hannahaichelman/Dropbox/BU/TVE/Protein/plots/protein_lineage_CIonly.pdf", width=5, height=4, units=c("in"), useDingbats=FALSE)

##### Host Carbohydrate Concentrations #####
# first, subset physiology data add in standard curve data needed to correct for separate carbohydrate plates
hcarb_phys = phys %>%
  select(frag,treat,hcarbplate_ha,hcarb1_ha,hcarb2_ha,hcarb3_ha,
         survivedtoend,blastvol,gen_site,origsitecode,sitename,fragid,reef,genet,SAcm2,dominant_type,lineage)

head(hcarb_phys)

# used this to sort out who we are still missing carbs from
missing_hcarbs = hcarb_phys %>%
  dplyr::filter(survivedtoend=="yes") %>%
  dplyr::filter(is.na(hcarb1_ha))
# O2C8, I3B4, I4E7

# For both host carbs and symbiont carbs, for all time points,
# we are missing 3, but both of these do not have blaster or blastvol info so not much we can do

# Add in standard curve info for each plate
hcarb_phys$hcarbxcoef <- as.numeric(ifelse(hcarb_phys$hcarbplate_ha == 'O1', 0.1678088,
                                           ifelse(hcarb_phys$hcarbplate_ha == 'O2', 0.187208,
                                                  ifelse(hcarb_phys$hcarbplate_ha == 'O3', 0.229014,
                                                         ifelse(hcarb_phys$hcarbplate_ha == 'O4', 0.231081,
                                                                ifelse(hcarb_phys$hcarbplate_ha == 'O5', 0.196008,
                                                                       ifelse(hcarb_phys$hcarbplate_ha == 'O6', 0.162171,
                                                                              ifelse(hcarb_phys$hcarbplate_ha == 'O7', 0.1765377,
                                                                                     ifelse(hcarb_phys$hcarbplate_ha == 'O8', 0.141797,
                                                                                            ifelse(hcarb_phys$hcarbplate_ha == 'O9', 0.152009,
                                                                                                   ifelse(hcarb_phys$hcarbplate_ha == 'O10', 0.1547928,
                                                                                                          ifelse(hcarb_phys$hcarbplate_ha == 'O11', 0.170403,
                                                                                                                 ifelse(hcarb_phys$hcarbplate_ha == 'O12', 0.142706,
                                                                                                                        ifelse(hcarb_phys$hcarbplate_ha == 'O13', 0.148280,
                                                                                                                               ifelse(hcarb_phys$hcarbplate_ha == 'O14', 0.164326,
                                                                                                                                      ifelse(hcarb_phys$hcarbplate_ha == 'O15', 0.146469,
                                                                                                                                             ifelse(hcarb_phys$hcarbplate_ha == 'O16', 0.175768,
                                                                                                                                                    ifelse(hcarb_phys$hcarbplate_ha == 'O17', 0.163126,
                                                                                                                                                           ifelse(hcarb_phys$hcarbplate_ha == 'O35', 0.1536346,
                                                                                                                                                                  ifelse(hcarb_phys$hcarbplate_ha == 'O36', 0.1717878,
                                                                                                                                                                         'NA'))))))))))))))))))))

hcarb_phys$hcarbint <- as.numeric(ifelse(hcarb_phys$hcarbplate_ha == 'O1', 0.0005952,
                                         ifelse(hcarb_phys$hcarbplate_ha == 'O2', -0.011863,
                                                ifelse(hcarb_phys$hcarbplate_ha == 'O3', -0.005677,
                                                       ifelse(hcarb_phys$hcarbplate_ha == 'O4', -0.001425,
                                                              ifelse(hcarb_phys$hcarbplate_ha == 'O5', -0.001076,
                                                                     ifelse(hcarb_phys$hcarbplate_ha == 'O6', -0.007433,
                                                                            ifelse(hcarb_phys$hcarbplate_ha == 'O7', 0.0005785,
                                                                                   ifelse(hcarb_phys$hcarbplate_ha == 'O8', -0.003569,
                                                                                          ifelse(hcarb_phys$hcarbplate_ha == 'O9', -0.001455,
                                                                                                 ifelse(hcarb_phys$hcarbplate_ha == 'O10', -0.0003914,
                                                                                                        ifelse(hcarb_phys$hcarbplate_ha == 'O11', -0.001060,
                                                                                                               ifelse(hcarb_phys$hcarbplate_ha == 'O12', 0.002935,
                                                                                                                      ifelse(hcarb_phys$hcarbplate_ha == 'O13', 0.002679,
                                                                                                                             ifelse(hcarb_phys$hcarbplate_ha == 'O14', -0.004464,
                                                                                                                                    ifelse(hcarb_phys$hcarbplate_ha == 'O15', 0.008575,
                                                                                                                                           ifelse(hcarb_phys$hcarbplate_ha == 'O16', 0.002615,
                                                                                                                                                  ifelse(hcarb_phys$hcarbplate_ha == 'O17', -0.006809,
                                                                                                                                                         ifelse(hcarb_phys$hcarbplate_ha == 'O35', -0.0004174,
                                                                                                                                                                ifelse(hcarb_phys$hcarbplate_ha == 'O36', 0.0003559,
                                                                                                                                                                       'NA'))))))))))))))))))))

head(hcarb_phys)

# get the carbohydrate values in mg/mL (concentration of D-glucose standard)
# the standard curve equation for carbs is linear:  y = m*x + b
# where y is the calculated concentration in mg/mL, m is the slope, or x coefficient, and b is the intercept
hcarb_phys$calchcarb1 <- ((hcarb_phys$hcarbxcoef*hcarb_phys$hcarb1_ha) + hcarb_phys$hcarbint)
hcarb_phys$calchcarb2 <- ((hcarb_phys$hcarbxcoef*hcarb_phys$hcarb2_ha) + hcarb_phys$hcarbint)
hcarb_phys$calchcarb3 <- ((hcarb_phys$hcarbxcoef*hcarb_phys$hcarb3_ha) + hcarb_phys$hcarbint)

# calculate average carbohydrate concentration, mg/mL
# this code allows us to take into account when there was an outlier replaced with an NA
hcarb_phys = hcarb_phys %>%
  mutate(avghcarb=rowMeans(.[ , c("calchcarb1","calchcarb2","calchcarb3")], na.rm=TRUE))

# [carb concentration in mg/cm2] = ([mean concentration mg/mL] * [dilution factor] * number of mls in slurry / [surface area of fragment]
# dilution factor is 5 here because in the protocol, each reaction has 10 uL of sample and 40 uL of seawater
hcarb_phys$hcarb_mgcm2 <- (hcarb_phys$avghcarb*5*hcarb_phys$blastvol)/hcarb_phys$SAcm2

# make data subsets for stats and plotting
hcarb_phys_all_lin = hcarb_phys %>%
  drop_na(hcarb_mgcm2) %>%
  dplyr::filter(treat!="Control 2") %>%
  dplyr::filter(treat!="Low Var") %>%
  dplyr::filter(treat!="High Var") %>%
  dplyr::filter(gen_site!="I4G") %>% # remove clone
  select(frag, gen_site, treat, sitename, dominant_type, lineage, hcarb_mgcm2)

hcarb_phys_2_lin = hcarb_phys_all_lin %>%
  dplyr::filter(is.na(lineage) | lineage!="L3") %>%
  drop_na(sitename, lineage, treat, dominant_type) # doing this for plotting

initial_hcarb_phys_all_lin = hcarb_phys_all_lin %>%
  dplyr::filter(treat=="Initial")

initial_hcarb_phys_2_lin = initial_hcarb_phys_all_lin %>%
  dplyr::filter(is.na(lineage) | lineage!="L3") # want to keep NA values for lineage here since they still have other info, will remove na's for lineage specific plots

post_hcarb_phys_all_lin = hcarb_phys_all_lin %>%
  dplyr::filter(treat!="Initial")

post_hcarb_phys_2_lin = post_hcarb_phys_all_lin %>%
  dplyr::filter(is.na(lineage) | lineage!="L3")


# STATS
## NEED TO RUN STATS SEPARATELY FOR INITIAL AND FINAL TIME POINTS

model.hcarb.initial = lm(hcarb_mgcm2 ~ lineage, data = initial_hcarb_phys_2_lin)
summary(model.hcarb.initial)
anova(model.hcarb.initial)
# Response: hcarb_mgcm2
#         Df   Sum Sq  Mean Sq F value  Pr(>F)
# lineage    1 0.052428 0.052428  7.1386 0.01077 *
# Residuals 41 0.301115 0.007344

check_model(model.hcarb.initial)


model.hcarb.final <- lmer(hcarb_mgcm2 ~ lineage+treat + (1|gen_site), data = post_hcarb_phys_2_lin)
summary(model.hcarb.final)
anova(model.hcarb.final)
# Type III Analysis of Variance Table with Satterthwaite's method
#          Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)
# lineage 0.29784 0.29784     1 45.068  16.462 0.0001947 ***
# treat   0.37173 0.37173     1 44.136  20.546 4.412e-05 ***

check_model(model.hcarb.final)


#specify model (because we are interested in pairwise, have to include the interaction)
m.emm<- lmer(hcarb_mgcm2 ~ treat*lineage + (1|gen_site), data = post_hcarb_phys_2_lin, REML=FALSE)

emms<-emmeans(m.emm, ~lineage|treat) #, adjust="Bonferoni"
pairs(emms, interaction = "pairwise") %>% rbind(adjust="fdr")
# treat   lineage_pairwise estimate     SE   df t.ratio p.value
# Control L1 - L2            0.2076 0.0473 83.0   4.393  0.0001
# Mod Var L1 - L2            0.0752 0.0441 82.8   1.706  0.0917


## PLOTS
#SummarySE to format data for plotting
hcarb_means_site_2_lin <- summarySE(hcarb_phys_2_lin, measurevar="hcarb_mgcm2", groupvars=c("treat","sitename"))

# plot, treatment x axis colored by site data figure
hcarb_plot_site <- ggplot(hcarb_means_site_2_lin,aes(x = treat, y = hcarb_mgcm2, fill = sitename))+
  theme_bw()+
  geom_errorbar(aes(x = treat, ymax = hcarb_mgcm2+se, ymin = hcarb_mgcm2-se), width = .2, position = position_dodge(width=0.3)) +
  geom_point(size = 3, pch = 21, position = position_dodge(width=0.3))+
  scale_fill_manual(name = "Site",
                     values = cols_site)+
  xlab("Treatment")+
  ylab(bquote("Total Host Carbohydrate (mg" ~cm^-2~')'))+
  ylim(0,0.8) +
  geom_vline(xintercept = 1.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        legend.position = "none"
)
hcarb_plot_site

ggsave(hcarb_plot_site, filename = "/Users/hannahaichelman/Dropbox/BU/TVE/Carbohydrates/Plots/hcarb_site_2_lin.pdf", width=5, height=4, units=c("in"), useDingbats=FALSE)

#SummarySE to format data for plotting lineage
hcarb_means_lineage_2_lin <- summarySE(hcarb_phys_2_lin, measurevar="hcarb_mgcm2", groupvars=c("treat","lineage"))

# plot, treatment x axis colored by lineage data figure
hcarb_plot_lineage <- ggplot(hcarb_means_lineage_2_lin,aes(x = treat, y = hcarb_mgcm2, fill = lineage))+
  theme_bw()+
  geom_errorbar(aes(x = treat, ymax = hcarb_mgcm2+se, ymin = hcarb_mgcm2-se), width = .2, position = position_dodge(width=0.3)) +
  geom_point(size = 3, pch = 21, position = position_dodge(width=0.3))+
  scale_fill_manual(name = "Lineage",
                     values = cols_lineage)+
  xlab("Treatment")+
  ylab(bquote("Total Host Carbohydrate (mg" ~cm^-2~')'))+
  ylim(0,0.8) +
  geom_vline(xintercept = 1.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        legend.position = "none"
)
hcarb_plot_lineage

ggsave(hcarb_plot_lineage, filename = "/Users/hannahaichelman/Dropbox/BU/TVE/Carbohydrates/Plots/hcarb_lineage_2lin.pdf", width=5, height=4, units=c("in"), useDingbats=FALSE)


# plot, treatment x axis and symtype color
#SummarySE to format data for plotting
hcarb_means_sym <- summarySE(hcarb_phys_2_lin, measurevar="hcarb_mgcm2", groupvars=c("treat","dominant_type"))

# plot, treatment x axis colored by site data figure
hcarb_plot_sym <- ggplot(hcarb_means_sym,aes(x = treat, y = hcarb_mgcm2, fill = dominant_type))+
  theme_bw()+
  geom_errorbar(aes(x = treat, ymax = hcarb_mgcm2+se, ymin = hcarb_mgcm2-se), width = .2, position = position_dodge(width=0.3)) +
  geom_point(size = 3, pch = 21, position = position_dodge(width=0.3))+
  scale_fill_manual(name = "Dominant Symbiont",
                    values = its2_cols_greens)+
  xlab("Treatment")+
  ylab(bquote("Total Host Carbohydrate (mg" ~cm^-2~')'))+
  ylim(0,0.8) +
  geom_vline(xintercept = 1.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        legend.position = "none"
)
hcarb_plot_sym

ggsave(hcarb_plot_sym, filename = "/Users/hannahaichelman/Dropbox/BU/TVE/Carbohydrates/Plots/protein_lineage_sym.pdf", width=6, height=4, units=c("in"), useDingbats=FALSE)


# Combine these host carbohydrate plots
hcarb_plots_all = ggarrange(hcarb_plot_site, hcarb_plot_lineage, hcarb_plot_sym,
                           ncol = 3, nrow = 1)
ggsave(hcarb_plots_all, filename = "/Users/hannahaichelman/Dropbox/BU/TVE/Carbohydrates/Plots/hcarb_all_plots.pdf", width=12, height=4, units=c("in"), useDingbats=FALSE)



# plot just CI corals to see if lineage difference holds
hcarb_phys_2_lin_CI = hcarb_phys_2_lin %>%
  subset(sitename == "CI")

hcarb_means_lineage_CI <- summarySE(hcarb_phys_2_lin_CI, measurevar="hcarb_mgcm2", groupvars=c("treat","lineage"))

hcarb_plot_lineage_CI <- ggplot(hcarb_means_lineage_CI,aes(x = treat, y = hcarb_mgcm2, color = lineage, pch = lineage))+
  theme_bw()+
  geom_point(size = 3, position = position_dodge(width=0.3))+
  geom_errorbar(aes(x = treat, ymax = hcarb_mgcm2+se, ymin = hcarb_mgcm2-se), width = .2, position = position_dodge(width=0.3)) +
  scale_color_manual(name = "Lineage",
                     breaks = c("L1","L2"),
                     values = cols_lineage)+
  scale_shape_manual(name = "Lineage",
                     breaks = c("L1","L2"),
                     values=c(19,17))+
  xlab("Treatment")+
  ylab(bquote("Total Host Carbohydrate (mg" ~cm^-2~')'))+
  #ylim(0,0.4) +
  geom_vline(xintercept = 1.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))
hcarb_plot_lineage_CI

ggsave(hcarb_plot_lineage_CI, filename = "/Users/hannahaichelman/Dropbox/BU/TVE/Carbohydrates/Plots/hcarb_lineage_CIonly.pdf", width=5, height=4, units=c("in"), useDingbats=FALSE)


##### Sym Carb Concentrations #####
scarb_phys = phys %>%
  select(frag,treat,scarbplate_ha,scarb1_ha,scarb2_ha,scarb3_ha,
         survivedtoend,blastvol,gen_site,origsitecode,sitename,fragid,reef,genet,SAcm2,dominant_type,lineage)

head(scarb_phys)

# Add in standard curve info for each plate
scarb_phys$scarbxcoef <- as.numeric(ifelse(scarb_phys$scarbplate_ha == 'O18', 0.149980,
                                           ifelse(scarb_phys$scarbplate_ha == 'O19', 0.168115,
                                                  ifelse(scarb_phys$scarbplate_ha == 'O20', 0.1412873,
                                                         ifelse(scarb_phys$scarbplate_ha == 'O21', 0.142275,
                                                                ifelse(scarb_phys$scarbplate_ha == 'O22', 0.123317,
                                                                       ifelse(scarb_phys$scarbplate_ha == 'O23', 0.1539610,
                                                                              ifelse(scarb_phys$scarbplate_ha == 'O24', 0.1672786,
                                                                                     ifelse(scarb_phys$scarbplate_ha == 'O25', 0.2442510,
                                                                                            ifelse(scarb_phys$scarbplate_ha == 'O26', 0.123947,
                                                                                                   ifelse(scarb_phys$scarbplate_ha == 'O27', 0.133494,
                                                                                                          ifelse(scarb_phys$scarbplate_ha == 'O28', 0.136688,
                                                                                                                 ifelse(scarb_phys$scarbplate_ha == 'O29', 0.166897,
                                                                                                                        ifelse(scarb_phys$scarbplate_ha == 'O30', 0.1466,
                                                                                                                               ifelse(scarb_phys$scarbplate_ha == 'O31', 0.209305,
                                                                                                                                      ifelse(scarb_phys$scarbplate_ha == 'O32', 0.1308968,
                                                                                                                                             ifelse(scarb_phys$scarbplate_ha == 'O33', 0.153948,
                                                                                                                                                    ifelse(scarb_phys$scarbplate_ha == 'O34', 0.1407179,
                                                                                                                                                           ifelse(scarb_phys$scarbplate_ha == 'O35', 0.1536346,
                                                                                                                                                                  ifelse(scarb_phys$scarbplate_ha == 'O36', 0.1717878,
                                                                                                                                                                         'NA'))))))))))))))))))))

scarb_phys$scarbint <- as.numeric(ifelse(scarb_phys$scarbplate_ha == 'O18', -0.004060,
                                         ifelse(scarb_phys$scarbplate_ha == 'O19', 0.001898,
                                                ifelse(scarb_phys$scarbplate_ha == 'O20', -0.0003309,
                                                       ifelse(scarb_phys$scarbplate_ha == 'O21', 0.001065,
                                                              ifelse(scarb_phys$scarbplate_ha == 'O22', -0.001284,
                                                                     ifelse(scarb_phys$scarbplate_ha == 'O23', -0.0005789,
                                                                            ifelse(scarb_phys$scarbplate_ha == 'O24', 0.0012109,
                                                                                   ifelse(scarb_phys$scarbplate_ha == 'O25', -0.0003185,
                                                                                          ifelse(scarb_phys$scarbplate_ha == 'O26', 0.002024,
                                                                                                 ifelse(scarb_phys$scarbplate_ha == 'O27', -0.000288,
                                                                                                        ifelse(scarb_phys$scarbplate_ha == 'O28', -0.001351,
                                                                                                               ifelse(scarb_phys$scarbplate_ha == 'O29', -0.004318,
                                                                                                                      ifelse(scarb_phys$scarbplate_ha == 'O30', -0.00004093,
                                                                                                                             ifelse(scarb_phys$scarbplate_ha == 'O31', -0.001261,
                                                                                                                                    ifelse(scarb_phys$scarbplate_ha == 'O32', -0.0001588,
                                                                                                                                           ifelse(scarb_phys$scarbplate_ha == 'O33', -0.005824,
                                                                                                                                                  ifelse(scarb_phys$scarbplate_ha == 'O34', 0.000427,
                                                                                                                                                         ifelse(scarb_phys$scarbplate_ha == 'O35', -0.0004174,
                                                                                                                                                                ifelse(scarb_phys$scarbplate_ha == 'O36', 0.0003559,
                                                                                                                                                                       'NA'))))))))))))))))))))

# get the carbohydrate values in mg/mL (concentration of D-glucose standard)
# the standard curve equation for carbs is linear - y = m*x + b -
# where y is the calculated concentration in mg/mL, m is the slope, or x coefficient, and b is the intercept
scarb_phys$calcscarb1 <- ((scarb_phys$scarbxcoef*scarb_phys$scarb1_ha) + scarb_phys$scarbint)
scarb_phys$calcscarb2 <- ((scarb_phys$scarbxcoef*scarb_phys$scarb2_ha) + scarb_phys$scarbint)
scarb_phys$calcscarb3 <- ((scarb_phys$scarbxcoef*scarb_phys$scarb3_ha) + scarb_phys$scarbint)

# calculate average protein concentration, mg/mL
scarb_phys = scarb_phys %>%
  mutate(avgscarb=rowMeans(.[ , c("calcscarb1","calcscarb2","calcscarb3")], na.rm=TRUE))

# [carb concentration in in mg/cm2] = ([mean concentration mg/mL] * [dilution factor] * number of mls in slurry / [surface area of fragment]
# dilution factor is 5 here because in the protocol, each reaction has 10 uL of sample and 40 uL of seawater
scarb_phys$scarb_mgcm2 <- (scarb_phys$avgscarb*5*scarb_phys$blastvol)/scarb_phys$SAcm2

head(scarb_phys)

# make data subsets for stats and plotting
scarb_phys_all_lin = scarb_phys %>%
  drop_na(scarb_mgcm2) %>%
  dplyr::filter(treat!="Control 2") %>%
  dplyr::filter(treat!="Low Var") %>%
  dplyr::filter(treat!="High Var") %>%
  dplyr::filter(gen_site!="I4G") %>% # remove clone
  select(frag, gen_site, treat, sitename, dominant_type, lineage, scarb_mgcm2)

scarb_phys_2_lin = scarb_phys_all_lin %>%
  dplyr::filter(is.na(lineage) | lineage!="L3") %>%
  drop_na(sitename, lineage, treat, dominant_type) # doing this for plotting

initial_scarb_phys_all_lin = scarb_phys_all_lin %>%
  dplyr::filter(treat=="Initial")

initial_scarb_phys_2_lin = initial_scarb_phys_all_lin %>%
  dplyr::filter(is.na(lineage) | lineage!="L3") # want to keep NA values for lineage here since they still have other info, will remove na's for lineage specific plots

post_scarb_phys_all_lin = scarb_phys_all_lin %>%
  dplyr::filter(treat!="Initial")

post_scarb_phys_2_lin = post_scarb_phys_all_lin %>%
  dplyr::filter(is.na(lineage) | lineage!="L3")


# STATS
## NEED TO RUN STATS SEPARATELY FOR INITIAL AND FINAL TIME POINTS
# use lmer()
model.scarb.initial = lm(scarb_mgcm2 ~ lineage, data = initial_scarb_phys_2_lin)
summary(model.scarb.initial)
anova(model.scarb.initial)
# Response: scarb_mgcm2
#           Df   Sum Sq   Mean Sq F value   Pr(>F)
# lineage    1 0.021801 0.0218014    9.01 0.004556 **
# Residuals 41 0.099208 0.0024197

check_model(model.scarb.initial)


model.scarb.final <- lmer(scarb_mgcm2 ~ lineage+treat + (1|gen_site), data = post_scarb_phys_2_lin)
summary(model.scarb.final)
anova(model.scarb.final)
# Type III Analysis of Variance Table with Satterthwaite's method
#           Sum Sq  Mean Sq NumDF  DenDF F value    Pr(>F)
# lineage 0.012491 0.012491     1 46.485  3.1981 0.0802413 .
# treat   0.048970 0.048970     1 45.600 12.5380 0.0009318 ***

check_model(model.scarb.final)


#specify model (because we are interested in pairwise, have to include the interaction)
m.emm<- lmer(scarb_mgcm2 ~ treat*lineage + (1|gen_site), data = post_scarb_phys_2_lin, REML=FALSE)

emms<-emmeans(m.emm, ~lineage|treat) #, adjust="Bonferoni"
pairs(emms, interaction = "pairwise") %>% rbind(adjust="fdr")
# treat   lineage_pairwise estimate     SE   df t.ratio p.value
# Control L1 - L2           0.05020 0.0223 83.0   2.252  0.0540
# Mod Var L1 - L2           0.00884 0.0208 82.8   0.425  0.6718

# PLOTS
#SummarySE to format data for plotting
scarb_means_site_2_lin <- summarySE(scarb_phys_2_lin, measurevar="scarb_mgcm2", groupvars=c("treat","sitename"))

# plot, treatment x axis colored by site data figure
scarb_plot_site <- ggplot(scarb_means_site_2_lin,aes(x = treat, y = scarb_mgcm2, fill = sitename))+
  theme_bw()+
  geom_errorbar(aes(x = treat, ymax = scarb_mgcm2+se, ymin = scarb_mgcm2-se), width = .2, position = position_dodge(width=0.3)) +
  geom_point(size = 3, pch = 21, position = position_dodge(width=0.3))+
  scale_fill_manual(name = "Site",
                    values = cols_site)+
  xlab("Treatment")+
  ylab(bquote("Total Symbiont Carbohydrate (mg" ~cm^-2~')'))+
  ylim(0,0.35) +
  geom_vline(xintercept = 1.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        legend.position = "none"
)
scarb_plot_site

ggsave(scarb_plot_site, filename = "/Users/hannahaichelman/Dropbox/BU/TVE/Carbohydrates/Plots/scarb_site_2_lin.pdf", width=5, height=4, units=c("in"), useDingbats=FALSE)

#SummarySE to format data for plotting lineage
scarb_means_lineage_2_lin <- summarySE(scarb_phys_2_lin, measurevar="scarb_mgcm2", groupvars=c("treat","lineage"))

# plot, treatment x axis colored by lineage data figure
scarb_plot_lineage <- ggplot(scarb_means_lineage_2_lin,aes(x = treat, y = scarb_mgcm2, fill = lineage))+
  theme_bw()+
  geom_errorbar(aes(x = treat, ymax = scarb_mgcm2+se, ymin = scarb_mgcm2-se), width = .2, position = position_dodge(width=0.3)) +
  geom_point(size = 3, pch = 21, position = position_dodge(width=0.3))+
  scale_fill_manual(name = "Lineage",
                     values = cols_lineage)+
  xlab("Treatment")+
  ylab(bquote("Total Symbiont Carbohydrate (mg" ~cm^-2~')'))+
  ylim(0,0.35) +
  geom_vline(xintercept = 1.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        legend.position = "none"
)

scarb_plot_lineage

ggsave(scarb_plot_lineage, filename = "/Users/hannahaichelman/Dropbox/BU/TVE/Carbohydrates/Plots/scarb_lineage_2lin.pdf", width=5, height=4, units=c("in"), useDingbats=FALSE)


# plot, treatment x axis and symtype color
#SummarySE to format data for plotting
scarb_means_sym <- summarySE(scarb_phys_2_lin, measurevar="scarb_mgcm2", groupvars=c("treat","dominant_type"))

# plot, treatment x axis colored by site data figure
scarb_plot_sym <- ggplot(scarb_means_sym,aes(x = treat, y = scarb_mgcm2, fill = dominant_type))+
  theme_bw()+
  geom_errorbar(aes(x = treat, ymax = scarb_mgcm2+se, ymin = scarb_mgcm2-se), width = .2, position = position_dodge(width=0.3)) +
  geom_point(size = 3, pch = 21, position = position_dodge(width=0.3))+
  scale_fill_manual(name = "Dominant Symbiont",
                    values = its2_cols_greens)+
  xlab("Treatment")+
  ylab(bquote("Total Symbiont Carbohydrate (mg" ~cm^-2~')'))+
  ylim(0,0.35) +
  geom_vline(xintercept = 1.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        legend.position = "none"
)
scarb_plot_sym

ggsave(scarb_plot_sym, filename = "/Users/hannahaichelman/Dropbox/BU/TVE/Carbohydrates/Plots/scarb_lineage_sym.pdf", width=6, height=4, units=c("in"), useDingbats=FALSE)


# Combine these symbiont carbohydrate plots
scarb_plots_all = ggarrange(scarb_plot_site, scarb_plot_lineage, scarb_plot_sym,
                            ncol = 3, nrow = 1)
ggsave(scarb_plots_all, filename = "/Users/hannahaichelman/Dropbox/BU/TVE/Carbohydrates/Plots/scarb_all_plots.pdf", width=12, height=4, units=c("in"), useDingbats=FALSE)


##### Symbiont Density #####
sym_phys = phys %>%
  select(frag,treat,symcount1,symcount2,symcount3,
         survivedtoend,blastvol,gen_site,origsitecode,sitename,fragid,reef,genet,SAcm2,dominant_type,lineage)

head(sym_phys)

# Calculate average per surface area
# make data subsets for stats and plotting
sym_phys_all_lin = sym_phys %>%
  mutate(avgsymcount=rowMeans(.[ , c("symcount1","symcount2","symcount3")], na.rm=TRUE)) %>%
  mutate(sym_cellsmL=avgsymcount*10000) %>% # this is based on the following website http://insilico.ehu.eus/counting_chamber/neubauer_improved.php
  mutate(sym_cm2=(sym_cellsmL*blastvol)/SAcm2) %>%
  mutate(sym_cm2_div = sym_cm2/1000000) %>%
  dplyr::filter(treat!="Control 2") %>%
  dplyr::filter(treat!="Low Var") %>%
  dplyr::filter(treat!="High Var") %>%
  dplyr::filter(frag!="I3E6_v2") %>% # same data for sym count for original and v2 of I3E6
  dplyr::filter(frag!="I3A4") %>%
  dplyr::filter(gen_site != "I4G") %>% # remove clone
  drop_na(sym_cm2) %>%
  select(frag, gen_site, treat, sitename, dominant_type, lineage, sym_cm2, sym_cm2_div)


sym_phys_2_lin = sym_phys_all_lin %>%
  dplyr::filter(is.na(lineage) | lineage!="L3") %>%
  drop_na(sitename, lineage, treat, dominant_type) # doing this for plotting

initial_sym_phys_all_lin = sym_phys_all_lin %>%
  dplyr::filter(treat=="Initial")

initial_sym_phys_2_lin = initial_sym_phys_all_lin %>%
  dplyr::filter(is.na(lineage) | lineage!="L3") # want to keep NA values for lineage here since they still have other info, will remove na's for lineage specific plots

post_sym_phys_all_lin = sym_phys_all_lin %>%
  dplyr::filter(treat!="Initial")

post_sym_phys_2_lin = post_sym_phys_all_lin %>%
  dplyr::filter(is.na(lineage) | lineage!="L3")


# STATS
## NEED TO RUN STATS SEPARATELY FOR INITIAL AND FINAL TIME POINTS
# use lmer()
model.sym.initial = lm(sym_cm2 ~ lineage, data = initial_sym_phys_2_lin)
summary(model.sym.initial)
anova(model.sym.initial)
#           Df     Sum Sq    Mean Sq F value    Pr(>F)
# lineage    1 9.3732e+12 9.3732e+12  18.371 0.0001074 ***
# Residuals 41 2.0919e+13 5.1022e+11

check_model(model.sym.initial)


model.sym.final <- lmer(sym_cm2 ~ lineage+treat + (1|gen_site), data = post_sym_phys_2_lin)
summary(model.sym.final)
anova(model.sym.final)
# Type III Analysis of Variance Table with Satterthwaite's method
#             Sum Sq    Mean Sq NumDF  DenDF F value   Pr(>F)
# lineage 3.0830e+12 3.0830e+12     1   40.8 11.7416 0.001408 **
# treat   1.7649e+12 1.7649e+12     1 4967.4  6.7216 0.009553 **

check_model(model.sym.final)


#specify model (because we are interested in pairwise, have to include the interaction)
m.emm<- lmer(sym_cm2 ~ treat*lineage + (1|gen_site), data = post_sym_phys_2_lin, REML=FALSE)

emms<-emmeans(m.emm, ~lineage|treat) #, adjust="Bonferoni"
pairs(emms, interaction = "pairwise") %>% rbind(adjust="fdr")
# treat   lineage_pairwise estimate     SE   df t.ratio p.value
# Control L1 - L2            495566 188601 83.0   2.628  0.0198
# Mod Var L1 - L2            417690 175851 82.8   2.375  0.0198


# PLOTS
#SummarySE to format data for plotting
sym_means_site_2_lin <- summarySE(sym_phys_2_lin, measurevar="sym_cm2_div", groupvars=c("treat","sitename"))

# plot, treatment x axis colored by site data figure
sym_plot_site <- ggplot(sym_means_site_2_lin,aes(x = treat, y = sym_cm2_div, fill = sitename))+
  theme_bw()+
  geom_errorbar(aes(x = treat, ymax = sym_cm2_div+se, ymin = sym_cm2_div-se), width = .2, position = position_dodge(width=0.3)) +
  geom_point(size = 3, pch = 21, position = position_dodge(width=0.3))+
  scale_fill_manual(name = "Site",
                    values = cols_site)+
  xlab("Treatment")+
  ylab(bquote("Symbiont Density ("~x10^6~ 'cells' ~cm^-2~')'))+
  ylim(0,2.55) +
  geom_vline(xintercept = 1.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        legend.position = "none"
)
sym_plot_site

ggsave(sym_plot_site, filename = "/Users/hannahaichelman/Dropbox/BU/TVE/SymbiontDensity/plots/sym_site_2_lin.pdf", width=5, height=4, units=c("in"), useDingbats=FALSE)

#SummarySE to format data for plotting
sym_means_lineage_2_lin <- summarySE(sym_phys_2_lin, measurevar="sym_cm2_div", groupvars=c("treat","lineage"))

# plot, treatment x axis colored by lineage data figure
sym_plot_lineage <- ggplot(sym_means_lineage_2_lin,aes(x = treat, y = sym_cm2_div, fill = lineage))+
  theme_bw()+
  geom_errorbar(aes(x = treat, ymax = sym_cm2_div+se, ymin = sym_cm2_div-se), width = .2, position = position_dodge(width=0.3)) +
  geom_point(size = 3, pch = 21, position = position_dodge(width=0.3))+
  scale_fill_manual(name = "Lineage",
                    values = cols_lineage)+
  xlab("Treatment")+
  ylab(bquote("Symbiont Density ("~x10^6~ 'cells' ~cm^-2~')'))+
  ylim(0,2.55) +
  geom_vline(xintercept = 1.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        legend.position = "none"
  )
sym_plot_lineage

ggsave(sym_plot_lineage, filename = "/Users/hannahaichelman/Dropbox/BU/TVE/SymbiontDensity/plots/sym_lineage_2_lin.pdf", width=5, height=4, units=c("in"), useDingbats=FALSE)


# plot, treatment x axis and symtype color
#SummarySE to format data for plotting
sym_means_sym <- summarySE(sym_phys_2_lin, measurevar="sym_cm2_div", groupvars=c("treat","dominant_type"))

# plot, treatment x axis colored by site data figure
sym_plot_sym <- ggplot(sym_means_sym,aes(x = treat, y = sym_cm2_div, fill = dominant_type))+
  theme_bw()+
  geom_errorbar(aes(x = treat, ymax = sym_cm2_div+se, ymin = sym_cm2_div-se), width = .2, position = position_dodge(width=0.3)) +
  geom_point(size = 3, pch = 21, position = position_dodge(width=0.3))+
  scale_fill_manual(name = "Dominant Symbiont",
                    values = its2_cols_greens)+
  xlab("Treatment")+
  ylab(bquote("Symbiont Density ("~x10^6~ 'cells' ~cm^-2~')'))+
  ylim(0,2.55) +
  geom_vline(xintercept = 1.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        legend.position = "none"
  )
sym_plot_sym

ggsave(sym_plot_sym, filename = "/Users/hannahaichelman/Dropbox/BU/TVE/SymbiontDensity/plots/sym_lineage_2_lin_symtype.pdf", width=6, height=4, units=c("in"), useDingbats=FALSE)


# Combine these symbiont density plots
sym_plots_all = ggarrange(sym_plot_site, sym_plot_lineage, sym_plot_sym,
                            ncol = 3, nrow = 1)
ggsave(sym_plots_all, filename = "/Users/hannahaichelman/Dropbox/BU/TVE/SymbiontDensity/plots/sym_all_plots.pdf", width=12, height=4, units=c("in"), useDingbats=FALSE)



## Checking for outliers - used this to remove I3A4
#plot, facet by treatment
plot <- ggplot(sym_phys_nona,aes(x = sitename, y = sym_cm2, color = sitename, pch = sitename, text = paste("Nubbin:",frag)))+
  theme_bw()+
  geom_boxplot()+
  ggtitle("Symbiont Density (cells/cm2)")+
  scale_fill_manual(values = cols_site)+
  geom_jitter(aes(color = sitename),width = 0.3)+
  scale_color_manual(values = cols_site)+
  theme(axis.title.x = element_blank())+
  theme(axis.text.x=element_text(angle = 45, vjust = 0.5))+
  facet_wrap(~treat, nrow=5)
#  labs(y=expression(paste("Protein (mg cm"^-2*' )')))
ggplotly(subplot(list(plot),titleY=F) %>% layout(showlegend=F))

# facet by genotype
plot <- ggplot(phys,aes(x = treat, y = sym_cm2, color = sitename, pch = sitename, text = paste("Nubbin:",frag)))+
  theme_bw()+
  geom_jitter(width = 0.3, aes(fill = sitename))+
  scale_fill_manual(values = palsite)+
  scale_color_manual(values = palsite)+
  theme(axis.title.x = element_blank())+
  theme(axis.text.x=element_text(angle = 45, vjust = 0.5))+
  facet_wrap(~genet+origsitecode, nrow=4)+
  ggtitle("Symbiont Density (cells/cm2)")
#  labs(y=expression(paste("Protein (mg cm"^-2*' )')))
ggplotly(subplot(list(plot),titleY=T) %>% layout(showlegend=T))


##### Chlorophyll Concentration #####

## read in additional chlorophyll data that Olivia added
chl <- read.csv('Physiology_Data/data_files/TVEChlorphyll_ON.csv')

phys_for_chl = phys %>%
  select(frag, survivedtoend, treat, blastvol, blaster, SAcm2, reef, sitename, gen_site, dominant_type, lineage) %>%
  mutate(treat = as.factor(treat), blastvol = as.numeric(blastvol))

#want to keep all of the individuals in the original dtv spreadsheet and add NA's for any that don't have chl info for
chl_phys <- left_join(phys_for_chl, chl, by = "frag")

# find mean absorbance values for each wavelength
chl_phys = chl_phys %>%
  mutate(chl630_avg=rowMeans(.[ , c("blk_630_1","blk_630_2","blk_630_3")], na.rm=TRUE)) %>%
  mutate(chl663_avg=rowMeans(.[ , c("blk_663_1","blk_663_2","blk_663_3")], na.rm=TRUE))

# calculate chlorophyll A and C2 from values in ug/mL (equations come from Carly Kenkel's method, originally a Jeffrey and Haxo 1968 equation)
# 20 cancels out mL so you end up with ug/cm2
chl_phys$chlA <- ((((13.31*chl_phys$chl663_avg) - (0.27*chl_phys$chl630_avg))*20)/chl_phys$SAcm2)
chl_phys$chlC2 <- ((((-8.37*chl_phys$chl663_avg) + (51.72*chl_phys$chl630_avg))*20)/chl_phys$SAcm2)

head(chl_phys)

# make data subsets for stats and plotting
chl_phys_all_lin = chl_phys %>%
  drop_na(chlA) %>%
  dplyr::filter(treat!="Control 2") %>%
  dplyr::filter(treat!="Low Var") %>%
  dplyr::filter(treat!="High Var") %>%
  dplyr::filter(frag!="I2A12", frag!="I2A7", frag!="O4D8", frag!="I3E10", frag!="I3I10", frag!="I3E6_v2") %>% #these frags are being removed because they were duplicate genotypes within treatment and had the most complete information of the two
  dplyr::filter(frag!="I3A4") %>% # unexplained high outlier
  dplyr::filter(gen_site != "I4G") %>% # clone with I4F, remove from dataset
  select(frag, gen_site, treat, sitename, dominant_type, lineage, chlA, chlC2)

chl_phys_2_lin = chl_phys_all_lin %>%
  dplyr::filter(is.na(lineage) | lineage!="L3") %>%
  drop_na(sitename, lineage, treat, dominant_type) # doing this for plotting

initial_chl_phys_all_lin = chl_phys_all_lin %>%
  dplyr::filter(treat=="Initial")

initial_chl_phys_2_lin = initial_chl_phys_all_lin %>%
  dplyr::filter(is.na(lineage) | lineage!="L3") # want to keep NA values for lineage here since they still have other info, will remove na's for lineage specific plots

post_chl_phys_all_lin = chl_phys_all_lin %>%
  dplyr::filter(treat!="Initial")

post_chl_phys_2_lin = post_chl_phys_all_lin %>%
  dplyr::filter(is.na(lineage) | lineage!="L3")


# STATS
## NEED TO RUN STATS SEPARATELY FOR INITIAL AND FINAL TIME POINTS
# use lmer()
model.chl.initial = lm(chlA ~ lineage, data = initial_chl_phys_2_lin)
summary(model.chl.initial)
anova(model.chl.initial)
#           Df Sum Sq Mean Sq F value    Pr(>F)
# lineage    1 297.51 297.505   20.45 5.138e-05 ***
# Residuals 41 596.46  14.548

check_model(model.chl.initial)


model.chl.final <- lmer(chlA ~ treat+lineage + (1|gen_site), data = post_chl_phys_2_lin)
summary(model.chl.final)
anova(model.chl.final)
# Type III Analysis of Variance Table with Satterthwaite's method
#          Sum Sq Mean Sq NumDF DenDF F value   Pr(>F)
# treat    0.0048  0.0048     1    78  0.0016 0.968064
# lineage 29.0579 29.0579     1    78  9.7166 0.002557 **

check_model(model.chl.final)


#specify model (because we are interested in pairwise, have to include the interaction)
m.emm<- lmer(chlA ~ treat*lineage + (1|gen_site), data = post_chl_phys_2_lin, REML=FALSE)

emms<-emmeans(m.emm, ~lineage|treat) #, adjust="Bonferoni"
pairs(emms, interaction = "pairwise") %>% rbind(adjust="fdr")
# treat   lineage_pairwise estimate    SE   df t.ratio p.value
# Control L1 - L2             1.561 0.601 85.2   2.597  0.0222
# Mod Var L1 - L2             0.983 0.537 85.2   1.829  0.0709

# PLOTS
#SummarySE to format data for plotting
chla_means_site_2_lin <- summarySE(chl_phys_2_lin, measurevar="chlA", groupvars=c("treat","sitename"))

# plot, treatment x axis colored by site data figure
chla_plot_site <- ggplot(chla_means_site_2_lin,aes(x = treat, y = chlA, fill = sitename))+
  theme_bw()+
  geom_errorbar(aes(x = treat, ymax = chlA+se, ymin = chlA-se), width = .2, position = position_dodge(width=0.3)) +
  geom_point(size = 3, pch = 21, position = position_dodge(width=0.3))+
  scale_fill_manual(name = "Site",
                    values = cols_site)+
  xlab("Treatment")+
  ylab(bquote("Chl a ("*mu*"g" ~cm^-2~')'))+
  ylim(0,15) +
  geom_vline(xintercept = 1.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        legend.position = "none")
chla_plot_site

ggsave(chla_plot_site, filename = "/Users/hannahaichelman/Dropbox/BU/TVE/Chlorophyll/plots/chla_site_2_lin.pdf", width=5, height=4, units=c("in"), useDingbats=FALSE)

#SummarySE to format data for plotting
chla_means_lineage_2_lin <- summarySE(chl_phys_2_lin, measurevar="chlA", groupvars=c("treat","lineage"))

# plot, treatment x axis colored by lineage data figure
chla_plot_lineage <- ggplot(chla_means_lineage_2_lin,aes(x = treat, y = chlA, fill = lineage))+
  theme_bw()+
  geom_errorbar(aes(x = treat, ymax = chlA+se, ymin = chlA-se), width = .2, position = position_dodge(width=0.3)) +
  geom_point(size = 3, pch = 21, position = position_dodge(width=0.3))+
  scale_fill_manual(name = "Lineage",
                    values = cols_lineage)+
  xlab("Treatment")+
  ylab(bquote("Chl a ("*mu*"g" ~cm^-2~')'))+
  ylim(0,15) +
  geom_vline(xintercept = 1.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        legend.position = "none")
chla_plot_lineage

ggsave(chla_plot_lineage, filename = "/Users/hannahaichelman/Dropbox/BU/TVE/Chlorophyll/plots/chla_lineage_2_lin.pdf", width=5, height=4, units=c("in"), useDingbats=FALSE)

# plot, treatment x axis and symtype color
#SummarySE to format data for plotting
chla_means_sym <- summarySE(chl_phys_2_lin, measurevar="chlA", groupvars=c("treat","dominant_type"))

# plot, treatment x axis colored by site data figure
chla_plot_sym <- ggplot(chla_means_sym,aes(x = treat, y = chlA, fill = dominant_type))+
  theme_bw()+
  geom_errorbar(aes(x = treat, ymax = chlA+se, ymin = chlA-se), width = .2, position = position_dodge(width=0.3)) +
  geom_point(size = 3, pch = 21, position = position_dodge(width=0.3))+
  scale_fill_manual(name = "Dominant Symbiont",
                    values = its2_cols_greens)+
  xlab("Treatment")+
  ylab(bquote("Chl a ("*mu*"g" ~cm^-2~')'))+
  ylim(0,15) +
  geom_vline(xintercept = 1.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        legend.position = "none")
chla_plot_sym

ggsave(chla_plot_sym, filename = "/Users/hannahaichelman/Dropbox/BU/TVE/Chlorophyll/plots/chla_symtype_2_lin.pdf", width=6, height=4, units=c("in"), useDingbats=FALSE)



# Combine these chlorophyll A plots
chla_plots_all = ggarrange(chla_plot_site, chla_plot_lineage, chla_plot_sym,
                           ncol = 3, nrow = 1)
ggsave(chla_plots_all, filename = "/Users/hannahaichelman/Dropbox/BU/TVE/Chlorophyll/plots/chla_all_plots.pdf", width=12, height=4, units=c("in"), useDingbats=FALSE)



# plot just CI corals to see if lineage difference holds
chl_phys_lineage_CI = chl_phys_2_lin %>%
  subset(sitename == "CI")

chla_means_lineage_CI <- summarySE(chl_phys_lineage_CI, measurevar="chlA", groupvars=c("treat","lineage"))

chla_plot_lineage_CI <- ggplot(chla_means_lineage_CI,aes(x = treat, y = chlA, color = lineage, pch = lineage))+
  theme_bw()+
  geom_point(size = 3, position = position_dodge(width=0.3))+
  geom_errorbar(aes(x = treat, ymax = chlA+se, ymin = chlA-se), width = .2, position = position_dodge(width=0.3)) +
  scale_color_manual(name = "Lineage",
                     breaks = c("L1","L2"),
                     values = cols_lineage)+
  scale_shape_manual(name = "Lineage",
                     breaks = c("L1","L2"),
                     values=c(19,17))+
  xlab("Treatment")+
  ylab("Chlorophyll A")+
  #ylim(0,0.4) +
  geom_vline(xintercept = 1.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))
chla_plot_lineage_CI

ggsave(chla_plot_lineage_CI, filename = "/Users/hannahaichelman/Dropbox/BU/TVE/Chlorophyll/plots/chla_lineage_CIonly.pdf", width=5, height=4, units=c("in"), useDingbats=FALSE)



# look for outliers with ggplotly
# chl A
# violin plot with individual data point overlay
plot <- ggplot(chl_phys,aes(x = sitename, y = chlA, color = sitename, pch = sitename,text = paste("Frag:",frag)))+
  theme_bw()+
  geom_boxplot(alpha = 0.2, aes(fill = sitename))+
  scale_fill_manual(values = cols_site)+
  geom_jitter(aes(color = sitename),width = 0.3)+
  ggtitle("Chl A")+
  scale_color_manual(values = cols_site)+
  theme(axis.title.x = element_blank())+
  theme(axis.text.x=element_text(angle = 45, vjust = 0.5))
#  labs(y=expression(paste("Chl A (ug cm"^-2*' )')))
ggplotly(subplot(list(plot),nrows=1,titleY=T) %>% layout(showlegend=T))

# chl C2
plot <- ggplot(chl_phys,aes(x = sitename, y = chlC2, color = sitename, pch = sitename,text = paste("Frag:",frag)))+
  theme_bw()+
  geom_boxplot(alpha = 0.2, aes(fill = sitename))+
  scale_fill_manual(values = palsite)+
  geom_jitter(aes(color = sitename),width = 0.3)+
  ggtitle("Chl C2")+
  scale_color_manual(values = palsite)+
  theme(axis.title.x = element_blank())+
  theme(axis.text.x=element_text(angle = 45, vjust = 0.5))
#  labs(y=expression(paste("Chl C2 (ug cm"^-2*' )')))
ggplotly(subplot(list(plot),nrows=1,titleY=T) %>% layout(showlegend=T))

# there are lots of negative chlorophyll c2 values...seems like something was off with these readings so omitting from the manuscript.

##### Growth #####
# need to re-read in the data sheet for calcification only
post_phys_forcalc <- read.csv('Physiology_Data/data_files/dtvmaster.csv')

# remove unexplained outliers, and subset by only the data we need here
calc_phys = post_phys_forcalc %>%
  select(frag, survivedtoend, treat, blastvol, blaster,
         t2sastan1, t2sastan2, t2sastan3, t2sarec1, t2sarec2, t2sarec3,
         t0sastan1, t0sastan2, t0sastan3, t0sarec1, t0sarec2, t0sarec3,
         t3sastan1, t3sastan2, t3sastan3, t3sarec1, t3sarec2, t3sarec3,
         pabwrec1,pabwrec2,pabwrec3,t0bwrec1, t0bwrec2,t0bwrec3,t1bwrec1, t1bwrec2,t1bwrec3,t2bwrec1, t2bwrec2,t2bwrec3,t3bwrec1, t3bwrec2,t3bwrec3) %>%
  mutate(treat = as.factor(treat), blastvol = as.numeric(blastvol)) %>%
  mutate_at(c(6:23), as.numeric) %>%
  dplyr::filter(frag!="I2A12", frag!="I2A7", frag!="O4D8", frag!="I3E10", frag!="I3I10") %>% #these frags are being removed because they were duplicate genotypes within treatment and had the most complete information of the two
  dplyr::filter(frag!="O2F1", frag!="O4G7", frag!="O2I3", frag!="O3H3", frag!="I3D8", frag!="O3F4", frag!="O3G10", frag!="O3F3", frag!="I4F8") #these frags are being removed because they were unexplainable outliers with weird weight values

#calculate surface area of the coral fragments
calc_phys$SAcm2_t0 <- (4*((calc_phys$t0sarec1+calc_phys$t0sarec2+calc_phys$t0sarec3)/3))/((calc_phys$t0sastan1+calc_phys$t0sastan2+calc_phys$t0sastan3)/3)
calc_phys$SAcm2_t2 <- (4*((calc_phys$t2sarec1+calc_phys$t2sarec2+calc_phys$t2sarec3)/3))/((calc_phys$t2sastan1+calc_phys$t2sastan2+calc_phys$t2sastan3)/3)
calc_phys$SAcm2_t3 <- (4*((calc_phys$t3sarec1+calc_phys$t3sarec2+calc_phys$t3sarec3)/3))/((calc_phys$t3sastan1+calc_phys$t3sastan2+calc_phys$t3sastan3)/3)

#add in descriptive information for coral samples
calc_phys$origsitecode <- substr(calc_phys$frag, 1, 2)
calc_phys$origsitecode <- as.factor(calc_phys$origsitecode)

# add in descriptive site name with inshore/offshore indicator
calc_phys$sitename <- ifelse(calc_phys$origsitecode == 'I2', 'SP',
                             ifelse(calc_phys$origsitecode == 'I3', 'CI',
                                    ifelse(calc_phys$origsitecode == 'I4', 'PD',
                                           ifelse(calc_phys$origsitecode == 'O2', 'BS',
                                                  ifelse(calc_phys$origsitecode == 'O3', 'CA',
                                                         'BN')))))
calc_phys$sitename <- as.factor(calc_phys$sitename)

# extract the genotype and fragment number
calc_phys$fragid <- substr(calc_phys$frag,3,5)

# add inshore/offshore designation
calc_phys$reef <- substr(calc_phys$frag,1,1)
calc_phys$reef <- ifelse(calc_phys$reef == 'O', 'Outer Reef', 'Inner Reef')
calc_phys$reef <- as.factor(calc_phys$reef)

# add in genotype
calc_phys$genet <- substr(calc_phys$fragid,1,1)

#create a new column of combined genotype that will be used as a random effect in stats later
calc_phys = calc_phys %>%
  unite(gen_site, c(origsitecode,genet), sep = "", remove = FALSE) %>%
  mutate(gen_site = as.factor(gen_site))

# re-level and re-name treatment
calc_phys$treat <- factor(calc_phys$treat, levels = c("init", "1", "2", "3","4","5"))
levels(calc_phys$treat) <- c("Initial","Control","Low Var","Mod Var","High Var","Control 2")

# calculate average weight for each time point
# want mg/cm2/day (weight was recorded in grams)
# PA buoyant weight on 9.6.2016 & 9.7.2016
# T0 buoyant weight on 9.22.2016 & 9.23.2016
# T1 buoyant weight on 10.17.2016 & 10.18.2016
# T2 buoyant weight on 11.8.2016 & 11.9.2016 (taken just before heat stress)
# T3 buoyant weight on 12.10.2016 & 12.11.2016 (taken after stress, at end of recovery)
# T0-PA = 16 days, T1-T0 = 25 days, T2-T1 = 22 days, T3-T2 = 32 days
# T3-T0 =

# calculate average weight at each time point
calc_phys = calc_phys %>%
  mutate(paavgbw=rowMeans(.[ , c("pabwrec1","pabwrec2","pabwrec3")], na.rm=TRUE)) %>%
  mutate(t0avgbw=rowMeans(.[ , c("t0bwrec1","t0bwrec2","t0bwrec3")], na.rm=TRUE)) %>%
  mutate(t1avgbw=rowMeans(.[ , c("t1bwrec1","t1bwrec2","t1bwrec3")], na.rm=TRUE)) %>%
  mutate(t2avgbw=rowMeans(.[ , c("t2bwrec1","t2bwrec2","t2bwrec3")], na.rm=TRUE)) %>%
  mutate(t3avgbw=rowMeans(.[ , c("t3bwrec1","t3bwrec2","t3bwrec3")], na.rm=TRUE))

# calculate growth rate at the time point intervals - can't use this though because we do not have dry weights
calc_phys = calc_phys %>%
  mutate(T0_PA_g_cm2_day=((t0avgbw-paavgbw)/16)/SAcm2_t0) %>%
  mutate(T2_T0_g_cm2_day=((t2avgbw-t0avgbw)/57)/SAcm2_t2) %>%
  mutate(T2_T1_g_cm2_day=((t2avgbw-t1avgbw)/22)/SAcm2_t2) %>%
  mutate(T3_T2_g_cm_day2=((t3avgbw-t2avgbw)/32)/SAcm2_t3) %>%
  mutate(T3_T0_g_cm_day2=((t3avgbw-t0avgbw)/79)/SAcm2_t3)

# calculate percent change in growth
calc_phys = calc_phys %>%
  mutate(T0_PA_perc=((t0avgbw-paavgbw)/paavgbw)*100) %>%
  mutate(T1_T0_perc=((t1avgbw-t0avgbw)/t0avgbw)*100) %>%
  mutate(T2_T1_perc=((t2avgbw-t1avgbw)/t1avgbw)*100) %>%
  mutate(T3_T2_perc=((t3avgbw-t2avgbw)/t2avgbw)*100) %>%  # this is from before to after stress
  mutate(T2_T0_perc=((t2avgbw-t0avgbw)/t0avgbw)*100) %>%  # this is time 0 til before the stress
  mutate(T3_T0_perc=((t3avgbw-t0avgbw)/t0avgbw)*100)      # this is time 0 til after the stress

# calculate relative growth rate (RGR)
calc_phys = calc_phys %>%
  mutate(T0_PA_rgr=((ln(t0avgbw)-ln(paavgbw))/16)) %>%
  mutate(T2_T0_rgr=((ln(t2avgbw)-ln(t0avgbw))/57)) %>%
  mutate(T2_T1_rgr=((ln(t2avgbw)-ln(t1avgbw))/22)) %>%
  mutate(T3_T2_rgr=((ln(t3avgbw)-ln(t2avgbw))/32)) %>%
  mutate(T3_T0_rgr=((ln(t3avgbw)-ln(t0avgbw))/79))

# take a look at the dataset
str(calc_phys)

# un-comment the drop_na() corresponding to the time point of data you want to look at.
calc_phys2 = calc_phys %>%
  #drop_na(T2_T0_rgr) %>%
  drop_na(T3_T2_rgr) %>%
  dplyr::filter(treat!="Control 2") %>%
  dplyr::filter(treat!="Low Var") %>%
  dplyr::filter(treat!="High Var") %>%
  dplyr::filter(gen_site != "I4G") # clone with I4F, remove from dataset

# merge with lineage info for later plotting
lineages = read.csv("Physiology_Data/data_files/tve_lineages_noclones.csv")

calc_phys_all_lin <- left_join(calc_phys2, lineages, by = "gen_site")
calc_phys_all_lin$lineage = as.factor(calc_phys_all_lin$lineage)

# merge with its2 types for plotting
its2_types = read.csv("Physiology_Data/data_files/ITS2.dominanttype.prestress.csv") %>%
  select(frag, dominant_type)

calc_phys_all_lin <- left_join(calc_phys_all_lin, its2_types, by = "frag")
calc_phys_all_lin$dominant_type = as.factor(calc_phys_all_lin$dominant_type)

calc_phys_2_lin = calc_phys_all_lin %>%
  dplyr::filter(is.na(lineage) | lineage!="L3") # want to keep NA values for lineage here since they still have other info, will remove na's for lineage specific plots

# exploratory figure
#scatter plot with linear regression and confidence interval
ggplot(calc_phys, aes(treat, T2_T0_perc, color = sitename))+
  geom_point()+
  geom_smooth(aes(group=sitename), method=lm)+
  theme_classic()


## Mixed Model
# interested in the effects of dtv treatment and lineage
#m1 <- lmer(T2_T0_rgr ~ treat+lineage + (1|gen_site), data = calc_phys_2_lin, REML=TRUE)
m1 <- lmer(T3_T2_rgr ~ treat+lineage + (1|gen_site), data = calc_phys_2_lin, REML=TRUE)
summary(m1)
# T2_T0_rgr:
# Fixed effects:
#                   Estimate Std. Error         df t value Pr(>|t|)
# (Intercept)   2.431e-04  4.225e-05  6.487e+01   5.753 2.58e-07 ***
# treatMod Var  8.780e-05  3.690e-05  4.107e+01   2.379   0.0221 *
# lineageL2    -1.092e-04  6.032e-05  4.546e+01  -1.811   0.0768 .

# T3_T2_rgr:
# Fixed effects:
#               Estimate Std. Error         df t value Pr(>|t|)
# (Intercept)   7.066e-04  6.926e-05  6.762e+01  10.202 2.55e-15 ***
# treatMod Var  7.741e-05  6.930e-05  4.279e+01   1.117 0.270251
# lineageL2    -3.469e-04  9.709e-05  4.603e+01  -3.573 0.000841 ***

anova(m1)
# T2_T0_rgr:
# Type III Analysis of Variance Table with Satterthwaite's method
#             Sum Sq    Mean Sq NumDF  DenDF F value  Pr(>F)
# treat   1.4621e-07 1.4621e-07     1 41.065  5.6615 0.02207 *
# lineage 8.4674e-08 8.4674e-08     1 45.464  3.2786 0.07680 .

# T3_T2_rgr:
# Type III Analysis of Variance Table with Satterthwaite's method
#             Sum Sq    Mean Sq NumDF  DenDF F value    Pr(>F)
# treat   1.1085e-07 1.1085e-07     1 42.789  1.2476 0.2702508
# lineage 1.1343e-06 1.1343e-06     1 46.026 12.7665 0.0008411 ***


# now let's do some more looking into the model
plot_model(m1, "eff", terms="treat")
plot_model(m1, "eff", terms="lineage")
plot_model(m1, type="re")

#Check model fit
r2(m1) #get r2 and adjusted r2
check_model(m1) #check assumptions and model fit

# Now look at custom contrasts with emmmeans

#specify model (because we are interested in pairwise, have to include the interaction)
m.emm<- lmer(T2_T0_rgr ~ treat*lineage + (1|gen_site), data = calc_phys_2_lin, REML=FALSE)
#m.emm<- lmer(T3_T2_rgr ~ treat*lineage + (1|gen_site), data = calc_phys_2_lin, REML=FALSE)

emms<-emmeans(m.emm, ~lineage|treat) #, adjust="Bonferoni"
pairs(emms, interaction = "pairwise") %>% rbind(adjust="fdr")
# for T2_T0_RGR:
# treat   lineage_pairwise estimate       SE   df t.ratio p.value
# Control L1 - L2          8.42e-05 7.32e-05 78.7   1.151  0.2533
# Mod Var L1 - L2          1.30e-04 6.95e-05 74.3   1.873  0.1301

# for T3_T2_RGR:
# treat   lineage_pairwise estimate       SE   df t.ratio p.value
# Control L1 - L2          0.000325 0.000127 80.9   2.560  0.0123
# Mod Var L1 - L2          0.000364 0.000116 77.3   3.125  0.0050


# Make growth figures in manuscript:
#SummarySE to format data for plotting by treatment and lineage

calc_phys_2_lin_nona = calc_phys_2_lin %>%
  drop_na(lineage)

growth_means_2_lin <- summarySE(calc_phys_2_lin_nona, measurevar="T2_T0_rgr", groupvars=c("treat","lineage"))
#growth_means_2_lin <- summarySE(calc_phys_2_lin_nona, measurevar="T3_T2_rgr", groupvars=c("treat","lineage"))

# to calculate percent increase in growth by treatment:
#    treat  N    T2_T0_rgr           sd           se           ci
#1 Control 39 0.0001938333 0.0002046579 3.277149e-05 6.634242e-05
#2 Mod Var 43 0.0002890555 0.0002539979 3.873432e-05 7.816902e-05

# lineage  N    T2_T0_rgr           sd           se           ci
# 1      L1 49 0.0002862905 0.0002483428 3.547754e-05 7.133237e-05
# 2      L2 33 0.0001806258 0.0002020336 3.516953e-05 7.163800e-05

# treat  N    T3_T2_rgr           sd           se           ci
# 1 Control 37 0.0005853828 0.0003648244 5.997677e-05 0.0001216385
# 2 Mod Var 42 0.0006407149 0.0004373524 6.748494e-05 0.0001362886

# lineage  N    T3_T2_rgr           sd           se           ci
# 1      L1 49 0.0007427781 0.0003885958 5.551369e-05 0.0001116178
# 2      L2 30 0.0004057687 0.0003389837 6.188967e-05 0.0001265786

# plot, treatment x axis colored by lineage data figure
calc_plot_lineage <- ggplot(calc_phys_2_lin_nona,aes(x = treat, y = T2_T0_rgr))+
  theme_bw()+
  geom_jitter(aes(color = lineage, fill = lineage),
              position=position_dodge(width=0.3),
              alpha=0.2, pch = 21,
              color = "black") +
  geom_errorbar(data = growth_means_2_lin, aes(x = treat, ymax = T2_T0_rgr+se, ymin = T2_T0_rgr-se, color = lineage), width = .2, position = position_dodge(width=0.4)) +
  geom_point(data = growth_means_2_lin, mapping = aes(x=treat, y=T2_T0_rgr, color = lineage, fill = lineage), size = 3.5, pch = 21, color = "black", position = position_dodge(width=0.4))+
  scale_fill_manual(name = "Lineage",
                    breaks = c("L1","L2"),
                    values = cols_lineage)+
  scale_color_manual(name = "Lineage",
                     breaks = c("L1","L2"),
                     values = cols_lineage)+
  xlab("Treatment")+
  ylab(bquote("Specific growth rate ("~day^-1~')'))+
  ylim(-0.0006,0.0015) +
  geom_hline(yintercept=0, linetype='dotted', color = 'gray')+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))
calc_plot_lineage

ggsave(calc_plot_lineage, filename = "/Users/hannahaichelman/Dropbox/BU/TVE/Growth/Perc_T2T0_lineage_2_lin.pdf", width=5, height=4, units=c("in"), useDingbats=FALSE)


## More exploratory figures below
#SummarySE to format data for plotting by treatment and lineage_symtype
calc_phys_2_lin_nona = calc_phys_2_lin %>%
  drop_na(lineage,dominant_type)

# Plot based on dominant symbiont type
growth_means_lineage_sym <- summarySE(calc_phys_2_lin_nona, measurevar="T2_T0_rgr", groupvars=c("treat","lineage","dominant_type"))
# sample size:
#     treat lineage dominant_type  N     T2_T0_rgr           sd           se           ci
# 1  Control      L1            C1 11  2.698986e-04 0.0002392437 7.213469e-05 0.0001607261
# 2  Control      L1          C3af  1  9.277553e-05           NA           NA           NA
# 3  Control      L1            D1 12  1.930527e-04 0.0001922653 5.550221e-05 0.0001221595
# 4  Control      L2            C1 10  1.680577e-04 0.0001734497 5.484960e-05 0.0001240784
# 5  Control      L2            C3  2 -8.081519e-05 0.0001301622 9.203860e-05 0.0011694612
# 6  Control      L2            D1  2  1.736633e-04 0.0002794062 1.975700e-04 0.0025103651
# 7  Mod Var      L1            C1 12  3.488881e-04 0.0002850910 8.229867e-05 0.0001811382
# 8  Mod Var      L1          C3af  1  3.413739e-04           NA           NA           NA
# 9  Mod Var      L1            D1 11  3.587772e-04 0.0002867804 8.646755e-05 0.0001926617
# 10 Mod Var      L2            C1 14  2.489225e-04 0.0002116028 5.655322e-05 0.0001221758
# 11 Mod Var      L2            C3  1 -7.401003e-05           NA           NA           NA
# 12 Mod Var      L2            D1  2  7.496903e-05 0.0002190365 1.548822e-04 0.0019679651

growth_means_lineage_sym <- summarySE(calc_phys_2_lin_nona, measurevar="T3_T2_rgr", groupvars=c("treat","lineage","dominant_type"))
# sample size:
# treat lineage dominant_type  N     T3_T2_rgr           sd           se           ci
# 1  Control      L1            C1 11  7.705977e-04 3.755353e-04 1.132281e-04 0.0002522880
# 2  Control      L1          C3af  1  7.778082e-04           NA           NA           NA
# 3  Control      L1            D1 12  6.243384e-04 2.151228e-04 6.210062e-05 0.0001366825
# 4  Control      L2            C1  9  4.832587e-04 3.816676e-04 1.272225e-04 0.0002933757
# 5  Control      L2            C3  2 -1.233795e-04 2.451063e-05 1.733163e-05 0.0002202193
# 6  Control      L2            D1  2  4.050748e-04 3.317223e-04 2.345631e-04 0.0029804066
# 7  Mod Var      L1            C1 12  8.335776e-04 3.793931e-04 1.095213e-04 0.0002410548
# 8  Mod Var      L1          C3af  1 -5.833342e-04           NA           NA           NA
# 9  Mod Var      L1            D1 11  8.753668e-04 3.868618e-04 1.166432e-04 0.0002598973
# 10 Mod Var      L2            C1 14  4.928040e-04 2.874217e-04 7.681669e-05 0.0001659524
# 11 Mod Var      L2            C3  1  5.802842e-05           NA           NA           NA
# 12 Mod Var      L2            D1  2  1.515286e-04 1.102259e-04 7.794147e-05 0.0009903403

# plot, treatment x axis colored by site data figure
calc_plot_lineage_sym <- ggplot(calc_phys_2_lin_nona, aes(x = treat, y = T3_T2_rgr))+
  theme_bw()+
  geom_jitter(aes(color = dominant_type, fill = dominant_type, shape = dominant_type),
              position=position_dodge(width=0.3),
              alpha=0.2, pch = 21,
              color = "black") +
  geom_errorbar(data = growth_means_lineage_sym, aes(x = treat, ymax = T3_T2_rgr+se, ymin = T3_T2_rgr-se, color = dominant_type), width = .2, position = position_dodge(width=0.4)) +
  geom_point(data = growth_means_lineage_sym, mapping = aes(x=treat, y=T3_T2_rgr, color = dominant_type, fill = dominant_type, shape = dominant_type), size = 3.5, color = "black", position = position_dodge(width=0.4))+
  scale_fill_manual(name = "Dominant ITS2",
                    #breaks = c("L1","L2"),
                    values = its2_cols_greens)+
  scale_color_manual(name = "Dominant ITS2",
                     #breaks = c("L1","L2"),
                     values = its2_cols_greens)+
  scale_shape_manual(name = "Dominant ITS2",
                     values = c(21,22,23,24))+
  xlab("Treatment")+
  ylab(bquote("Specific growth rate ("~day^-1~')'))+
  ylim(-0.0006,0.0015) +
  geom_hline(yintercept=0, linetype='dotted', color = 'gray')+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  facet_wrap(~lineage)
calc_plot_lineage_sym

ggsave(calc_plot_lineage_sym, filename = "/Users/hannahaichelman/Dropbox/BU/TVE/Growth/calcification_T3T2_lineage_dominantType.pdf", width=8.5, height=4, units=c("in"), useDingbats=FALSE)

# plot just CI corals to see if lineage difference holds
calc_phys_2_lin_nona_CI = calc_phys_2_lin_nona %>%
  subset(sitename == "CI")

growth_means_lineage_CI <- summarySE(calc_phys_2_lin_nona_CI, measurevar="T2_T0_g_cm2_day", groupvars=c("treat","lineage"))

calc_plot_lineage_CI <- ggplot(growth_means_lineage_CI,aes(x = treat, y = T2_T0_g_cm2_day, color = lineage, pch = lineage))+
  theme_bw()+
  geom_point(size = 3, position = position_dodge(width=0.3))+
  geom_smooth(aes(group=lineage), position = position_dodge(width=0.3), size = 0.5, method='lm', se = FALSE)+
  #geom_line(aes(group = sitename), size = 0.5, linetype="dashed", position = position_dodge(width=0.3))+
  geom_errorbar(aes(x = treat, ymax = T2_T0_g_cm2_day+se, ymin = T2_T0_g_cm2_day-se), width = .2, position = position_dodge(width=0.3)) +
  scale_color_manual(name = "Lineage",
                     breaks = c("L1","L2"),
                     values = cols_lineage)+
  scale_shape_manual(name = "Lineage",
                     breaks = c("L1","L2"),
                     values=c(19,17))+
  xlab("Treatment")+
  #ylab("Percent Change in Weight")+
  ylab(bquote("Calcification (g" ~cm^-2~ ~day^-1~')'))+
  #ggtitle("Growth during variability")+
  #ylim(-0.5,4) +
  geom_hline(yintercept=0, linetype='dotted', color = 'gray')+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))
calc_plot_lineage_CI

ggsave(calc_plot_lineage_CI, filename = "/Users/hannahaichelman/Documents/BU/TVE/Growth/calcificationo_T2T0_lineage_CIonly.pdf", width=5, height=4, units=c("in"), useDingbats=FALSE)

#SummarySE to format data for plotting by treatment and site
growth_means_site_all_lin <- summarySE(calc_phys_all_lin, measurevar="T2_T0_perc", groupvars=c("treat","sitename"))
growth_means_site_2_lin <- summarySE(calc_phys_2_lin, measurevar="T2_T0_perc", groupvars=c("treat","sitename"))

# plot, treatment x axis colored by site data figure
calc_plot_site <- ggplot(growth_means_site_2_lin,aes(x = treat, y = T2_T0_perc, color = sitename, pch = sitename))+
  theme_bw()+
  geom_point(size = 3, position = position_dodge(width=0.3))+
  geom_smooth(aes(group=sitename), position = position_dodge(width=0.3), size = 0.5, method='lm', se = FALSE)+
  #geom_line(aes(group = sitename), size = 0.5, linetype="dashed", position = position_dodge(width=0.3))+
  geom_errorbar(aes(x = treat, ymax = T2_T0_perc+se, ymin = T2_T0_perc-se), width = .2, position = position_dodge(width=0.3)) +
  scale_color_manual(name = "Site",
                     labels = c("CI","PD","SP","BN","BS","CA"),
                     values = cols_site)+
  scale_shape_manual(name = "Site",
                     labels = c("CI","PD","SP","BN","BS","CA"),
                     values=c(19,19,19,17,17,17))+
  xlab("Treatment")+
  ylab("Percent Change in Weight")+
  #ylab(bquote("Calcification (g" ~cm^-2~ ~day^-1~')'))+
  #ggtitle("Growth during variability")+
  #ylim(-0.5,4) +
  geom_hline(yintercept=0, linetype='dotted', color = 'gray')+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))
calc_plot_site

ggsave(calc_plot_site, filename = "/Users/hannahaichelman/Documents/BU/TVE/Growth/calcification_T2T0_site_2_lin.pdf", width=5, height=4, units=c("in"), useDingbats=FALSE)

#plot, faceted by treatment to looking for outliers
plot <- ggplot(calc_phys,aes(x = sitename, y = T2_T0_perc, color = sitename, pch = sitename, text = paste("Nubbin:",frag)))+
  theme_bw()+
  geom_boxplot()+
  ggtitle("Percent Change in Weight T0-T2")+
  scale_fill_manual(values = palsite)+
  geom_jitter(aes(color = sitename),width = 0.3)+
  scale_color_manual(values = palsite)+
  theme(axis.title.x = element_blank())+
  theme(axis.text.x=element_text(angle = 45, vjust = 0.5))+
  facet_wrap(~treat, nrow=5)
#  labs(y=expression(paste("Protein (mg cm"^-2*' )')))
ggplotly(subplot(list(plot),titleY=F) %>% layout(showlegend=F))


# manipulate the data to create a data frame we can use to look at the time course of percent change in growth
calc_phys_timecourse = calc_phys %>%
  select(frag, treat, gen_site, sitename, reef,
         T0_PA_perc,T1_T0_perc,T2_T1_perc,T3_T2_perc)

head(calc_phys_timecourse)

calc_phys_timecourse.melt <- melt(calc_phys_timecourse, id.vars = c('frag','treat','gen_site','reef','sitename'), measure.vars = c('T0_PA_perc','T1_T0_perc','T2_T1_perc','T3_T2_perc'))
head(calc_phys_timecourse.melt)
str(calc_phys_timecourse.melt)

colnames(calc_phys_timecourse.melt)[6] <- "time"
levels(calc_phys_timecourse.melt$time)
levels(calc_phys_timecourse.melt$time) <- c("PA-T0","T0-T1","T1-T2","T2-T3")

# Stats for time course changes:
# first try with lmer()
phys_calc_time_nona = calc_phys_timecourse.melt %>%
  drop_na(value) %>%
  filter(treat!="Control 2")

model.calc.time <- lmer(value ~ treat + sitename + time + (1|gen_site), data = phys_calc_time_nona)
summary(model.calc.time)

# SummarySE and plot for time course
growth_time_means <- summarySE(phys_calc_time_nona, measurevar="value", groupvars=c("treat","sitename","time"))

# plot, treatment x axis colored by site data figure
plot.time <- ggplot(growth_time_means,aes(x = time, y = value, color = sitename, pch = treat))+
  theme_bw()+
  geom_point(size = 3, position = position_dodge(width=0.3))+
  #geom_smooth(aes(group=sitename), position = position_dodge(width=0.3), size = 0.5, method='lm', se = FALSE)+
  #geom_line(aes(group = sitename), size = 0.5, linetype="dashed", position = position_dodge(width=0.3))+
  geom_errorbar(aes(x = time, ymax = value+se, ymin = value-se), width = .2, position = position_dodge(width=0.3)) +
  scale_color_manual(name = "Site",
                     labels = c("I-Cristobal","I-Punta Donato","I-STRI Point","O-Bastimentos N","O-Bastimentos S","O-Cayo de Agua"),
                     values = cols_site)+
  scale_shape_manual(name = "Treatment",
                     labels = c("Control 1","Low Var","Mod Var","High Var"),
                     values=c(15,16,17,18))+
  xlab("Time")+
  ylab("Percent Change in Weight")+
  #ggtitle("Growth during variability")+
  #ylim(-0.5,4) +
  geom_hline(yintercept=0, linetype='dotted', color = 'gray')+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))
plot.time
ggsave(file="/Users/hannahaichelman/Documents/BU/TVE/Growth/PercentChange_timecourse.pdf", plot.time, width=8, height=5, units=c("in"), useDingbats=FALSE)

# plot, treatment x axis colored by site data figure
plot.time.2 <- ggplot(growth_time_means,aes(x = time, y = value, color = treat, pch = sitename))+
  theme_bw()+
  geom_point(size = 3, position = position_dodge(width=0.3))+
  geom_errorbar(aes(x = time, ymax = value+se, ymin = value-se), width = .2, position = position_dodge(width=0.3)) +
  scale_color_manual(name = "Treatment",
                     labels = c("Control","Low Var","Mod Var","High Var"),
                     values = cols_treat)+
  scale_shape_manual(name = "Site",
                     labels = c("CI","PD","SP","BN","BS","CA"),
                     values=c(15,16,17,18,4,8))+
  xlab("Time")+
  ylab("Percent Change in Weight")+
  #ggtitle("Growth during variability")+
  #ylim(-0.5,4) +
  geom_hline(yintercept=0, linetype='dotted', color = 'gray')+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  facet_wrap(~sitename, nrow=2, ncol=3)
plot.time.2
ggsave(file="/Users/hannahaichelman/Documents/BU/TVE/Growth/PercentChange_timecourse_treatfacet.pdf", plot.time.2, width=8, height=5, units=c("in"), useDingbats=FALSE)


## Plot 90% daily temperature range by % change in growth
# an estimate of the magnitude of high-frequency temp fluctuations

# logger          percent90
#  Punta.IR1 (PD)           0.999
#  STRI.IR2  (SP)           0.833
#  Cristo.IR3.arr1 (CI)     1.92
#  Cayo.OR3.arr3 (CA)      1.49
#  Drago.OR4 (DM)          1.24

# add in column for the 90% DTR
phys_calc_nona_lineage$dtr_90 <- ifelse(phys_calc_nona_lineage$sitename == 'PD', '0.999',
                                        ifelse(phys_calc_nona_lineage$sitename == 'SP', '0.833',
                                               ifelse(phys_calc_nona_lineage$sitename == 'CI', '1.92',
                                                      ifelse(phys_calc_nona_lineage$sitename == 'CA', '1.49',
                                                             'NA'))))
phys_calc_nona_lineage$dtr_90 <- as.numeric(phys_calc_nona_lineage$dtr_90)
phys_calc_nona_lineage_dtr = phys_calc_nona_lineage %>%
  drop_na(dtr_90)

#SummarySE to format data for plotting by treatment and site
growth_means_dtr <- summarySE(phys_calc_nona_lineage_dtr, measurevar="T2_T0_perc", groupvars=c("dtr_90","treat"))

# plot, treatment x axis colored by site data figure

calc_plot_dtr <- ggplot(growth_means_dtr,aes(x = dtr_90, y = T2_T0_perc, color = treat, pch = treat))+
  theme_bw()+
  geom_point(size = 3, position = position_dodge(width=0.1))+
  geom_smooth(aes(group=treat), position = position_dodge(width=0.1), size = 0.5, method='lm', se = FALSE)+
  #geom_line(aes(group = sitename), size = 0.5, linetype="dashed", position = position_dodge(width=0.3))+
  geom_errorbar(aes(x = dtr_90, ymax = T2_T0_perc+se, ymin = T2_T0_perc-se), width = .2, position = position_dodge(width=0.1)) +
  scale_color_manual(name = "Treatment",
                     labels = c("Control","Low Var","Mod Var","High Var"),
                     values = cols_treat)+
  scale_shape_manual(name = "Treatment",
                     labels = c("Control","Low Var","Mod Var","High Var"),
                     values=c(15,16,17,18))+
  xlab("90% Daily Temperature Range (°C)")+
  ylab("Percent Change in Weight")+
  #ggtitle("Growth during variability")+
  #ylim(-0.5,4) +
  geom_hline(yintercept=0, linetype='dotted', color = 'gray')+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))
calc_plot_dtr

ggsave(calc_plot_dtr, filename = "/Users/hannahaichelman/Documents/BU/TVE/Growth/percentchange_T2T0_dtr.pdf", width=5, height=4, units=c("in"), useDingbats=FALSE)

# seems to not be an effect of reef zone so let's check that out
library(jtools)
m.dtr <- lmer(T2_T0_perc ~ dtr_90 + treat + (1|gen_site), data = phys_calc_nona_lineage_dtr, REML=TRUE)
summ(m.dtv)


##### PAM (Fv/Fm) #####
# need to re-read in the data sheet for pam only until I can organize better
post_phys_forpam <- read.csv('Physiology_Data/data_files/dtvmaster.csv')

# analyzing separately because not sure of timing yet
pam_phys = post_phys_forpam %>%
  dplyr::select(frag, survivedtoend, treat, t3sastan1, t3sastan2, t3sastan3, t3sarec1, t3sarec2, t3sarec3, blastvol, blaster,
                papamrec1,papamrec2,papamrec3,t0pamrec1,t0pamrec2,t0pamrec3,t1pamrec1,t1pamrec2,t1pamrec3,t2pamrec1,t2pamrec2,t2pamrec3,t3pamrec1,t3pamrec2,t3pamrec3,
                t4pamrec1,t4pamrec2,t4pamrec3,t5pamrec1,t5pamrec2,t5pamrec3,t6pamrec1,t6pamrec2,t6pamrec3,t7pamrec1,t7pamrec2,t7pamrec3,t8pamrec1,t8pamrec2,t8pamrec3,t9pamrec1,t9pamrec2,t9pamrec3) %>%
  dplyr::rename(sastan1 = t3sastan1, sastan2 = t3sastan2, sastan3 = t3sastan3, sarec1 = t3sarec1, sarec2 = t3sarec2, sarec3 = t3sarec3) %>%
  mutate(treat = as.factor(treat), sarec3 = as.numeric(sarec3), blastvol = as.numeric(blastvol)) %>%
  dplyr::filter(frag!="I2A12", frag!="I2A7", frag!="O4D8", frag!="I3E10", frag!="I3I10") #these frags are being removed because they were duplicate genotypes within treatment and had the most complete information of the two

# add identifying data
pam_phys$origsitecode <- substr(pam_phys$frag, 1, 2)

# add in site name with inshore/offshore indicator
pam_phys$sitename <- ifelse(pam_phys$origsitecode == 'I2', 'SP',
                            ifelse(pam_phys$origsitecode == 'I3', 'CI',
                                   ifelse(pam_phys$origsitecode == 'I4', 'PD',
                                          ifelse(pam_phys$origsitecode == 'O2', 'BS',
                                                 ifelse(pam_phys$origsitecode == 'O3', 'CA',
                                                        'BN')))))
pam_phys$sitename <- as.factor(pam_phys$sitename)

# make new nubbin IDs based on the new, more informative site codes
# extract the genotype and frag number
pam_phys$fragid <- substr(pam_phys$frag,3,5)

# add inshore/offshore designation
pam_phys$reef <- substr(pam_phys$frag,1,1)
pam_phys$reef <- ifelse(pam_phys$reef == 'O', 'Outer Reef', 'Inner Reef')
pam_phys$reef <- as.factor(pam_phys$reef)
pam_phys$genet <- substr(pam_phys$fragid,1,1)

#create a new column of combined genotype and site for stats later
pam_phys = pam_phys %>%
  unite(gen_site, c(origsitecode,genet), sep = "", remove = FALSE) %>%
  mutate(gen_site = as.factor(gen_site))

# re-level and re-name treatment
pam_phys$treat <- factor(pam_phys$treat, levels = c("init", "1", "2", "3","4","5"))
levels(pam_phys$treat) <- c("Initial","Control","Low Var","Mod Var","High Var","Control 2")
pam_phys$treat <- as.factor(pam_phys$treat)

# calculate average fv/fm for each time point
# PA PAM on 9.2.2016 & 9.3.2016 (day -18)
# T0 PAM on 9.20.2016 (day 0)
# T1 PAM on 10.5.2016 (day 15)
# T2 PAM on 10.25.2016 (day 35)
# T3 PAM on 11.5.2016 (taken 4 days before start of heat stress ramp) (day 45)
# T4 PAM on 11.14.2016 (first day corals were at 32°C) (day 54) - temp increase started day 51
# T5 PAM on 11.21.2016 (end of 32 hold, which lasted 8 days) (day 61)
# T6 PAM on 11.25.2016 (after cool down ramp, back at 29.5°C) (day 65)
# T7 PAM on 11.30.2016 (on day 5 of recovery at 29.5°C) (day 70)
# T8 PAM on 12.4.2016 (on day 9 of recovery at 29.5°C) (day 74)
# T9 PAM on 12.9.2016 (on day 14 of recovery at 29.5°C) (day 79)
# the experiment ended on 12.13.2016

# calculate average PAM at each time point
pam_phys = pam_phys %>%
  mutate(paavgpam=rowMeans(.[ , c("papamrec1","papamrec2","papamrec3")], na.rm=TRUE)) %>%
  mutate(t0avgpam=rowMeans(.[ , c("t0pamrec1","t0pamrec2","t0pamrec3")], na.rm=TRUE)) %>%
  mutate(t1avgpam=rowMeans(.[ , c("t1pamrec1","t1pamrec2","t1pamrec3")], na.rm=TRUE)) %>%
  mutate(t2avgpam=rowMeans(.[ , c("t2pamrec1","t2pamrec2","t2pamrec3")], na.rm=TRUE)) %>%
  mutate(t3avgpam=rowMeans(.[ , c("t3pamrec1","t3pamrec2","t3pamrec3")], na.rm=TRUE)) %>%
  mutate(t4avgpam=rowMeans(.[ , c("t4pamrec1","t4pamrec2","t4pamrec3")], na.rm=TRUE)) %>%
  mutate(t5avgpam=rowMeans(.[ , c("t5pamrec1","t5pamrec2","t5pamrec3")], na.rm=TRUE)) %>%
  mutate(t6avgpam=rowMeans(.[ , c("t6pamrec1","t6pamrec2","t6pamrec3")], na.rm=TRUE)) %>%
  mutate(t7avgpam=rowMeans(.[ , c("t7pamrec1","t7pamrec2","t7pamrec3")], na.rm=TRUE)) %>%
  mutate(t8avgpam=rowMeans(.[ , c("t8pamrec1","t8pamrec2","t8pamrec3")], na.rm=TRUE)) %>%
  mutate(t9avgpam=rowMeans(.[ , c("t9pamrec1","t9pamrec2","t9pamrec3")], na.rm=TRUE)) %>%
  mutate(pamdiff=t9avgpam-t3avgpam)

# explore difference in PAM values from the end of recovery period to start of the heat stress ramp
str(pam_phys)
plot(x=pam_phys$sitename,y=pam_phys$pamdiff, data = pam_phys)

lm02 = aov(pamdiff ~ sitename + treat, data = pam_phys)
summary(lm02)
par(mfrow=c(2,2))
plot(lm02)
TukeyHSD(lm02)
# not a lot interesting happening here, on to stats and plotting

# summarySE doesn't work with NA's, so use this!!
phys_pam_wide = pam_phys %>%
  dplyr::filter(complete.cases(paavgpam,t0avgpam,t1avgpam,t2avgpam,t3avgpam,t4avgpam,t5avgpam,t6avgpam,t7avgpam,t8avgpam,t9avgpam)) %>% #drop any row that has an NA for any time point
  dplyr::filter(treat!="Control 2") %>%
  dplyr::filter(treat!="Low Var") %>%
  dplyr::filter(treat!="High Var") %>%
  select(frag,treat,reef,gen_site,sitename,paavgpam,t0avgpam,t1avgpam,t2avgpam,t3avgpam,t4avgpam,t5avgpam,t6avgpam,t7avgpam,t8avgpam,t9avgpam) %>%
  dplyr::filter(gen_site != "I4G") # remove clone

# transform the data to long format so time point is its own column
phys_pam_long = phys_pam_wide %>%
  gather(time, pam, paavgpam:t9avgpam)

# re-level and re-name treatment
phys_pam_long$time <- as.factor(phys_pam_long$time)
levels(phys_pam_long$time) <- c("-18","0", "15","35","45","54","61","65","70","74","79")

# merge with lineage info for later plotting
lineages = read.csv("Physiology_Data/data_files/tve_lineages_noclones.csv")

pam_phys_all_lin <- left_join(phys_pam_long, lineages, by = "gen_site")
pam_phys_all_lin$lineage = as.factor(pam_phys_all_lin$lineage)

# merge with its2 types for plotting
its2_types = read.csv("Physiology_Data/data_files/ITS2.dominanttype.prestress.csv") %>%
  select(frag, dominant_type)

pam_phys_all_lin <- left_join(pam_phys_all_lin, its2_types, by = "frag")
pam_phys_all_lin$dominant_type = as.factor(pam_phys_all_lin$dominant_type)

# filter to only include T3-T9 for plotting and stats
phys_pam_all_lin_plots = pam_phys_all_lin %>%
  dplyr::filter(time != "-18" & time != "0" & time != "15" & time != "35")

phys_pam_2_lin_plots = phys_pam_all_lin_plots %>%
  dplyr::filter(is.na(lineage) | lineage!="L3") # want to keep NA values for lineage here since they still have other info, will remove na's for lineage specific plots

phys_pam_lin1 = phys_pam_2_lin_plots %>%
  dplyr::filter(lineage == "L1")

phys_pam_lin2 = phys_pam_2_lin_plots %>%
  dplyr::filter(lineage == "L2")

# Stats
# Using the full model with interactions of all our parameters of interest
m.full <- lmer(pam ~ time*lineage*treat + (1|gen_site), data = phys_pam_2_lin_plots, REML=TRUE)
summary(m.full)
anova(m.full)
# Type III Analysis of Variance Table with Satterthwaite's method
#                      Sum Sq   Mean Sq NumDF  DenDF F value    Pr(>F)
# time               0.034591 0.0057652     6 305.29  4.1089 0.0005536 ***
# lineage            0.009568 0.0095677     1  38.09  6.8190 0.0128313 *
# treat              0.002252 0.0022521     1 335.72  1.6051 0.2060558
# time:lineage       0.014009 0.0023348     6 305.29  1.6641 0.1294167
# time:treat         0.005664 0.0009440     6 305.29  0.6728 0.6717092
# lineage:treat      0.002250 0.0022496     1 335.72  1.6033 0.2063112
# time:lineage:treat 0.008302 0.0013836     6 305.29  0.9861 0.4346687

# check interactions
m.emm<- lmer(pam ~ time*lineage + (1|gen_site), data = phys_pam_2_lin_plots, REML=FALSE)

emms<-emmeans(m.emm, ~lineage|time) #, adjust="Bonferoni"
pairs(emms, interaction = "pairwise") %>% rbind(adjust="fdr")
# time lineage_pairwise estimate     SE   df t.ratio p.value
# 45   L1 - L2            0.0254 0.0197 72.9   1.289  0.2016
# 54   L1 - L2            0.0451 0.0197 72.9   2.290  0.0581
# 61   L1 - L2            0.0610 0.0197 72.9   3.098  0.0107
# 65   L1 - L2            0.0603 0.0197 72.9   3.063  0.0107
# 70   L1 - L2            0.0402 0.0197 72.9   2.043  0.0782
# 74   L1 - L2            0.0342 0.0197 72.9   1.736  0.1013
# 79   L1 - L2            0.0361 0.0197 72.9   1.835  0.0988

m.emm2<- lmer(pam ~ time*treat + (1|gen_site), data = phys_pam_2_lin_plots, REML=FALSE)
emms2<-emmeans(m.emm2, ~treat) #, adjust="Bonferoni"
pairs(emms2, interaction = "pairwise") %>% rbind(adjust="fdr")
# treat_pairwise    estimate     SE  df t.ratio p.value
# Control - Mod Var  -0.0075 0.0052 397  -1.441  0.1504

m.emm3<- lmer(pam ~ dominant_type*lineage + (1|gen_site), data = phys_pam_2_lin_plots, REML=FALSE)
m.emm3<-emmeans(m.emm3, ~dominant_type|lineage) #, adjust="Bonferoni"
pairs(m.emm3, interaction = "pairwise") %>% rbind(adjust="fdr")


#Check model fit
r2(m.full) #get r2 and adjusted r2
check_model(m.full) #check assumptions and model fit

# Now plot Fv/Fm data
#SummarySE to format data for plotting - lineage
phys_pam_all_lin_plots_nona = phys_pam_all_lin_plots %>%
  drop_na(lineage)

phys_pam_2_lin_plots_nona = phys_pam_2_lin_plots %>%
  drop_na(lineage)

pam_means_all_lin <- summarySE(phys_pam_all_lin_plots_nona, measurevar="pam", groupvars=c("lineage","time"))
pam_means_all_lin$time <- as.numeric(as.character(pam_means_all_lin$time))
pam_means_2_lin <- summarySE(phys_pam_2_lin_plots_nona, measurevar="pam", groupvars=c("lineage","time"))
pam_means_2_lin$time <- as.numeric(as.character(pam_means_2_lin$time))

# lineage   N       pam         sd          se          ci
# 1      L1 224 0.5492359 0.05034007 0.003363488 0.006628287
# 2      L2 147 0.5061236 0.07849098 0.006473827 0.012794519

# plot, treatment x axis colored by site data figure
pam_plot_lineage <- ggplot(pam_means_2_lin,aes(x = time, y = pam, color = lineage, fill = lineage))+
  theme_bw()+
  annotate("rect", xmin = 51, xmax = 64, ymin = - Inf, ymax = Inf, fill = "red4", alpha = 0.15)+
  annotate("rect", xmin = 64, xmax = Inf, ymin = - Inf, ymax = Inf, fill = "royalblue4", alpha = 0.15)+
  geom_errorbar(aes(x = time, ymax = pam+se, ymin = pam-se), width = .2, color = "black", position = position_dodge(width=0.3)) +
  geom_line(aes(group = lineage), size = 1, linetype="dashed", position = position_dodge(width=0.3))+
  geom_point(size = 3.5, color = "black", pch = 21, position = position_dodge(width=0.3))+
  #geom_smooth(aes(group=sitename), position = position_dodge(width=0.3), size = 0.5, method='lm', se = FALSE)+
  scale_fill_manual(name = "Lineage",
                    breaks = c("L1","L2"),
                    #breaks = c("L1","L2","L3"),
                    values = cols_lineage)+
  scale_color_manual(name = "Lineage",
                     breaks = c("L1","L2"),
                     #breaks = c("L1","L2","L3"),
                     values = cols_lineage)+
  ylab("Fv/Fm")+
  scale_x_continuous(name = "Time Point", breaks = c(45,54,61,65,70,74,79)) +
  #ylim(-0.5,4) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))
pam_plot_lineage

ggsave(pam_plot_lineage, filename = "/Users/hannahaichelman/Dropbox/BU/TVE/PAM/pam_lineage_2.pdf", width=5, height=4, units=c("in"), useDingbats=FALSE)

# plot just CI corals to see if lineage difference holds
phys_pam_lineage_plots_CI = phys_pam_2_lin_plots_nona %>%
  subset(sitename == "CI")

pam_lineagemeans_CI = summarySE(phys_pam_lineage_plots_CI, measurevar="pam", groupvars=c("lineage","time"))
pam_lineagemeans_CI$time = as.numeric(as.character(pam_lineagemeans_CI$time))

# plot, time x axis, shape = symbiont type, color = lineage
pam_plot_lineage_CI <- ggplot(pam_lineagemeans_CI,aes(x = time, y = pam, color = lineage, fill = lineage))+
  theme_bw()+
  annotate("rect", xmin = 51, xmax = 64, ymin = - Inf, ymax = Inf, fill = "red4", alpha = 0.15)+
  annotate("rect", xmin = 64, xmax = Inf, ymin = - Inf, ymax = Inf, fill = "royalblue4", alpha = 0.15)+
  geom_point(size = 3, color = "black", pch = 21, position = position_dodge(width=0.3))+
  geom_line(aes(color = lineage), size = 0.5, linetype="dashed", position = position_dodge(width=0.3))+
  geom_errorbar(aes(x = time, ymax = pam+se, ymin = pam-se), width = .2, position = position_dodge(width=0.3)) +
  scale_fill_manual(name = "Lineage",
                    breaks = c("L1","L2"),
                    values = cols_lineage)+
  scale_color_manual(name = "Lineage",
                     breaks = c("L1","L2"),
                     values = cols_lineage)+
  ylab("Fv/Fm")+
  scale_x_continuous(name = "Time Point", breaks = c(45,54,61,65,70,74,79)) +
  #ylim(-0.5,4) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))
pam_plot_lineage_CI

ggsave(pam_plot_lineage_CI, filename = "/Users/hannahaichelman/Dropbox/BU/TVE/PAM/pam_lineage_CIonly.pdf", width=6, height=4, units=c("in"), useDingbats=FALSE)

#SummarySE to format data for plotting - dtv treatment
pam_means_treat_all_lin <- summarySE(phys_pam_all_lin_plots, measurevar="pam", groupvars=c("treat","time"))
pam_means_treat_all_lin$time <- as.numeric(as.character(pam_means_treat_all_lin$time))
pam_means_treat_2_lin <- summarySE(phys_pam_2_lin_plots, measurevar="pam", groupvars=c("treat","time"))
pam_means_treat_2_lin$time <- as.numeric(as.character(pam_means_treat_2_lin$time))

pam_plot_treatment <- ggplot(pam_means_treat_2_lin,aes(x = time, y = pam, color = treat, fill = treat))+
  theme_bw()+
  annotate("rect", xmin = 51, xmax = 64, ymin = - Inf, ymax = Inf, fill = "red4", alpha = 0.15)+
  annotate("rect", xmin = 64, xmax = 80, ymin = - Inf, ymax = Inf, fill = "royalblue4", alpha = 0.15)+
  geom_errorbar(aes(x = time, ymax = pam+se, ymin = pam-se), width = .2, position = position_dodge(width=0.3)) +
  geom_point(size = 3, pch = 21, color = "black", position = position_dodge(width=0.3))+
  geom_line(aes(group = treat), size = 0.5, linetype="dashed", position = position_dodge(width=0.3))+
  scale_fill_manual(name = "Treatment",
                    labels = c("Control","Mod Var"),
                    values = cols_treat)+
  scale_color_manual(name = "Treatment",
                     labels = c("Control","Mod Var"),
                     values = cols_treat)+
  ylab("Fv/Fm")+
  scale_x_continuous(name = "Time Point", breaks = c(45,54,61,65,70,74,79))
#ylim(0.4,0.7) +
pam_plot_treatment

ggsave(pam_plot_treatment, file="/Users/hannahaichelman/Dropbox/BU/TVE/PAM/PAM_treat_2lin.pdf", width=8, height=5, units=c("in"), useDingbats=FALSE)

#SummarySE to format data for plotting - dtv treatment faceted by lineage
pam_means_treat_all_lin <- summarySE(phys_pam_all_lin_plots, measurevar="pam", groupvars=c("treat","time","lineage"))
pam_means_treat_all_lin$time <- as.numeric(as.character(pam_means_treat_all_lin$time))
pam_means_treat_2_lin <- summarySE(phys_pam_2_lin_plots, measurevar="pam", groupvars=c("treat","time","lineage"))
pam_means_treat_2_lin$time <- as.numeric(as.character(pam_means_treat_2_lin$time))

pam_means_treat_2_lin = pam_means_treat_2_lin %>%
  drop_na(lineage)

pam_plot_treatment_lin <- ggplot(pam_means_treat_2_lin,aes(x = time, y = pam, color = treat, fill = treat))+
  theme_bw()+
  annotate("rect", xmin = 51, xmax = 64, ymin = - Inf, ymax = Inf, fill = "red4", alpha = 0.15)+
  annotate("rect", xmin = 64, xmax = 80, ymin = - Inf, ymax = Inf, fill = "royalblue4", alpha = 0.15)+
  geom_errorbar(aes(x = time, ymax = pam+se, ymin = pam-se), width = .2, position = position_dodge(width=0.3)) +
  geom_point(size = 3, pch = 21, color = "black", position = position_dodge(width=0.3))+
  geom_line(aes(group = treat), size = 0.5, linetype="dashed", position = position_dodge(width=0.3))+
  scale_fill_manual(name = "Treatment",
                    labels = c("Control","Mod Var"),
                    values = cols_treat)+
  scale_color_manual(name = "Treatment",
                     labels = c("Control","Mod Var"),
                     values = cols_treat)+
  ylab("Fv/Fm")+
  scale_x_continuous(name = "Time Point", breaks = c(45,54,61,65,70,74,79))+
  facet_wrap(~lineage)
#ylim(0.4,0.7) +
pam_plot_treatment_lin

ggsave(pam_plot_treatment_lin, file="/Users/hannahaichelman/Dropbox/BU/TVE/PAM/PAM_treat_lin_facet.pdf", width=6, height=4, units=c("in"), useDingbats=FALSE)


#SummarySE to format data for plotting - dominant symtype
phys_pam_2_linsym_nona = phys_pam_2_lin_plots %>%
  drop_na(lineage, dominant_type)

phys_pam_all_linsym_nona = phys_pam_all_lin_plots %>%
  drop_na(lineage, dominant_type)

pam_2_lineagesymtype_means <- summarySE(phys_pam_2_linsym_nona, measurevar="pam", groupvars=c("time","dominant_type","lineage"))
pam_2_lineagesymtype_means$time <- as.numeric(as.character(pam_2_lineagesymtype_means$time))

pam_all_lineagesymtype_means <- summarySE(phys_pam_all_linsym_nona, measurevar="pam", groupvars=c("time","dominant_type","lineage"))
pam_all_lineagesymtype_means$time <- as.numeric(as.character(pam_all_lineagesymtype_means$time))

its2_cols_greens = c("#edf8e9","#a1d99b",  "#238b45", "#00441b")

# plot, time x axis colored by dominant sym type data figure
pam_plot_lineage_symtype = pam_2_lineagesymtype_means %>%
  #filter(lineage == "L2") %>%
  ggplot(aes(x = time, y = pam, fill = dominant_type, color = dominant_type, shape = dominant_type))+
  theme_bw()+
  annotate("rect", xmin = 51, xmax = 64, ymin = - Inf, ymax = Inf, fill = "red4", alpha = 0.15)+
  annotate("rect", xmin = 64, xmax = Inf, ymin = - Inf, ymax = Inf, fill = "royalblue4", alpha = 0.15)+
  geom_errorbar(aes(x = time, ymax = pam+se, ymin = pam-se), width = .2, color = "black", position = position_dodge(width=1.5)) +
  geom_line(aes(group = dominant_type), size = 1, linetype="dashed", position = position_dodge(width=1.5))+
  geom_point(size = 3.5, color = "black", position = position_dodge(width=1.5))+
  scale_color_manual(name = "Dominant ITS2",
                     values = its2_cols_greens)+
  scale_fill_manual(name = "Dominant ITS2",
                    values = its2_cols_greens)+
  scale_shape_manual(name = "Dominant ITS2",
                     values = c(21,22,23,24))+
  ylab("Fv/Fm")+
  scale_x_continuous(name = "Time Point", breaks = c(45,54,61,65,70,74,79)) +
  #ylim(-0.5,4) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  facet_wrap(~lineage)
pam_plot_lineage_symtype

ggsave(pam_plot_lineage_symtype, filename = "/Users/hannahaichelman/Dropbox/BU/TVE/PAM/pam_lineage_sym_2_lin_facet_dominanttype.pdf", width=8.5, height=4, units=c("in"), useDingbats=FALSE)

#SummarySE to format data for plotting - CI only dominant symtype
phys_pam_2_linsym_CI = phys_pam_2_linsym_nona %>%
  filter(sitename == "CI")

phys_pam_2_linsym_CI_means <- summarySE(phys_pam_2_linsym_CI, measurevar="pam", groupvars=c("time","dominant_type","lineage"))
phys_pam_2_linsym_CI_means$time <- as.numeric(as.character(phys_pam_2_linsym_CI_means$time))

its2_cols_greens_CI = c("#edf8e9", "#00441b")

# plot, time x axis colored by dominant sym type data figure
pam_plot_lineage_symtype = phys_pam_2_linsym_CI_means %>%
  #filter(lineage == "L2") %>%
  ggplot(aes(x = time, y = pam, fill = dominant_type, color = dominant_type, shape = dominant_type))+
  theme_bw()+
  annotate("rect", xmin = 51, xmax = 64, ymin = - Inf, ymax = Inf, fill = "red4", alpha = 0.15)+
  annotate("rect", xmin = 64, xmax = Inf, ymin = - Inf, ymax = Inf, fill = "royalblue4", alpha = 0.15)+
  geom_errorbar(aes(x = time, ymax = pam+se, ymin = pam-se), width = .2, color = "black", position = position_dodge(width=1.5)) +
  geom_line(aes(group = dominant_type), size = 1, linetype="dashed", position = position_dodge(width=1.5))+
  geom_point(size = 3.5, color = "black", position = position_dodge(width=1.5))+
  scale_color_manual(name = "Dominant ITS2",
                     values = its2_cols_greens_CI)+
  scale_fill_manual(name = "Dominant ITS2",
                    values = its2_cols_greens_CI)+
  scale_shape_manual(name = "Dominant ITS2",
                     values = c(21,24))+
  ylab("Fv/Fm")+
  scale_x_continuous(name = "Time Point", breaks = c(45,54,61,65,70,74,79)) +
  #ylim(-0.5,4) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  facet_wrap(~lineage)
pam_plot_lineage_symtype

ggsave(pam_plot_lineage_symtype, filename = "/Users/hannahaichelman/Documents/BU/TVE/PAM/pam_lineage_sym_2_lin_facet_dominanttype_CIonly.pdf", width=8.5, height=4, units=c("in"), useDingbats=FALSE)


#SummarySE to format data for plotting - combine lineage and dominant sym type
phys_pam_2_lin_plots2 = phys_pam_2_lin_plots %>%
  drop_na(lineage, dominant_type) %>%
  unite(lin_sym, c(lineage,dominant_type), sep = "_", remove = FALSE) %>%
  mutate(lin_sym = as.factor(lin_sym))


pam_2_lineagediv_means <- summarySE(phys_pam_2_lin_plots2, measurevar="pam", groupvars=c("time","lin_sym","sitename"))
pam_2_lineagediv_means$time <- as.numeric(as.character(pam_2_lineagediv_means$time))

# plot, treatment x axis colored by site data figure
pam_plot_lineage_symtype <- ggplot(pam_2_lineagediv_means,aes(x = time, y = pam, color = lin_sym))+
  theme_bw()+
  annotate("rect", xmin = 46, xmax = 63, ymin = - Inf, ymax = Inf, fill = "red4", alpha = 0.15)+
  annotate("rect", xmin = 63, xmax = 80, ymin = - Inf, ymax = Inf, fill = "royalblue4", alpha = 0.15)+
  geom_point(size = 3, position = position_dodge(width=0.3))+
  geom_line(aes(group = lin_sym), size = 0.5, linetype="dashed", position = position_dodge(width=0.3))+
  geom_errorbar(aes(x = time, ymax = pam+se, ymin = pam-se), width = .2, position = position_dodge(width=0.3)) +
  #scale_color_manual(name = "Lineage-Sym",
  #                   values = c("#D53E4F", "#F46D43", "#3288BD", "#FEE08B"))+
  #values = c("#3f007d","#3f007d","#807dba","#807dba","#bcbddc","#bcbddc"))+
  #scale_fill_manual(name = "Lineage-Sym",
  #                  values = c("#D53E4F", "#F46D43", "#3288BD", "#FEE08B"))+
  #values = c("#3f007d","#3f007d","#807dba","#807dba","#bcbddc","#bcbddc"))+
  #scale_shape_manual(name = "Lineage-Sym",
  #                   values=c(21,24,21,24,21,24))+
  ylab("Fv/Fm")+
  scale_x_continuous(name = "Time Point", breaks = c(46,55,62,66,71,75,80)) +
  #ylim(-0.5,4) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  facet_wrap(~sitename)
pam_plot_lineage_symtype

ggsave(pam_plot_lineage_symtype, filename = "/Users/hannahaichelman/Documents/BU/TVE/PAM/pam_lineage_sym_2_lin_DIVs.pdf", width=8, height=4, units=c("in"), useDingbats=FALSE)

#SummarySE to format data for plotting
pam_means_site_all_lin <- summarySE(phys_pam_all_lin_plots, measurevar="pam", groupvars=c("sitename","time"))
pam_means_site_all_lin$time <- as.numeric(as.character(pam_means_site_all_lin$time))
pam_means_site_2_lin <- summarySE(phys_pam_2_lin_plots, measurevar="pam", groupvars=c("sitename","time"))
pam_means_site_2_lin$time <- as.numeric(as.character(pam_means_site_2_lin$time))

# plot, treatment x axis colored by site data figure
pam_plot_site <- ggplot(pam_means_site_2_lin,aes(x = time, y = pam, color = sitename, fill = sitename))+
  theme_bw()+
  annotate("rect", xmin = 51, xmax = 64, ymin = - Inf, ymax = Inf, fill = "red4", alpha = 0.15)+
  annotate("rect", xmin = 64, xmax = 80, ymin = - Inf, ymax = Inf, fill = "royalblue4", alpha = 0.15)+
  geom_point(size = 3, pch = 21, color = "black", position = position_dodge(width=0.3))+
  #geom_smooth(aes(group=sitename), position = position_dodge(width=0.3), size = 0.5, method='lm', se = FALSE)+
  geom_line(aes(group = sitename), size = 0.5, linetype="dashed", position = position_dodge(width=0.3))+
  geom_errorbar(aes(x = time, ymax = pam+se, ymin = pam-se), width = .2, position = position_dodge(width=0.3)) +
  scale_fill_manual(name = "Site",
                    labels = c("CI","PD","SP","BN","BS","CA"),
                    values = cols_site)+
  scale_color_manual(name = "Site",
                     labels = c("CI","PD","SP","BN","BS","CA"),
                     values = cols_site)+
  ylab("Fv/Fm")+
  scale_x_continuous(name = "Time Point", breaks = c(45,54,61,65,70,74,79)) +
  #ylim(-0.5,4) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))
pam_plot_site
ggsave(pam_plot_site, file="/Users/hannahaichelman/Documents/BU/TVE/PAM/PAM_site_2_lin.pdf", width=6, height=4, units=c("in"), useDingbats=FALSE)

# Plot separated by site and treatment to compare with ITS2 types
phys_pam_sites = phys_pam_lineage_plots %>%
  subset(sitename == "O-Cayo de Agua") %>%
  subset(treat == "High Var")

str(phys_pam_sites)
phys_pam_sites$frag = as.factor(phys_pam_sites$frag)
phys_pam_sites$time = as.numeric(as.character(phys_pam_sites$time))
phys_pam_sites$lineage = as.factor(phys_pam_sites$lineage)

# plot, time x axis, shape = symbiont type, color = lineage
pam_plot_lineage <- ggplot(phys_pam_sites,aes(x = time, y = pam, shape = frag, color = lineage))+
  theme_bw()+
  geom_point(size = 3, position = position_dodge(width=0.3))+
  geom_line(aes(color = lineage), size = 0.5, linetype="dashed", position = position_dodge(width=0.3))+
  #geom_errorbar(aes(x = time, ymax = pam+se, ymin = pam-se), width = .2, position = position_dodge(width=0.3)) +
  scale_color_manual(name = "Lineage",
                     labels = c("L1","L2"),
                     values = c("#fc8d59", "#99d594"))+
  scale_shape_manual(name = "Frag",
                     values=c(7,8,9,10, 11, 12, 13, 14))+
  ylab("Fv/Fm")+
  scale_x_continuous(name = "Time Point", breaks = c(46,55,62,66,71,75,80)) +
  #ylim(-0.5,4) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))
pam_plot_lineage

ggsave(pam_plot_lineage, filename = "/Users/hannahaichelman/Documents/BU/TVE/PAM/pam_CA_HighVar.pdf", width=6, height=4, units=c("in"), useDingbats=FALSE)


# include temperatures during PAM measurements
pam_temps <- read.csv("/Users/hannahaichelman/Documents/BU/TVE/PAM/pam_temp_timeline.csv")
str(pam_temps)

pam_temps$day <- as.numeric(as.character(pam_temps$day))

# plot
tempplot <-
  ggplot(pam_temps, aes(x = day, y = temp)) +
  theme_bw() +
  geom_point(size = 3) +
  geom_line() +
  ylim(27,33) +
  ylab("Temperature (°C)") +
  scale_x_continuous(name = "Time Point", breaks = c(46,55,62,66,71,75,80))
tempplot

quartz() # to be able to squish and save temp plot figure:
tempplot


#### Corallite Surface Area ####

# corallite surface area has a few extra samples from other initial phys dataframes, but
# it is still all taken from initial phys data

corrsa_phys <- read_excel("Physiology_Data/data_files/Corallite_SA_Measurements.xlsx", sheet = "data")

head(corrsa_phys)

# remove test fragments and column of notes to look at data
corrsa_phys = corrsa_phys %>%
  dplyr::filter(sample.ID != "I2A6 (test)") %>%
  dplyr::filter(sample.ID != "I2A6 (test2)") %>%
  select(sample.ID, std.area, corallite.avg.area, corallite.avg.poly) %>%
  dplyr::rename(frag = sample.ID)

# re-name treatment
corrsa_phys$treat <- c("Initial")

# add identifying data
corrsa_phys$origsitecode <- substr(corrsa_phys$frag, 1, 2)

# add in site name with inshore/offshore indicator
corrsa_phys$sitename <- ifelse(corrsa_phys$origsitecode == 'I2', 'SP',
                               ifelse(corrsa_phys$origsitecode == 'I3', 'CI',
                                      ifelse(corrsa_phys$origsitecode == 'I4', 'PD',
                                             ifelse(corrsa_phys$origsitecode == 'O2', 'BS',
                                                    ifelse(corrsa_phys$origsitecode == 'O3', 'CA',
                                                           'BN')))))
corrsa_phys$sitename <- as.factor(corrsa_phys$sitename)

# make new nubbin IDs based on the new, more informative site codes
# add inshore/offshore designation
corrsa_phys$reef <- substr(corrsa_phys$frag,1,1)
corrsa_phys$reef <- ifelse(corrsa_phys$reef == 'O', 'Outer Reef', 'Inner Reef')
corrsa_phys$reef <- as.factor(corrsa_phys$reef)

corrsa_phys$genet <- substr(corrsa_phys$frag,3,3)

#create a new column of combined genotype and site for stats later
corrsa_phys = corrsa_phys %>%
  unite(gen_site, c(origsitecode,genet), sep = "", remove = FALSE) %>%
  mutate(gen_site = as.factor(gen_site)) %>%
  dplyr::filter(gen_site != "I4G") # clone with I4F, remove from dataset

# corallite.avg.area = Area of the corallite calculated by using area of ellipse formula (corallite.min/2*corallite.max/2*pi)
# corallite.avg.poly = Area of the corallite measured in ImageJ using the polygon tool

# look at relationship between the two ways that Laura measured corallite surface area
# correlate the two measurements

plot(corrsa_phys$corallite.avg.area, corrsa_phys$corallite.avg.poly, pch = 16)
abline(lm(corrsa_phys$corallite.avg.poly ~ corrsa_phys$corallite.avg.area), col = "red", lwd = 3) # add regression line
text(paste("Correlation:", round(cor(corrsa_phys$corallite.avg.area, corrsa_phys$corallite.avg.poly), 2)), x = 10000, y = 4500) # add Pearson correlation

# they are highly correlated, so will stick with the measurement taken using the polygon tool in ImageJ: corallite.avg.poly

# now need to convert between pixel # and mm^2 for the corallite.avg.poly measurement
# cross multiply to solve for unknown number of mm^2
corrsa_phys = corrsa_phys %>%
  mutate(corallite.avg.poly.mm2 = ((corallite.avg.poly*400)/std.area))

# merge with lineage info for later plotting
lineages = read.csv("Physiology_Data/data_files/tve_lineages_noclones.csv")

corrsa_phys_all_lin <- left_join(corrsa_phys, lineages, by = "gen_site")
corrsa_phys_all_lin$lineage = as.factor(corrsa_phys_all_lin$lineage)

corrsa_phys_2_lin = corrsa_phys_all_lin %>%
  filter(is.na(lineage) | lineage!="L3") # want to keep NA values for lineage here since they still have other info, will remove na's for lineage specific plots

# Stats by lineage
str(corrsa_phys_2_lin)
corrsa_phys_2_lin$gen_site = as.factor(corrsa_phys_2_lin$gen_site)
#m1 <- lmer(corallite.avg.poly.mm2 ~ lineage + (1|gen_site), data = corrsa_phys_all_lin, REML=TRUE)
m1 <- lm(corallite.avg.poly.mm2 ~ lineage, data = corrsa_phys_2_lin)
summary(m1)
#               Estimate Std. Error t value Pr(>|t|)
# (Intercept)   8.8205     0.3375  26.134  < 2e-16 ***
# lineageL2     5.0749     0.5660   8.966 2.13e-11 ***

anova(m1)
# Response: corallite.avg.poly.mm2
#           Df Sum Sq Mean Sq F value    Pr(>F)
# lineage    1 265.56 265.558  80.386 2.129e-11 ***
# Residuals 43 142.05   3.304

#SummarySE to format data for plotting with lineage
corrsa_phys_all_lin_nona = corrsa_phys_all_lin %>%
  drop_na(lineage)

corrsa_phys_2_lin_nona = corrsa_phys_2_lin %>%
  drop_na(lineage)

corrsa_means_all_lin <- summarySE(corrsa_phys_all_lin_nona, measurevar="corallite.avg.poly.mm2", groupvars=c("lineage"))
corrsa_means_2_lin <- summarySE(corrsa_phys_2_lin_nona, measurevar="corallite.avg.poly.mm2", groupvars=c("lineage"))

# plot, lineage x axis
corrsa_plot_lineage <- ggplot(corrsa_means_2_lin,aes(x = lineage, y = corallite.avg.poly.mm2, color = lineage, fill = lineage))+
  theme_bw()+
  geom_errorbar(aes(x = lineage, ymax = corallite.avg.poly.mm2+se, ymin = corallite.avg.poly.mm2-se, color = lineage), width = .2, position = position_dodge(width=0.3)) +
  geom_point(size = 3, position = position_dodge(width=0.3), shape = 21, color = "black")+
  scale_color_manual(name = "Lineage",
                     breaks = c("L1","L2"),
                     values = cols_lineage)+
  scale_fill_manual(name = "Lineage",
                    breaks = c("L1","L2"),
                    values=cols_lineage)+
  xlab("Lineage")+
  ylab(bquote("Corallite Area (mm" ^2~')'))+
  scale_y_continuous(limits = c(6,15), breaks = seq(6,15, by = 3))+
  #geom_vline(xintercept = 1.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  theme(legend.position = "none")
corrsa_plot_lineage
ggsave(corrsa_plot_lineage, filename = "/Users/hannahaichelman/Dropbox/BU/TVE/Corallite_SA/corrsa_lineage_2lin.pdf", width=3, height=3, units=c("in"), useDingbats=FALSE)

# plot just CI corals to see if lineage difference holds
corrsa_phys_lineage_CI = corrsa_phys_2_lin_nona %>%
  subset(sitename == "CI")

# stats for CI only
str(corrsa_phys_lineage_CI)
m2 <- lm(corallite.avg.poly.mm2 ~ lineage, data = corrsa_phys_lineage_CI)
summary(m2)
# Coefficients:
#               Estimate Std. Error t value Pr(>|t|)
# (Intercept)   7.5200     0.7259  10.359 4.73e-05 ***
# lineageL2     6.0098     1.0266   5.854   0.0011 **

anova(m2)
# Response: corallite.avg.poly.mm2
#           Df Sum Sq Mean Sq F value   Pr(>F)
# lineage    1 72.234  72.234   34.27 0.001097 **

lsmeans(m2, pairwise~lineage, adjust="tukey")
# $contrasts
# contrast estimate   SE df t.ratio p.value
# L1 - L2     -6.01 1.03  6  -5.854  0.0011

corrsa_means_lineage_CI <- summarySE(corrsa_phys_lineage_CI, measurevar="corallite.avg.poly.mm2", groupvars=c("lineage"))

corrsa_plot_lineage_CI <- ggplot(corrsa_means_lineage_CI,aes(x = lineage, y = corallite.avg.poly.mm2, color = lineage, fill = lineage))+
  theme_bw()+
  geom_errorbar(aes(x = lineage, ymax = corallite.avg.poly.mm2+se, ymin = corallite.avg.poly.mm2-se), width = .2, position = position_dodge(width=0.3)) +
  geom_point(size = 3, position = position_dodge(width=0.3), shape = 21, color = "black")+
  scale_color_manual(name = "Lineage",
                     breaks = c("L1","L2"),
                     values = cols_lineage)+
  scale_fill_manual(name = "Lineage",
                    breaks = c("L1","L2"),
                    values=cols_lineage)+
  xlab("Lineage")+
  ylab(bquote("Corallite Area (mm" ^2~')'))+
  scale_y_continuous(limits = c(6,15), breaks = seq(6,15, by = 3))+
  #geom_vline(xintercept = 1.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  theme(legend.position = "none")
corrsa_plot_lineage_CI
ggsave(corrsa_plot_lineage_CI, filename = "/Users/hannahaichelman/Dropbox/BU/TVE/Corallite_SA/corrsa_lineage_CIonly.pdf", width=3, height=3, units=c("in"), useDingbats=FALSE)

#SummarySE to format data for plotting site name
corrsa_means_site_all_lin <- summarySE(corrsa_phys_all_lin, measurevar="corallite.avg.poly.mm2", groupvars=c("sitename"))
corrsa_means_site_2_lin <- summarySE(corrsa_phys_2_lin, measurevar="corallite.avg.poly.mm2", groupvars=c("sitename"))

# plot, treatment x axis colored by site data figure
corrsa_plot_site <- ggplot(corrsa_means_site_2_lin,aes(x = sitename, y = corallite.avg.poly.mm2, color = sitename, pch = sitename))+
  theme_bw()+
  geom_point(size = 3, position = position_dodge(width=0.3))+
  geom_errorbar(aes(x = sitename, ymax = corallite.avg.poly.mm2+se, ymin = corallite.avg.poly.mm2-se), width = .2, position = position_dodge(width=0.3)) +
  scale_color_manual(name = "Site",
                     labels = c("CI","PD","SP","BN","BS","CA"),
                     values = cols_site)+
  scale_shape_manual(name = "Site",
                     labels = c("CI","PD","SP","BN","BS","CA"),
                     values=c(19,19,19,17,17,17))+
  xlab("Site Name")+
  ylab(bquote("Corallite Area (mm" ^2~')'))+
  #ylim(3.75,6) +
  #geom_vline(xintercept = 1.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), legend.position = "none")
corrsa_plot_site
ggsave(corrsa_plot_site, filename = "/Users/hannahaichelman/Dropbox/BU/TVE/Corallite_SA/corrsa_site_2_lin.pdf", width=4, height=4, units=c("in"), useDingbats=FALSE)

# SummarySE for reef zone
corrsa_means_reef_all_lin <- summarySE(corrsa_phys_all_lin, measurevar="corallite.avg.poly.mm2", groupvars=c("reef"))
corrsa_means_reef_2_lin <- summarySE(corrsa_phys_2_lin, measurevar="corallite.avg.poly.mm2", groupvars=c("reef"))

# plot, reef zone x axis
corrsa_plot_reef <- ggplot(corrsa_means_reef_2_lin,aes(x = reef, y = corallite.avg.poly.mm2, color = reef, pch = reef))+
  theme_bw()+
  geom_point(size = 3, position = position_dodge(width=0.3))+
  geom_errorbar(aes(x = reef, ymax = corallite.avg.poly.mm2+se, ymin = corallite.avg.poly.mm2-se), width = .2, position = position_dodge(width=0.3)) +
  scale_color_manual(name = "Reef Zone",
                     labels = c("Inshore","Offshore"),
                     values = c("red4","royalblue4"))+
  scale_shape_manual(name = "Reef Zone",
                     labels = c("Inshore","Offshore"),
                     values=c(19,17))+
  xlab("Reef Zone")+
  ylab(bquote("Corallite Area (mm" ^2~')'))+
  #ylim(3.75,6) +
  #geom_vline(xintercept = 1.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  theme(legend.position = "none")
corrsa_plot_reef
ggsave(corrsa_plot_reef, filename = "/Users/hannahaichelman/Dropbox/BU/TVE/Corallite_SA/corrsa_reef_2_lin.pdf", width=4, height=4, units=c("in"), useDingbats=FALSE)



#### Skeleton Morphometrics ####
skel_phys = read.csv("Physiology_Data/data_files/T0_morphology.csv")

head(skel_phys)
str(skel_phys)

skel_phys2 = skel_phys %>%
  mutate(treat = as.factor(treat), sitename = as.factor(sitename), reef = as.factor(reef), lineage = as.factor(lineage))
head(skel_phys2)

#create a new column of combined genotype and site for stats later
skel_phys2$origsitecode <- substr(skel_phys2$frag, 1, 2)
skel_phys2$genet <- substr(skel_phys2$frag,3,3)

skel_phys3 = skel_phys2 %>%
  unite(gen_site, c(origsitecode,genet), sep = "", remove = FALSE) %>%
  mutate(gen_site = as.factor(gen_site)) %>%
  dplyr::filter(gen_site != "I4G") # clone with I4F, remove from dataset

# subset data for plotting DLI and LEF
skel_phys_nona = skel_phys3 %>%
  drop_na(lef)

skel_phys_CI = skel_phys_nona %>%
  dplyr::filter(sitename == "CI")

skel_phys_SP = skel_phys_nona %>%
  dplyr::filter(sitename == "SP")

skel_phys_nona_2lin = skel_phys_nona %>%
  dplyr::filter(is.na(lineage) | lineage!="L3") # want to keep NA values for lineage here since they still have other info, will remove na's for lineage specific plots

# Stats of light enhancement factor by lineage
str(skel_phys2)

m1 <- lm(lef ~ lineage, data = skel_phys_nona_2lin)
#m1 <- lm(lef ~ lineage, data = skel_phys_CI)

summary(m1)
# ALL SITES:
# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)
# (Intercept)  1.33700    0.06241  21.423  < 2e-16 ***
# lineageL2    0.80445    0.10443   7.703 2.01e-09 ***

# CI ONLY:
#             Estimate Std. Error t value Pr(>|t|)
# (Intercept)   1.2715     0.1611   7.892 0.000525 ***
# lineageL2     0.6780     0.2461   2.755 0.040069 *

anova(m1)
# ALL SITES
# Response: lef
#           Df Sum Sq Mean Sq F value    Pr(>F)
# lineage    1 6.2404  6.2404  59.337 2.006e-09 ***
# Residuals 40 4.2067  0.1052

# CI ONLY
#           Df  Sum Sq Mean Sq F value  Pr(>F)
# lineage    1 0.78809 0.78809  7.5903 0.04007 *
# Residuals  5 0.51914 0.10383


# Stats of daily light integral by lineage - not included in manuscript since we don't have DLI data for each site.
str(skel_phys_nona_2lin)

m2 <- lm(dli_kd0.38 ~ lineage, data = skel_phys_nona_2lin)
summary(m2)
# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)
# (Intercept)   19.185      1.042  18.412  < 2e-16 ***
# lineageL2     -4.879      1.744  -2.798  0.00787 **

anova(m2)
# Response: dli_kd0.38
#         Df Sum Sq Mean Sq F value   Pr(>F)
# lineage    1  229.5 229.499  7.8283 0.007873 **
# Residuals 40 1172.7  29.317


# PLOTS
# light enhancement factor plots
# summarize data for plotting
lef_means <- summarySE(skel_phys_nona, measurevar="lef", groupvars=c("lineage"))
lef_means_2lin <- summarySE(skel_phys_nona_2lin, measurevar="lef", groupvars=c("lineage"))
lef_means_SP <- summarySE(skel_phys_SP, measurevar="lef", groupvars=c("lineage"))
lef_means_CI <- summarySE(skel_phys_CI, measurevar="lef", groupvars=c("lineage"))

# plot, colored by lineage data figure
lef_plot_lineage <- ggplot(skel_phys_nona_2lin, aes(x = lineage, y = lef))+
  theme_bw()+
  geom_jitter(aes(color = lineage, fill = lineage),
              position=position_dodge(width=0.3),
              alpha=0.2, pch = 21,
              color = "black") +
  geom_errorbar(data = lef_means_2lin, aes(x = lineage, ymax = lef+se, ymin = lef-se, color = lineage), width = .2, position = position_dodge(width=0.4)) +
  geom_point(data = lef_means_2lin, mapping = aes(x=lineage, y=lef, color = lineage, fill = lineage), size = 3.5, pch = 21, color = "black", position = position_dodge(width=0.4))+
  scale_fill_manual(name = "Lineage",
                    values = cols_lineage)+
  scale_color_manual(name = "Lineage",
                     values = cols_lineage)+
  xlab("Lineage")+
  ylab("Light Enhancement Factor (LEF)") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))
lef_plot_lineage
ggsave(lef_plot_lineage, filename = "/Users/hannahaichelman/Dropbox/BU/TVE/SkeletonMorphometry/LEF_2_lin.pdf", width=4, height=4, units=c("in"), useDingbats=FALSE)

# separate out CI and SP sites and re-plot LEF:
lef_plot_lineage_CI <- ggplot(skel_phys_CI,aes(x = lineage, y = lef))+
  theme_bw()+
  geom_errorbar(data = lef_means_CI, aes(x = lineage, ymax = lef+se, ymin = lef-se, color = lineage), width = .2, position = position_dodge(width=0.4)) +
  geom_point(data = lef_means_CI, mapping = aes(x=lineage, y=lef, color = lineage, fill = lineage), size = 3.5, pch = 21, color = "black", position = position_dodge(width=0.4))+
  scale_fill_manual(name = "Lineage",
                    values = cols_lineage)+
  scale_color_manual(name = "Lineage",
                     values = cols_lineage)+
  xlab("Treatment")+
  ylab("Light Enhancement Factor (LEF)") +
  #ggtitle("Cristobal Island") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1), legend.position = "none")
lef_plot_lineage_CI
ggsave(lef_plot_lineage_CI, filename = "/Users/hannahaichelman/Documents/BU/TVE/SkeletonMorphometry/LEF_2_lin_CI.pdf", width=3, height=3, units=c("in"), useDingbats=FALSE)

#SP
# plot, treatment x axis colored by lineage data figure
lef_plot_lineage_SP <- ggplot(skel_phys_SP,aes(x = lineage, y = lef))+
  theme_bw()+
  geom_jitter(aes(color = lineage, fill = lineage),
              position=position_dodge(width=0.3),
              alpha=0.2, pch = 21,
              color = "black") +
  geom_errorbar(data = lef_means_SP, aes(x = lineage, ymax = lef+se, ymin = lef-se, color = lineage), width = .2, position = position_dodge(width=0.4)) +
  geom_point(data = lef_means_SP, mapping = aes(x=lineage, y=lef, color = lineage, fill = lineage), size = 3.5, pch = 21, color = "black", position = position_dodge(width=0.4))+
  scale_fill_manual(name = "Lineage",
                    values = cols_lineage)+
  scale_color_manual(name = "Lineage",
                     values = cols_lineage)+
  xlab("Treatment")+
  ylab("Light enhancement factor (LEF)") +
  ggtitle("STRI Point") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))
lef_plot_lineage_SP
ggsave(lef_plot_lineage_SP, filename = "/Users/hannahaichelman/Documents/BU/TVE/SkeletonMorphometry/LEF_2_lin_SP.pdf", width=5, height=4, units=c("in"), useDingbats=FALSE)

#### Skeleton morphology PCAs ####

skel_phys_pca = skel_phys3 %>%
  select(frag,treat,sitename,lineage,area_cor_cm2,den_cor_cm2,diam_cal_cm,dist_cor_cm,s_length_cm,s_thick_cm,t_thick_cm,lef) %>%
  dplyr::filter(lineage != "L3") %>%
  dplyr::filter(complete.cases(.)) #drop any row that has an NA for any time point

str(skel_phys_pca)
#write.csv(skel_phys_pca, "Physiology_Data/data_files/skel_phys_full.csv",row.names=FALSE)

skel_phys_pca_log = skel_phys_pca %>%
  mutate(den_cor_cm2 = as.numeric(den_cor_cm2)) %>%
  mutate_if(is.numeric, log)

str(skel_phys_pca_log)
#write.csv(skel_phys_pca_log, "Physiology_Data/data_files/skel_phys_full_log.csv",row.names=FALSE)

colnames(skel_phys_pca_log)[colnames(skel_phys_pca_log)=="area_cor_cm2"] <-"cor_area" # corallite area
colnames(skel_phys_pca_log)[colnames(skel_phys_pca_log)=="den_cor_cm2"] <-"cor_den" # corallite density
colnames(skel_phys_pca_log)[colnames(skel_phys_pca_log)=="diam_cal_cm"] <-"cal_diam" # calyx diameter
colnames(skel_phys_pca_log)[colnames(skel_phys_pca_log)=="dist_cor_cm"] <-"cor_dist" # distance between corallites
colnames(skel_phys_pca_log)[colnames(skel_phys_pca_log)=="s_length_cm"] <-"septa_lenth" # septa length
colnames(skel_phys_pca_log)[colnames(skel_phys_pca_log)=="s_thick_cm"] <-"septa_width" # septa thickness
colnames(skel_phys_pca_log)[colnames(skel_phys_pca_log)=="t_thick_cm"] <-"theca_thick" # theca thickness

facto_t0_skel <- PCA(skel_phys_pca_log[,5:12], scale.unit = TRUE, ncp = 10, graph = TRUE)

# PCA for T0 Physiology - color = lineage, shape = sitename
pca_t0_lineage_skel <- fviz_pca_biplot(facto_t0_skel,
                                  label = "var",
                                  col.var = "black", labelsize = 4,
                                  alpha.ind = 0) + # makes individs transparent so they can be overwritten by geom_point()
  theme_bw()+
  geom_point(aes(colour=skel_phys_pca_log$lineage, shape = skel_phys_pca_log$sitename), size = 2, stroke = 1) +
  scale_color_manual(values = cols_lineage, breaks=c("L1","L2"), name = "Lineage") +
  scale_shape_manual(values = c(15,16,17,22,21,24),
                     breaks=c("BN", "BS", "CA", "CI", "PD", "SP"),
                     labels=c("BN", "BS", "CA", "CI", "PD", "SP"),
                     name = "Site") +
  stat_ellipse(aes(color=skel_phys_pca_log$lineage), type = "t", lwd = 1)+
  labs(x = "PC1 (50.01% explained variance)",
       y = "PC2 (14.12% explained variance)") +
  theme(plot.title = element_blank(),
        axis.title = element_text(face = "bold"),
        legend.title.align =  0.5, legend.text.align = 0,
        legend.title = element_text(face = "bold"))
pca_t0_lineage_skel
ggsave(pca_t0_lineage_skel, filename = "/Users/hannahaichelman/Dropbox/BU/TVE/SkeletonMorphometry/pca_t0_skeleton_lineage.pdf", width=5.5, height=4, units=c("in"), useDingbats=FALSE)

#### Prep data for PCA ####
# You will need to run the code sections for each individual physiology parameter before running this section
# separating data frames by metric and by time point to consider
calc_pca_end = calc_phys_all_lin %>%
  select(frag,treat,sitename,lineage,T3_T0_rgr) # want to consider growth over the whole experiment

hcarb_pca_t0 = hcarb_phys_all_lin %>%
  dplyr::filter(treat=="Initial") %>%
  select(frag,gen_site,treat,sitename,lineage,hcarb_mgcm2)

hcarb_pca_end = hcarb_phys_all_lin %>%
  dplyr::filter(treat!="Initial") %>%
  select(frag,treat,sitename,lineage,hcarb_mgcm2)

scarb_pca_t0 = scarb_phys_all_lin %>%
  dplyr::filter(treat=="Initial") %>%
  select(frag,gen_site,treat,sitename,lineage,scarb_mgcm2)

scarb_pca_end = scarb_phys_all_lin %>%
  dplyr::filter(treat!="Initial") %>%
  select(frag,treat,sitename,lineage,scarb_mgcm2)

hprot_pca_t0 = hprot_phys_all_lin %>%
  dplyr::filter(treat=="Initial") %>%
  select(frag,gen_site,treat,sitename,lineage,prot_mgcm2)

hprot_pca_end = hprot_phys_all_lin %>%
  dplyr::filter(treat!="Initial") %>%
  select(frag,treat,sitename,lineage,prot_mgcm2)

tissthick_pca_t0 = tiss_phys_all_lin %>%
  select(frag,gen_site,treat,sitename,lineage,avgtiss)

chla_pca_t0 = chl_phys_all_lin %>%
  dplyr::filter(treat=="Initial") %>%
  select(frag,gen_site,treat,sitename,lineage,chlA)

chla_pca_end = chl_phys_all_lin %>%
  dplyr::filter(treat!="Initial") %>%
  select(frag,treat,sitename,lineage,chlA)

corrsa_pca_t0 = corrsa_phys_all_lin %>%
  select(frag,gen_site,treat,sitename,lineage,corallite.avg.poly.mm2)

# Initial PAM measurements were not taken, so no t0 PAM
pam_pca_end = pam_phys_all_lin %>%
  dplyr::filter(time=="79") %>% # corresponds to the closest timepoint from when phenomic metrics were sampled
  select(frag,treat,sitename,lineage,pam)

sym_pca_t0 = sym_phys_all_lin %>%
  dplyr::filter(treat=="Initial") %>%
  select(frag,gen_site,treat,sitename,lineage,sym_cm2)

sym_pca_end = sym_phys_all_lin %>%
  dplyr::filter(treat!="Initial") %>%
  select(frag,treat,sitename,lineage,sym_cm2)

# Make T0 combined data frame
t0_full = join_all(list(hcarb_pca_t0,scarb_pca_t0,hprot_pca_t0,sym_pca_t0,tissthick_pca_t0,chla_pca_t0,corrsa_pca_t0), by = c("frag","treat","sitename","lineage"), type = "full")

str(t0_full)
#write.csv(t0_full, "Physiology_Data/data_files/t0_full.csv",row.names=FALSE)

t0_full_log = t0_full %>%
  dplyr::filter(complete.cases(.)) %>% #drop any row that has an NA for any time point
  mutate_if(is.numeric, log)

str(t0_full_log)

t0_full_log$sitename <- as.factor(t0_full_log$sitename)
t0_full_log$treat <- as.factor(t0_full_log$treat)
#write.csv(t0_full_log, "Physiology_Data/data_files/t0_full_log.csv", row.names=FALSE)

# Make end of variability combined data frame
end_full = join_all(list(calc_pca_end,hcarb_pca_end,scarb_pca_end,hprot_pca_end,pam_pca_end,sym_pca_end,chla_pca_end), by = c("frag","treat","sitename","lineage"), type = "full")

str(end_full)
#write.csv(end_full, "Physiology_Data/data_files/end_full.csv", row.names = FALSE)

end_full_log = end_full %>%
  dplyr::filter(complete.cases(.)) %>% #drop any row that has an NA for any time point
  mutate(T3_T0_rgr_2 = T3_T0_rgr + 2) %>%
  select(-T3_T0_rgr) %>% #get rid of the column with negative growth values before log transforming
  mutate_if(is.numeric, log) # log transform the numeric columns

str(end_full_log)

end_full_log$sitename <- as.factor(end_full_log$sitename)
end_full_log$treat <- as.factor(end_full_log$treat)
end_full_log$lineage <- as.factor(end_full_log$lineage)
#write.csv(end_full_log, "Physiology_Data/data_files/end_full_log.csv", row.names = FALSE)

##### Make our PCAs #####
# once you have made the output files in the section above you can skip right to here to finalize filtering and make pca's

library(ggpubr)
library(ggfortify)
library(ggplot2)
library(cluster)
library(FactoMineR) # lots of options for pca visuals and summary stats
library(factoextra)
library(corrplot)
library(dplyr)
#library(ggbiplot)
library(cowplot)
library(vegan)

# try exporting data files as csv and see what we get
t0_pca = read.csv("Physiology_Data/data_files/t0_full_log.csv")
end_pca = read.csv("Physiology_Data/data_files/end_full_log.csv")

# format t0 physiology data:
str(t0_pca)

t0_pca$sitename <- as.factor(t0_pca$sitename)
t0_pca$gen_site <- as.factor(t0_pca$gen_site)
t0_pca$treat <- as.factor(t0_pca$treat)
t0_pca$lineage <- as.factor(t0_pca$lineage)
# all of the samples that have NA's for lineage are ones that didn't prep so we don't have lineage data for those genets

colnames(t0_pca)[colnames(t0_pca)=="hcarb_mgcm2"] <-"hcarb"
colnames(t0_pca)[colnames(t0_pca)=="scarb_mgcm2"] <-"scarb"
colnames(t0_pca)[colnames(t0_pca)=="prot_mgcm2"] <-"prot"
colnames(t0_pca)[colnames(t0_pca)=="sym_cm2"] <-"syms"
colnames(t0_pca)[colnames(t0_pca)=="avgtiss"] <-"tiss"
colnames(t0_pca)[colnames(t0_pca)=="corallite.avg.poly.mm2"] <-"corr_sa"

# add in dominant symbiont type dataframe
its2_types_t0 = read.csv("Physiology_Data/data_files/ITS2.dominanttype.T0.csv") %>%
  select(frag, gen_site,dominant_type)

its2_types_t0$dominant_type = as.factor(its2_types_t0$dominant_type)

t0_pca_sym <- left_join(t0_pca, its2_types_t0, by = "gen_site")
str(t0_pca_sym)

t0_pca_all_lin = t0_pca_sym %>%
  drop_na(lineage) %>%
  select(-corr_sa, -frag.y) %>% # want to remove this here because it is included in the new skeleton morphology pca incorporated during the revision
  dplyr::rename(frag = frag.x) %>%
  select(frag, gen_site, treat, sitename, lineage, dominant_type, hcarb, scarb, prot, syms, tiss, chlA)

t0_pca_2_lin = t0_pca_all_lin %>%
  dplyr::filter(is.na(lineage) | lineage!="L3") # want to keep NA values for lineage here since they still have other info, will remove na's for lineage specific plots

t0_pca_lineage_CI = t0_pca_2_lin %>%
  dplyr::filter(sitename == "CI")


# now end of variability pca data:
str(end_pca)

end_pca$sitename <- as.factor(end_pca$sitename)
end_pca$treat <- as.factor(end_pca$treat)
end_pca$lineage <- as.factor(end_pca$lineage)

# re-name columns for legibility on plots
colnames(end_pca)[colnames(end_pca)=="hcarb_mgcm2"] <-"hcarb"
colnames(end_pca)[colnames(end_pca)=="scarb_mgcm2"] <-"scarb"
colnames(end_pca)[colnames(end_pca)=="prot_mgcm2"] <-"prot"
colnames(end_pca)[colnames(end_pca)=="pam"] <-"fv_fm"
colnames(end_pca)[colnames(end_pca)=="sym_cm2"] <-"syms"
colnames(end_pca)[colnames(end_pca)=="T3_T0_rgr_2"] <-"growth"

# add in dominant symbiont type dataframe
its2_types_ps = read.csv("Physiology_Data/data_files/ITS2.dominanttype.prestress.csv") %>%
  select(frag, dominant_type)

its2_types_ps$dominant_type = as.factor(its2_types_ps$dominant_type)

end_pca_sym <- left_join(end_pca, its2_types_ps, by = "frag")
str(end_pca_sym)

end_pca_all_lin = end_pca_sym %>%
  drop_na(lineage) %>%
  drop_na(dominant_type) %>%
  select(frag, treat, sitename, lineage, dominant_type, hcarb, scarb, prot, fv_fm, syms, chlA, growth)

end_pca_2_lin = end_pca_all_lin %>%
  dplyr::filter(lineage!="L3")

end_pca_lineage_CI = end_pca_2_lin %>%
  dplyr::filter(sitename == "CI")

###
# this part of the code uses the PCA function and stats from the FactoMineR package
# T0 both lineage df's
facto_t0_all_lin <- PCA(t0_pca_all_lin[,7:12], scale.unit = TRUE, ncp = 10, graph = TRUE)
facto_t0_2_lin <- PCA(t0_pca_2_lin[,7:12], scale.unit = TRUE, ncp = 10, graph = TRUE)

facto_end_all_lin <- PCA(end_pca_all_lin[,6:12], scale.unit = TRUE, ncp = 10, graph = TRUE) # no sym type info
facto_end_2_lin <- PCA(end_pca_2_lin[,6:12], scale.unit = TRUE, ncp = 10, graph = TRUE)

#facto_t0_CI <- PCA(t0_pca_lineage_CI[,5:11], scale.unit = TRUE, ncp = 10, graph = TRUE)
#facto_end_CI <- PCA(end_pca_lineage_CI[,5:11], scale.unit = TRUE, ncp = 10, graph = TRUE)

# correlations variables - dimensions
facto_end_2_lin$var$cor
# mean of the variables
facto_end_2_lin$call$centre
# standard error of variables
facto_end_2_lin$call$ecart.type
eig.val1 <- get_eigenvalue(facto_end_2_lin)
eig.val1
fviz_eig(facto_end_2_lin, addlabels = TRUE) # get percentage of explained variances per dimension
var1 <- get_pca_ind(facto_end_2_lin)
var1
var2 <- get_pca_var(facto_end_2_lin)
var2

# Contributions to the principal components
head(var2$contrib, 10)

# make a corrplot on contribution to PC
corrplot(var2$contrib, is.corr=FALSE)

## red line = expected average contribution... "In this case, the expected average contribution (cutoff) is calculated as follow: As mentioned above, if the contributions of the 10 variables were uniform, the expected average contribution on a given PC would be 1/10 = 10%. The expected average contribution of a variable for PC1 and PC2 is : [(10* Eig1) + (10 * Eig2)]/(Eig1 + Eig2)"
# bar graphs of contributions for all of the dimensions... for variables
c <- fviz_contrib(facto_end_2_lin, choice = "var", axes = 1, top = 10)
save_plot("pc1_contribution_of_var.png",c, base_aspect_ratio = 1.6)
d <- fviz_contrib(facto_end_2_lin, choice = "var", axes = 2, top = 10)
save_plot("pc2_contribution_of_var.png",d, base_aspect_ratio = 1.6)
fviz_contrib(facto_end_2_lin, choice = "var", axes = 3, top = 10)
fviz_pca_var(facto_end_2_lin, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")
)

### describe dimensions
dim_descript <- dimdesc(facto_t0_2_lin, axes = c(1,2,3,4), proba = 0.05)
# dim1
dim_descript$Dim.1
# dim2
dim_descript$Dim.2

### make the biplot
# this fuses the facto package and ggplot features (uses base facto biplot and overlays with geom_point so can have shape demark one factor and color another factor)
fviz_pca_biplot(facto_t0_all_lin)
fviz_pca_biplot(facto_t0_2_lin)

fviz_pca_biplot(facto_end_all_lin)
fviz_pca_biplot(facto_end_2_lin)

fviz_pca_biplot(facto_end_CI)
fviz_pca_biplot(facto_t0_CI)

# these are the color schemes, defined in first section of code
cols_treat
cols_site
cols_lineage
its2_cols_greens

# PCA for T0 Physiology - color = sitename
fviz_pca_biplot(facto_t0_2_lin)
pca_t0_site <- fviz_pca_biplot(facto_t0_2_lin,
                               label = "var",
                               col.var = "black", labelsize = 4,
                               alpha.ind = 0) + # makes individs transparent so they can be overwritten by geom_point()
  theme_bw()+
  geom_point(aes(colour=t0_pca_2_lin$sitename), size = 2, stroke = 1) +
  scale_color_manual(values = cols_site,
                     #breaks=c("CI","I-Punta Donato","I-STRI Point","O-Bastimentos N","O-Bastimentos S","O-Cayo de Agua"),
                     labels=c("CI","PD","SP","BN","BS","CA"),
                     name = "Site") +
  #  stat_ellipse(geom = "polygon", type = "t", alpha = 0.2,
  #               aes(fill= t0_pca_lineage$sitename), show.legend = FALSE) + scale_fill_manual(values=cols_site) + # ellipses assumes multivariate distribution using default confidence level (0.95)
  stat_ellipse(aes(color=t0_pca_2_lin$sitename), type = "t", lwd = 1)+
  labs(x = "PC1 (53.8% explained variance)",
       y = "PC2 (15.4% explained variance)") +
  theme(plot.title = element_blank(),
        axis.title = element_text(face = "bold"),
        legend.title.align =  0.5, legend.text.align = 0,
        legend.title = element_text(face = "bold"))
pca_t0_site

ggsave(pca_t0_site, filename = "/Users/hannahaichelman/Dropbox/BU/TVE/PCAs/pca_t0_site_2_lin.pdf", width=6, height=5, units=c("in"), useDingbats=FALSE)


# PCA for T0 Physiology - color = lineage, shape = sitename
pca_t0_lineage <- fviz_pca_biplot(facto_t0_2_lin,
                                  label = "var",
                                  col.var = "black", labelsize = 4,
                                  alpha.ind = 0) + # makes individs transparent so they can be overwritten by geom_point()
  theme_bw()+
  geom_point(aes(colour=t0_pca_2_lin$lineage, shape = t0_pca_2_lin$sitename), size = 2, stroke = 1) +
  scale_color_manual(values = cols_lineage, breaks=c("L1","L2"), name = "Lineage") +
  scale_shape_manual(values = c(15,16,17,22,21,24),
                     breaks=c("BN", "BS", "CA", "CI", "PD", "SP"),
                     labels=c("BN", "BS", "CA", "CI", "PD", "SP"),
                     name = "Site") +
  #  stat_ellipse(geom = "polygon", type = "t", alpha = 0.4,
  #               aes(fill= t0_pca_lineage$lineage), show.legend = FALSE) + scale_fill_manual(values=cols_lineage) + # ellipses assumes multivariate distribution using default confidence level (0.95)
  stat_ellipse(aes(color=t0_pca_2_lin$lineage), type = "t", lwd = 1)+
  labs(x = "PC1 (53.8% explained variance)",
       y = "PC2 (15.4% explained variance)") +
  theme(plot.title = element_blank(),
        axis.title = element_text(face = "bold"),
        legend.title.align =  0.5, legend.text.align = 0,
        legend.title = element_text(face = "bold"))
pca_t0_lineage

ggsave(pca_t0_lineage, filename = "/Users/hannahaichelman/Dropbox/BU/TVE/PCAs/pca_t0_lineage_2.pdf", width=5.5, height=4, units=c("in"), useDingbats=FALSE)

# PCA for T0 Physiology CI ONLY- color = lineage
pca_t0_lineage_CI <- fviz_pca_biplot(facto_t0_CI,
                                     label = "var",
                                     col.var = "black", labelsize = 4,
                                     alpha.ind = 0) + # makes individs transparent so they can be overwritten by geom_point()
  geom_point(aes(colour=t0_pca_lineage_CI$lineage), size = 2, stroke = 1) +
  scale_color_manual(values = cols_lineage, breaks=c("L1","L2"), name = "Lineage", labels = c("SSID-A","SSID-B")) +
  #  stat_ellipse(geom = "polygon", type = "t", alpha = 0.4,
  #               aes(fill= t0_pca_lineage_CI$lineage), show.legend = FALSE) + scale_fill_manual(values=cols_lineage) + # ellipses assumes multivariate distribution using default confidence level (0.95)
  stat_ellipse(aes(color=t0_pca_lineage_CI$lineage), type = "t", lwd = 1)+
  labs(x = "PC1 (57.98% explained variance)",
       y = "PC2 (19.99% explained variance)") +
  theme(plot.title = element_blank(),
        axis.title = element_text(face = "bold"),
        legend.title.align =  0.5, legend.text.align = 0,
        legend.title = element_text(face = "bold"))
pca_t0_lineage_CI

ggsave(pca_t0_lineage_CI, filename = "/Users/hannahaichelman/Documents/BU/TVE/PCAs/pca_t0_lineage_CIonly.pdf", width=5, height=5, units=c("in"), useDingbats=FALSE)


# PCA for End of Variability Physiology - color = site, shape = site
fviz_pca_biplot(facto_end_2_lin)
pca_end_site <- fviz_pca_biplot(facto_end_2_lin,
                                label = "var",
                                col.var = "black", labelsize = 4,
                                alpha.ind = 0) + # makes individs transparent so they can be overwritten by geom_point()
  theme_bw()+
  geom_point(aes(colour=end_pca_2_lin$sitename, shape =end_pca_2_lin$sitename), size = 2, stroke = 1) +
  scale_color_manual(values = cols_site,
                     breaks=c("BN", "BS", "CA", "CI", "PD", "SP"),
                     labels=c("BN", "BS", "CA", "CI", "PD", "SP"),
                     name = "Site") +
  scale_shape_manual(values = c(15,16,17,22,21,24),
                     breaks = c("BN", "BS", "CA", "CI", "PD", "SP"),
                     labels=c("BN", "BS", "CA", "CI", "PD", "SP"),
                     name = "Site") +
  stat_ellipse(aes(color=end_pca_2_lin$sitename), type = "t", lwd = 1, show.legend = FALSE)+
  labs(x = "PC1 (46.9% explained variance)",
       y = "PC2 (17% explained variance)") +
  theme(plot.title = element_blank(),
        axis.title = element_text(face = "bold"),
        legend.title.align =  0.5, legend.text.align = 0,
        legend.title = element_text(face = "bold"))
pca_end_site

ggsave(pca_end_site, filename = "/Users/hannahaichelman/Dropbox/BU/TVE/PCAs/pca_end_site_2.pdf", width=5, height=4, units=c("in"), useDingbats=FALSE)

# PCA for End of Variability Physiology - color = treatment, shape = site
fviz_pca_biplot(facto_end_2_lin)
pca_end_treat <- fviz_pca_biplot(facto_end_2_lin,
                                 label = "var",
                                 col.var = "black", labelsize = 4,
                                 alpha.ind = 0) + # makes individs transparent so they can be overwritten by geom_point()
  theme_bw()+
  geom_point(aes(colour=end_pca_2_lin$treat, shape = end_pca_2_lin$sitename), size = 2, stroke = 1) +
  scale_color_manual(values = cols_treat,
                     #breaks=c("Control", "Low Var", "Mod Var", "High Var"),
                     #labels=c("Control", "Low Var", "Mod Var", "High Var"),
                     name = "Treatment") +
  scale_shape_manual(values = c(15,16,17,22,21,24),
                     breaks = c("BN", "BS", "CA", "CI", "PD", "SP"),
                     labels=c("BN", "BS", "CA", "CI", "PD", "SP"),
                     name = "Site") +
  #  stat_ellipse(geom = "polygon", type = "t", alpha = 0.3,
  #               aes(fill= end_pca_lineage$treat), show.legend = FALSE) + scale_fill_manual(values=cols_treat) + # ellipses assumes multivariate distribution using default confidence level (0.95)
  stat_ellipse(aes(color=end_pca_2_lin$treat), type = "t", lwd = 1)+
  labs(x = "PC1 (46.9% explained variance)",
       y = "PC2 (17% explained variance)") +
  theme(plot.title = element_blank(),
        axis.title = element_text(face = "bold"),
        legend.title.align =  0.5, legend.text.align = 0,
        legend.title = element_text(face = "bold"))
pca_end_treat

ggsave(pca_end_treat, filename = "/Users/hannahaichelman/Dropbox/BU/TVE/PCAs/pca_end_treat_2_lin.pdf", width=5, height=4, units=c("in"), useDingbats=FALSE)

# PCA for End Variability Physiology - color = lineage, shape = treatment
fviz_pca_biplot(facto_end_2_lin)
pca_end_lineage <- fviz_pca_biplot(facto_end_2_lin,
                                   label = "var",
                                   col.var = "black", labelsize = 4,
                                   alpha.ind = 0) + # makes individs transparent so they can be overwritten by geom_point()
  theme_bw()+
  geom_point(aes(colour=end_pca_2_lin$lineage, shape = end_pca_2_lin$treat), size = 2, stroke = 1) +
  scale_color_manual(values = cols_lineage, breaks=c("L1","L2"), name = "Lineage") +
  scale_shape_manual(values = c(21,22),
                     breaks = c("Control", "Mod Var"),
                     labels=c("Control", "Mod Var"),
                     name = "Treatment") +
  #  stat_ellipse(geom = "polygon", type = "t", alpha = 0.4,
  #               aes(fill= end_pca_lineage$lineage), show.legend = FALSE) + scale_fill_manual(values=cols_lineage) + # ellipses assumes multivariate distribution using default confidence level (0.95)
  stat_ellipse(aes(color= end_pca_2_lin$lineage), type = "t", lwd = 1)+
  labs(x = "PC1 (46.9% explained variance)",
       y = "PC2 (17% explained variance)") +
  theme(plot.title = element_blank(),
        axis.title = element_text(face = "bold"),
        legend.title.align =  0.5, legend.text.align = 0,
        legend.title = element_text(face = "bold"))
pca_end_lineage

ggsave(pca_end_lineage, filename = "/Users/hannahaichelman/Dropbox/BU/TVE/PCAs/pca_end_lineage_2.pdf", width=5, height=4, units=c("in"), useDingbats=FALSE)

# PCA for End Variability Physiology - color = lineage+dominant sym type, shape = sitename
its2_cols_greens = c("C1" = "#edf8e9", "C3af" = "#238b45","C3" = "#a1d99b","D1" = "#00441b")

fviz_pca_biplot(facto_end_2_lin)
pca_end_lineagesym <- fviz_pca_biplot(facto_end_2_lin,
                                      label = "var",
                                      col.var = "black", labelsize = 4,
                                      alpha.ind = 0) + # makes individs transparent so they can be overwritten by geom_point()
  theme_bw()+
  geom_point(aes(fill=end_pca_2_lin$dominant_type, shape = end_pca_2_lin$lineage), size = 2, stroke = 1) +
  scale_fill_manual(values = its2_cols_greens,
                    name = "Majority ITS2") +
  scale_shape_manual(values = c(21,24),
                     name = "Lineage") +
  stat_ellipse(aes(color= end_pca_2_lin$dominant_type), type = "t", lwd = 1, show.legend = FALSE)+
  scale_color_manual(values = its2_cols_greens,
                     name = "Majority ITS2") +
  labs(x = "PC1 (46.9% explained variance)",
       y = "PC2 (17% explained variance)") +
  theme(plot.title = element_blank(),
        axis.title = element_text(face = "bold"),
        legend.title.align =  0.5, legend.text.align = 0,
        legend.title = element_text(face = "bold"))+
  guides(color = guide_legend(override.aes = list(color = its2_cols_greens)))
pca_end_lineagesym

ggsave(pca_end_lineagesym, filename = "/Users/hannahaichelman/Dropbox/BU/TVE/PCAs/pca_end_dominantsym_2lin.pdf", width=5, height=4, units=c("in"), useDingbats=FALSE)

# PCA for End Variability Physiology CI ONLY- color = lineage
pca_end_lineage_CI <- fviz_pca_biplot(facto_end_CI,
                                      label = "var",
                                      col.var = "black", labelsize = 4,
                                      alpha.ind = 0) + # makes individs transparent so they can be overwritten by geom_point()
  geom_point(aes(colour=end_pca_lineage_CI$lineage), size = 2, stroke = 1) +
  scale_color_manual(values = cols_lineage, breaks=c("L1","L2"), name = "Lineage", labels = c("SSID-A","SSID-B")) +
  #  stat_ellipse(geom = "polygon", type = "t", alpha = 0.4,
  #               aes(fill= end_pca_lineage_CI$lineage), show.legend = FALSE) + scale_fill_manual(values=cols_lineage) + # ellipses assumes multivariate distribution using default confidence level (0.95)
  stat_ellipse(aes(color= end_pca_lineage_CI$lineage), type = "t", lwd = 1)+
  labs(x = "PC1 (36.3% explained variance)",
       y = "PC2 (21.1% explained variance)") +
  theme(plot.title = element_blank(),
        axis.title = element_text(face = "bold"),
        legend.title.align =  0.5, legend.text.align = 0,
        legend.title = element_text(face = "bold"))
pca_end_lineage_CI

ggsave(pca_end_lineage_CI, filename = "/Users/hannahaichelman/Documents/BU/TVE/PCAs/pca_end_lineage_CIonly.pdf", width=5, height=5, units=c("in"), useDingbats=FALSE)


# try faceting by sitename using ggplot - this was used for data exploration
t0_df <- prcomp(t0_pca[,5:9], scale. = TRUE)
end_df <- prcomp(end_pca_2_lin[,9:15], scale. = TRUE)

facto_end_all_lin <- PCA(end_pca_all_lin[,9:15], scale.unit = TRUE, ncp = 10, graph = TRUE) # no sym type info
facto_end_2_lin <- PCA(end_pca_2_lin[,9:15], scale.unit = TRUE, ncp = 10, graph = TRUE)

# T0 autoplot faceted by sitename
autoplot(t0_df, data = t0_pca,
         colour = "sitename",
         loadings = TRUE, loadings.colour = "black", loadings.label = TRUE, loadings.label.size = 3, loadings.label.colour = "black",
         frame = TRUE, frame.type = "norm", frame.color = "sitename") +
  scale_color_manual(name = "Site",
                     labels = c("I-Cristobal","I-Punta Donato","I-STRI Point","O-Bastimentos N","O-Bastimentos S","O-Cayo de Agua"),
                     values = palsite)+
  scale_fill_manual(name = "Site",
                    labels = c("I-Cristobal","I-Punta Donato","I-STRI Point","O-Bastimentos N","O-Bastimentos S","O-Cayo de Agua"),
                    values = palsite)+
  facet_wrap(~ sitename, ncol = 2) +
  ggtitle("T0 Physiology")

# End of variability autoplot faceted by sitename and colored by treatment
autoplot(end_df, data = end_pca_2_lin,
         colour = "dominant_type",
         loadings = TRUE, loadings.colour = "black", loadings.label = TRUE, loadings.label.size = 3, loadings.label.colour = "black",
         frame = TRUE, frame.type = "norm", frame.color = "lineage") +
  scale_color_manual(name = "Dominant Type",
                     #labels = c("Control", "Low Var", "Mod Var", "High Var"),
                     values = c("#b8e186","#d73027","#fdae61","#74add1"))+
  scale_fill_manual(name = "Dominant Type",
                    #labels = c("Control", "Low Var", "Mod Var", "High Var"),
                    values = c("#b8e186","#d73027","#fdae61","#74add1"))+
  facet_wrap(~ lineage, ncol = 2) +
  ggtitle("End of Variability Physiology") +
  theme_bw()


##### PCA Adonis Tests #####
#Use an Adonis test to get significance of factors on holobiont physiology and skeleton morphology
library(vegan)
library(MCMC.OTU)
library(MicEco)
library(pairwiseAdonis)
library(tidyverse)

t0_full = read.csv("Physiology_Data/data_files/t0_full.csv")
end_full = read.csv("Physiology_Data/data_files/end_full.csv")
skel_full = read.csv("Physiology_Data/data_files/skel_phys_full.csv")

str(t0_full)

t0_full$sitename <- as.factor(t0_full$sitename)
t0_full$treat <- as.factor(t0_full$treat)
t0_full$lineage <- as.factor(t0_full$lineage)

str(end_full)

end_full$sitename <- as.factor(end_full$sitename)
end_full$treat <- as.factor(end_full$treat)
end_full$lineage <- as.factor(end_full$lineage)

str(skel_full)
skel_full$sitename <- as.factor(skel_full$sitename)
skel_full$treat <- as.factor(skel_full$treat)
skel_full$lineage <- as.factor(skel_full$lineage)

# add in dominant symbiont type dataframe
its2_types_t0 = read.csv("Physiology_Data/data_files/ITS2.dominanttype.T0.csv") %>%
  select(frag, gen_site,dominant_type) %>%
  mutate(dominant_type = as.factor(dominant_type))

t0_full_its2 <- left_join(t0_full, its2_types_t0, by = "gen_site")
str(t0_full_its2)

its2_types_ps = read.csv("Physiology_Data/data_files/ITS2.dominanttype.prestress.csv") %>%
  select(frag, dominant_type) %>%
  mutate(dominant_type = as.factor(dominant_type))

end_full_its2 <- left_join(end_full, its2_types_ps, by = "frag")
str(end_full_its2)

# Dont want to use log transformed data here
# but need to remove NAs
t0_full_adonis = t0_full_its2 %>%
  dplyr::filter(complete.cases(.)) %>% #drop any row that has an NA for any time point
  dplyr::filter(lineage != "L3") %>% # remove L3 individuals since these aren't in our PCAs
  select(-corallite.avg.poly.mm2, -frag.y) %>% # removing here because it is now included in the skeleton morphology pca (revision)
  dplyr::rename(frag = frag.x) %>% # removing here because it is now included in the skeleton morphology pca (revision)
  select(frag, gen_site, treat, sitename, lineage, dominant_type, hcarb_mgcm2, scarb_mgcm2, prot_mgcm2, sym_cm2, avgtiss, chlA)

end_full_adonis = end_full_its2 %>%
  dplyr::filter(complete.cases(.)) %>% #drop any row that has an NA for any time point
  mutate(T3_T0_rgr_2 = T3_T0_rgr + 2) %>%
  select(-T3_T0_rgr) %>% #get rid of the column with negative growth values before log transforming, just making sure to use the same data as we do in the PCAs
  dplyr::filter(lineage != "L3") %>% # remove L3 individuals since these aren't in our PCAs
  select(frag, treat, sitename, lineage, dominant_type, T3_T0_rgr_2, hcarb_mgcm2, scarb_mgcm2, prot_mgcm2, pam, sym_cm2, chlA)

skel_full_adonis = skel_full %>%
  dplyr::filter(complete.cases(.)) %>% #drop any row that has an NA for any time point
  dplyr::filter(lineage != "L3") # remove L3 individuals since these aren't in our PCAs

end_full_adonis_L1 = end_full_adonis %>%
  filter(lineage == "L1")

end_full_adonis_L2 = end_full_adonis %>%
  filter(lineage == "L2")

# Change dataframe here based on the comparison you are interested in
nl=startedLog(data=end_full_adonis,count.columns=6:12, logstart=1) # N=46
#nl=startedLog(data=t0_full_adonis,count.columns=7:12, logstart=1) # N=42
#nl=startedLog(data=skel_full_adonis,count.columns=5:12, logstart=1) # N=42

goods.dist=vegdist(nl, method="bray", na.rm = TRUE)
goods.pcoa=pcoa(goods.dist)

# PCA:
pcp=prcomp(nl, retx=TRUE, center=TRUE)
scores=goods.pcoa$vectors
summary(goods.pcoa)
conditions=end_full_adonis[, c("frag","treat","sitename","lineage","dominant_type")] #make sure to change dataframe here
#conditions=t0_full_adonis[, c("frag","treat","sitename","lineage","dominant_type")] #make sure to change dataframe here
#conditions=skel_full_adonis[, c("frag","treat","sitename","lineage")] #make sure to change dataframe here

# PERMANOVA
head(scores)
head(conditions)

t0_model = adonis(scores~lineage, data=conditions, method="euclidean", permutations = 10000)
end_model = adonis(scores~lineage+treat, data=conditions, method="euclidean", permutations = 10000)
skel_model = adonis(scores~lineage, data=conditions, method="euclidean", permutations = 10000)

t0_output = adonis_OmegaSq(t0_model, partial = TRUE)
t0_output$aov.tab
#           Df SumsOfSqs   MeanSqs F.Model      R2 parOmegaSq Pr(>F)
# lineage    1  0.008026 0.0080262  7.6458 0.16047    0.13662  3e-04 ***
# Residuals 40  0.041990 0.0010497         0.83953
# Total     41  0.050016                   1.00000

end_output = adonis_OmegaSq(end_model, partial = TRUE)
end_output$aov.tab

#           Df SumsOfSqs  MeanSqs F.Model      R2 parOmegaSq Pr(>F)
# lineage    1  0.035560 0.035560  12.571 0.18621    0.20099  2e-04 ***
# treat      1  0.033776 0.033776  11.941 0.17687    0.19214  3e-04 ***
# Residuals 43  0.121632 0.002829         0.63692
# Total     45  0.190968                  1.00000

pairwise.adonis(end_full_adonis[,6:12], end_full_adonis$treat)
#pairs Df SumsOfSqs  F.Model         R2 p.value p.adjusted sig
# Control vs Mod Var  1 0.5863816 4.675477 0.09605406   0.013      0.013   .

pairwise.adonis(end_full_adonis[,6:12], end_full_adonis$lineage)
#pairs Df SumsOfSqs  F.Model        R2 p.value p.adjusted sig
# L2 vs L1  1  1.173853 10.47477 0.1922866   0.002      0.002   *

skel_output = adonis_OmegaSq(skel_model, partial = TRUE)
skel_output$aov.tab

# Skeleton PCA:
#           Df SumsOfSqs   MeanSqs F.Model      R2 parOmegaSq    Pr(>F)
# lineage    1  0.017035 0.0170348  10.078 0.20125    0.17773 0.0006999 ***
# Residuals 40  0.067609 0.0016902         0.79875
# Total     41  0.084644                   1.00000


#### Correlation Matrices ####
# these data are not included in the manuscript #
#source the rquery.cormat function, will need it below
source("http://www.sthda.com/upload/rquery_cormat.r")

#Function to calculate p-value of correlations for corrplot()
# mat : is a matrix of data
# ... : further arguments to pass to the native R cor.test function
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}


# use same files as for adonis test - not log transformed
t0_full = read.csv("/Users/hannahaichelman/Documents/BU/TVE/PCAs/t0_full.csv")
end_full = read.csv("/Users/hannahaichelman/Documents/BU/TVE/PCAs/end_full.csv")

# format t0 physiology data:
str(t0_full)

t0_full$sitename <- as.factor(t0_full$sitename)
t0_full$treat <- as.factor(t0_full$treat)
t0_full$reef <- as.factor(t0_full$reef)

# add in genet id and combine with lineage dataframe to include 2bRAD population data in our T0 PCAs
t0_full$gen_site <- substr(t0_full$frag,1,3)

lineages = read.csv("/Users/hannahaichelman/Documents/BU/TVE/2bRAD/Analysis/tuftscustompipeline_denovo_nosyms/tve_lineages_noclones.csv")

t0_full_lineage <- left_join(t0_full, lineages, by = "gen_site")
# all of the samples that have NA's for lineage are ones that didn't prep so we don't have lineage data for those genets

colnames(t0_full_lineage)[colnames(t0_full_lineage)=="hcarb_mgcm2"] <-"hcarb"
colnames(t0_full_lineage)[colnames(t0_full_lineage)=="scarb_mgcm2"] <-"scarb"
colnames(t0_full_lineage)[colnames(t0_full_lineage)=="prot_mgcm2"] <-"prot"
colnames(t0_full_lineage)[colnames(t0_full_lineage)=="sym_cm2"] <-"syms"
colnames(t0_full_lineage)[colnames(t0_full_lineage)=="avgtiss"] <-"tiss"
colnames(t0_full_lineage)[colnames(t0_full_lineage)=="corallite.avg.poly.mm2"] <-"corr_sa"

t0_full_all_lin = t0_full_lineage %>%
  drop_na()

t0_full_2_lin = t0_full_all_lin %>%
  filter(lineage!="L3")

# now end of variability data:
str(end_full)

end_full$sitename <- as.factor(end_full$sitename)
end_full$treat <- as.factor(end_full$treat)
end_full$reef <- as.factor(end_full$reef)

# re-level and re-name treatment
end_full$treat <- factor(end_full$treat, levels = c("Control","Low Var","Mod Var","High Var"))

# add in genet id and combine with lineage dataframe to include 2bRAD population data in our T0 PCAs
end_full$gen_site <- substr(end_full$frag,1,3)
end_full_lineage <- left_join(end_full, lineages, by = "gen_site")

# re-name columns for legibility on plots
colnames(end_full_lineage)[colnames(end_full_lineage)=="T2_T0_perc"] <-"growth"
colnames(end_full_lineage)[colnames(end_full_lineage)=="hcarb_mgcm2"] <-"hcarb"
colnames(end_full_lineage)[colnames(end_full_lineage)=="scarb_mgcm2"] <-"scarb"
colnames(end_full_lineage)[colnames(end_full_lineage)=="prot_mgcm2"] <-"prot"
colnames(end_full_lineage)[colnames(end_full_lineage)=="sym_cm2"] <-"syms"

end_full_all_lin = end_full_lineage %>%
  drop_na()

end_full_2_lin = end_full_all_lin %>%
  filter(lineage!="L3")

# add in dominant symbiont type dataframe
its2_types = read.csv("/Users/hannahaichelman/Documents/BU/TVE/16S_ITS2/ITS_PreStress_Timepoint/ITS2.dominanttype.csv") %>%
  select(frag, dominant_type) #%>%
#mutate(D1_1 = D1_sum+1) %>%
#select(-D1_sum) %>%
#mutate(propD1 = log(D1_1)) %>%
#select(-D1_1)
its2_types$dominant_type = as.factor(its2_types$dominant_type)

end_full_all_lin_symtype <- left_join(end_full_all_lin, its2_types, by = "frag")
str(end_full_all_lin_symtype)

end_full_2_lin_symtype <- left_join(end_full_2_lin, its2_types, by = "frag")
str(end_full_2_lin_symtype)

end_full_all_lin_symtype$dominant_type = as.factor(end_full_all_lin_symtype$dominant_type)
end_full_2_lin_symtype$dominant_type = as.factor(end_full_2_lin_symtype$dominant_type)


# df's
#t0_full_all_lin
#t0_full_2_lin
#end_full_all_lin_symtype
#end_full_2_lin_symtype

# make subsets for making correlation plots
t0_full_L1 = t0_full_2_lin %>%
  filter(lineage == "L1") %>%
  select(hcarb, scarb, prot, syms, tiss, chlA, corr_sa)

t0_full_L2 = t0_full_2_lin %>%
  filter(lineage == "L2") %>%
  select(hcarb, scarb, prot, syms, tiss, chlA, corr_sa)

end_full_L1_D1 = end_full_2_lin_symtype %>%
  filter(lineage == "L1") %>%
  filter(dominant_type == "D1") %>%
  select(growth, hcarb, scarb, prot, pam, syms, chlA) %>%
  drop_na()

end_full_L1_C1 = end_full_2_lin_symtype %>%
  filter(lineage == "L1") %>%
  filter(dominant_type == "C1") %>%
  select(growth, hcarb, scarb, prot, pam, syms, chlA) %>%
  drop_na()

end_full_L2 = end_full_2_lin_symtype %>%
  filter(lineage == "L2") %>%
  select(growth, hcarb, scarb, prot, pam, syms, chlA) %>%
  drop_na()

## Lineage 1 correlations
#compute correlation matrix
end_L1_m <- cor(end_full_L1) #create a correlation matrix of data
end_full_L1_p.mat <- cor.mtest(end_full_L1) #calculate significance of correlations, confidence level = 0.95

end_L1_D1_m <- cor(end_full_L1_D1) #create a correlation matrix of data
end_full_L1_D1_p.mat <- cor.mtest(end_full_L1_D1) #calculate significance of correlations, confidence level = 0.95

end_L1_C1_m <- cor(end_full_L1_C1) #create a correlation matrix of data
end_full_L1_C1_p.mat <- cor.mtest(end_full_L1_C1) #calculate significance of correlations, confidence level = 0.95

#Use the rquery.cormat() function to get a table of correlation coefficients and p-values
end_full_L1_table <- rquery.cormat(end_full_L1, type="flatten", graph=FALSE)
end_full_L1_D1_table <- rquery.cormat(end_full_L1_D1, type="flatten", graph=FALSE)
end_full_L1_C1_table <- rquery.cormat(end_full_L1_C1, type="flatten", graph=FALSE)

#plot correlation matrix
#positive correlations in blue and negative correlations are in red. Color intensity and size of the circle are proportional to the correlation coefficients
pdf(file = "/Users/hannahaichelman/Documents/BU/TVE/CorrelationMatrix/L1_end_corrmat.pdf", width=4, height=4, useDingbats=FALSE)
par(mfrow=c(1,3))
corrplot(end_L1_m, method = "circle",
         type = "upper",
         order="alphabet", #hierarchical clustering order of factors
         tl.col="black", tl.srt=45, #text label color and rotation
         p.mat = end_full_L1_p.mat, sig.level = 0.05, insig = "blank", #combine with significance
         title = "L1-all" #can add title, but adds it really high for some reason
)

corrplot(end_L1_D1_m, method = "circle",
         type = "upper",
         order="alphabet", #hierarchical clustering order of factors
         tl.col="black", tl.srt=45, #text label color and rotation
         p.mat = end_full_L1_D1_p.mat, sig.level = 0.05, insig = "blank", #combine with significance
         title = "L1-D1" #can add title, but adds it really high for some reason
)

corrplot(end_L1_C1_m, method = "circle",
         type = "upper",
         order="alphabet", #hierarchical clustering order of factors
         tl.col="black", tl.srt=45, #text label color and rotation
         p.mat = end_full_L1_C1_p.mat, sig.level = 0.05, insig = "blank", #combine with significance
         title = "L1-C1" #can add title, but adds it really high for some reason
)

dev.off()

# T0
#compute correlation matrix
t0_L1_m <- cor(t0_full_L1) #create a correlation matrix of data
t0_full_L1_p.mat <- cor.mtest(t0_full_L1) #calculate significance of correlations, confidence level = 0.95

#Use the rquery.cormat() function to get a table of correlation coefficients and p-values
t0_full_L1_table <- rquery.cormat(t0_full_L1, type="flatten", graph=FALSE)

#plot correlation matrix
#positive correlations in blue and negative correlations are in red. Color intensity and size of the circle are proportional to the correlation coefficients
pdf(file = "/Users/hannahaichelman/Documents/BU/TVE/CorrelationMatrix/L1_t0_corrmat.pdf", width=4, height=4, useDingbats=FALSE)
corrplot(t0_L1_m, method = "circle",
         type = "upper",
         order="alphabet", #hierarchical clustering order of factors
         tl.col="black", tl.srt=45, #text label color and rotation
         p.mat = t0_full_L1_p.mat, sig.level = 0.05, insig = "blank" #combine with significance
         #title = "T30" #can add title, but adds it really high for some reason
)
dev.off()

## Lineage 2 correlations
#compute correlation matrix
end_L2_m <- cor(end_full_L2) #create a correlation matrix of data
end_full_L2_p.mat <- cor.mtest(end_full_L2) #calculate significance of correlations, confidence level = 0.95

#Use the rquery.cormat() function to get a table of correlation coefficients and p-values
end_full_L2_table <- rquery.cormat(end_full_L2, type="flatten", graph=FALSE)

#plot correlation matrix
#positive correlations in blue and negative correlations are in red. Color intensity and size of the circle are proportional to the correlation coefficients
pdf(file = "/Users/hannahaichelman/Documents/BU/TVE/CorrelationMatrix/L2_end_corrmat.pdf", width=4, height=4, useDingbats=FALSE)
corrplot(end_L2_m, method = "circle",
         type = "upper",
         order="alphabet", #hierarchical clustering order of factors
         tl.col="black", tl.srt=45, #text label color and rotation
         p.mat = end_full_L2_p.mat, sig.level = 0.05, insig = "blank" #combine with significance
         #title = "T30" #can add title, but adds it really high for some reason
)
dev.off()

#t0
#compute correlation matrix
t0_L2_m <- cor(t0_full_L2) #create a correlation matrix of data
t0_full_L2_p.mat <- cor.mtest(t0_full_L2) #calculate significance of correlations, confidence level = 0.95

#Use the rquery.cormat() function to get a table of correlation coefficients and p-values
t0_full_L2_table <- rquery.cormat(t0_full_L2, type="flatten", graph=FALSE)

#plot correlation matrix
#positive correlations in blue and negative correlations are in red. Color intensity and size of the circle are proportional to the correlation coefficients
pdf(file = "/Users/hannahaichelman/Documents/BU/TVE/CorrelationMatrix/L2_t0_corrmat.pdf", width=4, height=4, useDingbats=FALSE)
corrplot(t0_L2_m, method = "circle",
         type = "upper",
         order="alphabet", #hierarchical clustering order of factors
         tl.col="black", tl.srt=45, #text label color and rotation
         p.mat = t0_full_L2_p.mat, sig.level = 0.05, insig = "blank" #combine with significance
         #title = "T30" #can add title, but adds it really high for some reason
)
dev.off()

