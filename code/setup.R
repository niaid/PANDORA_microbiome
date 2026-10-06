## Setup from Markdown

library(ampvis2)
library(dplyr)
library(lmerTest)
library(tidyverse)
library(ggtext)
library(ggplot2)
library(maaslin3)

set.seed(123)

theme_set(theme_bw())
home <- c("~/OneDrive/Projects/Irini_HIV_Pandora/")
fig_loc <- paste0(home, "figs/")
phy <- readRDS(paste0(home, "data/proc/anly_data.rds"))

## likely need to identify some good colors:
set1_col <- RColorBrewer::brewer.pal(8, "Set1")
colors <- list(iris = c("0" = "cadetblue", "1" = "coral3"),
               time = c("Baseline" = "#EB8A7B", "Month 2" = "#B88B84", "Month 12" = "#857775"),
               age  = c("Younger" = "#bebada", "Older" = "#8dd3c7"),
               integ = c("0" = "#fb8072", "1" = "#80b1d3"),
               nnrti = c("0" = "#c65b76", "1" = "#5b76c6"),
               myco = c("0" = "#4D4D4D", "1" = "#B2182B"))

shapes <- list(iris = c("0" = 21, "1" = 23),
               time = c("Baseline" = 21, "Month 2" = 22, "Month 12" = 23),
               age  = c("Younger" = 21, "Older" = 22))

diversity <- list(metric = c("Chao1", "Shannon", "Evenness", "PhyDiv"),
                  name = c("Chao1 Richness","Shannon Diversity","Pielou's Evenness","Faith's\nPhylogenetic Diversity"))

## requested addition of Inverse Simpson to Figure 1 only
diversity2 <- list(metric = c("Chao1", "Shannon", "Evenness", "PhyDiv", "invSimpson"),
                  name = c("Chao1 Richness","Shannon Diversity","Pielou's Evenness","Faith's\nPhylogenetic Diversity", "Inverse Simpson"))

## Some variable clean-up was needed

sample_data <- readxl::read_xlsx(paste0(
  home,
  "data/raw/PANDORA microbiome metadata updated v2_deidentified.xlsx"
)) %>%
  dplyr::rename(
    SampleID = `sample ID`,
    Enrollment_Age = `enrollment age`,
    Project_Timepoint = `project timepoint`,
    Sample_Timepoint = `sample timepoint`,
    Hispanic_Latino = `Hispanic or Latino=1`,
    WK0_BMI = `WK0 BMI`,
    WK0_CD4 = `WK0 CD4`,
    WK0_HIV_VL = `WK0 HIV VL`,
    WK0_ART_Regimen = `WK0 ART regimen`,
    NNRTI_yr1 = `NNRTI (first year)=1`,
    PI_yr1 = `PI(first year)=1`,
    Integrase_yr1 = `Integrase(first year)=1`,
    SystemicSteroids_yr1 = `On Systemic  Steroids YR1=1`,
    TB = `Tb=1`,
    Latent_TB = `Latent TB=1`,
    PCP = `PCP=1`,
    CryptococcalMeningitis = `Cryptococcal Meningiitis=1`,
    Toxoplasmosis = `Toxoplasmosis=1`,
    Histoplasmosis = `Histoplasmosis=1`,
    KS = `KS=1`,
    HBV = `HBV=1`,
    HCV = `HCV=1`,
    HBV_HCV = `HBV or HCV=1`,
    CMV = `CMV disease=1`,
    Any_IRIS = `Any type of IRIS=1`,
    FU1_HIV_VL = `FU1 HIV VL`,
    FU1_CD4 = `FU1 CD4`,
    FU1_BMI = `FU1 BMI`,
    FU2_HIV_VL = `FU2 HIV VL`,
    FU2_CD4 = `FU2 CD4`,
    FU2_BMI = `FU2 BMI`
  ) %>%
  filter(!SampleID %in% c("P230069", "P213106")) %>%
  mutate(Any_IRIS = factor(Any_IRIS)) %>%
  mutate(Integrase_yr1 = factor(Integrase_yr1)) %>%
  mutate(NNRTI_yr1 = factor(NNRTI_yr1)) %>%
  mutate(mycobacterial_infection = factor(mycobacterial_infection)) %>% 
  mutate(FU2_CD4 = as.numeric(FU2_CD4)) %>%
  mutate(FU2_BMI = as.numeric(FU2_BMI)) %>%
  mutate(WK0_BMI = as.numeric(WK0_BMI)) %>%
  mutate(WK0_CD4 = as.numeric(WK0_CD4)) %>%
  mutate(delta_BMI_FU2 = FU2_BMI - WK0_BMI) %>%
  mutate(median_Age = ifelse(Enrollment_Age < median(Enrollment_Age[Project_Timepoint %in% "baseline"]), "Younger", "Older")) %>%
  mutate(median_Age = factor(median_Age, levels = c("Younger", "Older"))) %>%
  mutate(Project_Timepoint = factor(Project_Timepoint, labels = c("Baseline", "Month 2", "Month 12"), levels = c("baseline", "fu_1", "fu_2"))) %>%
  suppressWarnings()

amp <- amp_load(
  otutable = phyloseq::otu_table(phy),
  metadata = sample_data,
  taxonomy = phyloseq::tax_table(phy),
  tree = phyloseq::phy_tree(phy)
) %>%
  suppressWarnings()

print(amp)

adivs <- suppressWarnings(amp_alpha_diversity(amp, richness = T)) %>%
  mutate(
    Evenness = Shannon / log(uniqueOTUs),
    PhyDiv = picante::pd(t(amp$abund), amp$tree, include.root = T)$PD
  ) %>%
  select(!c(ACE, Simpson, uniqueOTUs))