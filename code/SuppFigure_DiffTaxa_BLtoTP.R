## Supplemental Figure 5
## Differential taxa analysis from Baseline to Month 2 and month 12

script_name <- "SuppFigure_DiffTaxa_BLtoTP"
saved_title <- get_figure_label(script_name)

source("~/OneDrive/projects/Irini_HIV_Pandora/figs/code/setup.R")

amp_bl_m2 <- amp_subset_samples(amp, Project_Timepoint %in% c("Baseline", "Month 2"))
amp_bl_m12 <- amp_subset_samples(amp, Project_Timepoint %in% c("Baseline", "Month 12"))

## Baseline to Month 2
run_maaslin_bl_m2 <- maaslin3(
  input_data = amp_bl_m2$abund,
  input_metadata = amp_bl_m2$metadata,
  output = tempfile(),
  fixed_effects = c("Project_Timepoint"),
  random_effects = c("subjectID"),
  normalization = "NONE",
  transform = "LOG",
  standardize = FALSE,
  min_prevalence = 0.1,
  warn_prevalence = F,
  plot_associations = FALSE,
  plot_summary_plot = FALSE,
  verbosity = "ERROR"
)
res_bl_m2 <- left_join(
  run_maaslin_bl_m2$fit_data_abundance$results,
  amp_bl_m2$tax,
  by = c("feature" = "OTU")
) %>%
  filter(is.na(error)) %>% 
  mutate(Genus = ifelse(Genus == "none", Family, Genus)) %>% 
  mutate(color = ifelse(qval_joint < 0.2, as.character(fct_lump_n(Genus, n = 15)), NA)) %>%
  mutate(color = ifelse(is.na(color), "Not Significant", as.character(color))) %>% 
  mutate(color = fct_relevel(color, c("Other", "Not Significant"), after = Inf))

genus_m2 <- RColorBrewer::brewer.pal(length(levels(res_bl_m2$color)), "Accent")
names(genus_m2) <- levels(res_bl_m2$color)

genus_m2["Not Significant"] <- "white"
genus_m2["Other"] <- "grey"

italic_labels <- ifelse(names(genus_m2) %in% c("Not Significant", "Other"), names(genus_m2), paste0("*", names(genus_m2), "*"))
names(italic_labels) <- names(genus_m2)

panelA <- ggplot(
  res_bl_m2,
  aes(x = coef, y = pval_joint, fill = color)
) +
  scale_y_continuous(trans = c("log10", "reverse")) + 
  scale_x_continuous(limits = c(-5, 5)) +
  geom_point(shape = 21, size = 2) +
  scale_fill_manual(values = genus_m2, labels = italic_labels) + 
  labs(x = "Coefficient", y = "P-value (log10 scale)", fill = "") + 
  theme(legend.text = element_markdown())


# baseline to month 12
run_maaslin_bl_m12 <- maaslin3(
  input_data = amp_bl_m12$abund,
  input_metadata = amp_bl_m12$metadata,
  output = tempfile(),
  fixed_effects = c("Project_Timepoint"),
  random_effects = c("subjectID"),
  normalization = "NONE",
  transform = "LOG",
  standardize = FALSE,
  min_prevalence = 0.1,
  warn_prevalence = F,
  plot_associations = FALSE,
  plot_summary_plot = FALSE,
  verbosity = "ERROR"
)
res_bl_m12 <- left_join(
  run_maaslin_bl_m12$fit_data_abundance$results,
  amp_bl_m12$tax,
  by = c("feature" = "OTU")
) %>%
  filter(is.na(error)) %>% 
  mutate(Genus = ifelse(Genus == "none", Family, Genus)) %>%
  mutate(Genus = ifelse(Family %in% "Ruminococcaceae", "Ruminococcaceae", Genus)) %>% 
  mutate(color = ifelse(qval_joint < 0.2, as.character(fct_lump_n(Genus, n = 10)), NA)) %>%
  mutate(color = ifelse(is.na(color), "Not Significant", as.character(color))) %>% 
  mutate(color = fct_relevel(color, c("Other", "Not Significant"), after = Inf))

# allows for two fewer colors to be requested, which get overwritten anyway
genus_m12 <- levels(res_bl_m12$color)
genus_m12 <- RColorBrewer::brewer.pal(length(genus_m12)-2, "Set2")
names(genus_m12) <- levels(res_bl_m12$color)[!levels(res_bl_m12$color) %in% c("Not Significant", "Other")]

genus_m12["Not Significant"] <- "white"
genus_m12["Other"] <- "grey"

italic_labels <- ifelse(names(genus_m12) %in% c("Not Significant", "Other"), names(genus_m12), paste0("*", names(genus_m12), "*"))
names(italic_labels) <- names(genus_m12)

panelB <- ggplot(
  res_bl_m12,
  aes(x = coef, y = pval_joint, fill = color, label = paste(Genus, Species), label2 = feature)
) +
  scale_y_continuous(trans = c("log10", "reverse")) +
  scale_x_continuous(limits = c(-5, 5)) +
  geom_point(shape = 21, size = 2) +
  scale_fill_manual(values = genus_m12, labels = italic_labels) + 
  labs(x = "Coefficient", y = "P-value (log10 scale)", fill = "") + 
  theme(legend.text = element_markdown())


full_plot <- ggpubr::ggarrange(panelA, 
                               panelB, 
                               ncol = 1, 
                               heights = c(1, 1),
                               labels = c("A","B"))

## Save plot
ggsave(paste0(fig_loc, saved_title[1], ".pdf"), full_plot, height = 9, width = 7, dpi = 300)

## publication columns for tables:
res_bl_m2 <- res_bl_m2 %>% 
  arrange(desc(qval_joint)) %>% 
  select(feature, coef, stderr, pval_joint, qval_joint, N, N_not_zero, Kingdom:Species)
names(res_bl_m2) <- diff_abund_table_names
res_bl_m12 <- res_bl_m12 %>% 
  arrange(desc(qval_joint)) %>%
  select(feature, coef, stderr, pval_joint, qval_joint, N, N_not_zero, Kingdom:Species)
names(res_bl_m12) <- diff_abund_table_names
writexl::write_xlsx(list("Baseline to Month 2" = res_bl_m2), path = paste0(fig_loc, saved_title[2], ".xlsx"))
writexl::write_xlsx(list("Baseline to Month 12" = res_bl_m12), path = paste0(fig_loc, saved_title[3], ".xlsx"))