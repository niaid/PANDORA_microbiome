## Supplemental Figure 12
# Month 12 alpha and beta diversity by mycobacterial coinfection

source("~/OneDrive/projects/Irini_HIV_Pandora/figs/code/setup.R")

script_name <- "SuppFigure_M12_MB_Myco"
saved_title <- get_figure_label(script_name)

## Panel A: alpha diversity

adivs_bl <- subset(adivs, Project_Timepoint %in% c("Month 12")) %>%
  pivot_longer(cols = diversity$metric) %>%
  mutate(mycobacterial_infection = factor(mycobacterial_infection)) %>%
  mutate(
    name = factor(name, levels = diversity$metric, labels = diversity$name)
  )

stats <- adivs_bl %>%
  group_by(name) %>%
  do(broom::tidy(lm(value ~ mycobacterial_infection, data = .))) %>%
  filter(!term %in% "(Intercept)") %>%
  mutate(
    sig_test = paste0(
      "Est = ",
      round(estimate, 2),
      "\nP ",
      ifelse(round(p.value, 2) < 0.01, "< 0.01", paste("= ", round(p.value, 2)))
    )
  )

panelA <- ggplot(
  adivs_bl,
  aes(
    x = mycobacterial_infection,
    y = value,
    fill = mycobacterial_infection,
    shape = mycobacterial_infection
  )
) +
  geom_boxplot(outliers = F) +
  geom_point() +
  facet_wrap(~name, scales = "free_y", nrow = 1) +
  geom_text(
    data = stats,
    aes(x = 1.5, y = Inf, label = sig_test),
    vjust = 1.25,
    size = 3,
    inherit.aes = F
  ) +
  theme(legend.position = "top") +
  scale_fill_manual(values = colors$myco) +
  scale_shape_manual(values = shapes$iris) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(
    x = "Mycobacterial Infection",
    y = "Alpha Diversity Value",
    fill = "Mycobacterial Infection",
    shape = "Mycobacterial Infection"
  )

## Panel B: beta diversity

amp_bl <- amp_subset_samples(amp, Project_Timepoint %in% c("Month 12"))

bray <- amp_ordinate(
  amp_bl,
  filter_species = 0,
  type = "PCOA",
  transform = "none",
  distmeasure = "bray",
  sample_color_by = "mycobacterial_infection",
  detailed_output = T
)

can <- amp_ordinate(
  amp_bl,
  filter_species = 0,
  type = "PCOA",
  transform = "none",
  distmeasure = "canberra",
  sample_color_by = "mycobacterial_infection",
  detailed_output = T
)

unifrac <- amp_ordinate(
  amp_bl,
  filter_species = 0,
  type = "PCOA",
  transform = "none",
  distmeasure = "unifrac",
  sample_color_by = "mycobacterial_infection",
  detailed_output = T
)

wunifrac <- amp_ordinate(
  amp_bl,
  filter_species = 0,
  type = "PCOA",
  transform = "none",
  distmeasure = "wunifrac",
  sample_color_by = "mycobacterial_infection",
  detailed_output = T
)

bray_perm <- vegan::adonis2(
  bray$distmatrix ~ amp_bl$metadata$mycobacterial_infection,
  permutations = 1999
)
can_perm <- vegan::adonis2(
  can$distmatrix ~ amp_bl$metadata$mycobacterial_infection,
  permutations = 1999
)
unif_perm <- vegan::adonis2(
  unifrac$distmatrix ~ amp_bl$metadata$mycobacterial_infection,
  permutations = 1999
)
wuf_perm <- vegan::adonis2(
  wunifrac$distmatrix ~ amp_bl$metadata$mycobacterial_infection,
  permutations = 1999
)

plt_func <- function(ordination, perm_result, title) {
  ggplot(
    ordination$plot@data,
    aes(
      x = PCo1,
      y = PCo2,
      fill = mycobacterial_infection,
      shape = mycobacterial_infection
    )
  ) +
    geom_point(size = 3) +
    stat_ellipse(aes(color = mycobacterial_infection), level = 0.2) +
    scale_shape_manual(values = shapes$iris) +
    scale_fill_manual(values = colors$myco) +
    scale_color_manual(values = colors$myco) +
    labs(
      x = ordination$plot@labels$x,
      y = ordination$plot@labels$y,
      fill = "Mycobacterial Infection",
      shape = "Mycobacterial Infection",
      color = "Mycobacterial Infection"
    ) +
    ggtitle(
      title,
      subtitle = paste0(
        "R^2^ = ",
        round(perm_result$R2[1], 3),
        "; P = ",
        perm_result$`Pr(>F)`[1]
      )
    ) +
    theme(plot.subtitle = element_markdown())
}

panelB <- ggpubr::ggarrange(
  plt_func(bray, bray_perm, "Bray Curtis"),
  plt_func(can, can_perm, "Canberra"),
  plt_func(unifrac, unif_perm, "Unweighted UniFrac"),
  plt_func(wunifrac, wuf_perm, "Weighted UniFrac"),
  ncol = 2,
  nrow = 2,
  common.legend = T
)


full_plot <- ggpubr::ggarrange(
  panelA,
  panelB,
  ncol = 1,
  heights = c(1, 1.3),
  labels = c("A", "B"),
  common.legend = TRUE
)


## Differential Abundance

run_maaslin_iris <- maaslin3(
  input_data = amp_bl$abund,
  input_metadata = amp_bl$metadata,
  output = tempfile(),
  fixed_effects = c("mycobacterial_infection"),
  normalization = "NONE",
  transform = "LOG",
  standardize = FALSE,
  min_prevalence = 0.1,
  warn_prevalence = F,
  plot_associations = FALSE,
  plot_summary_plot = FALSE,
  verbosity = "ERROR"
)
diff_abund_table <- left_join(
  run_maaslin_iris$fit_data_abundance$results,
  amp_bl$tax,
  by = c("feature" = "OTU")
) %>%
  select(
    feature,
    coef,
    stderr,
    pval_joint,
    qval_joint,
    N,
    N_not_zero,
    Kingdom:Species
  ) %>%
  arrange(desc(qval_joint))

names(diff_abund_table) <- diff_abund_table_names

ggsave(
  paste0(fig_loc, saved_title[1], ".pdf"),
  full_plot,
  height = 9,
  width = 8.5
)
writexl::write_xlsx(
  list("Mycobacterial Infection - MO12" = diff_abund_table),
  path = paste0(fig_loc, saved_title[2], ".xlsx")
)
