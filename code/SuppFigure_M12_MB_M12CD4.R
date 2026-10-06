source("~/OneDrive/projects/Irini_HIV_Pandora/figs/code/setup.R")

script_name <- "SuppFigure_M12_MB_M12CD4"
saved_title <- get_figure_label(script_name)

adivs_bl <- subset(adivs, Project_Timepoint %in% c("Month 12")) %>%
  pivot_longer(cols = diversity$metric) %>%
  mutate(
    name = factor(name, levels = diversity$metric, labels = diversity$name)
  )

stats <- adivs_bl %>%
  group_by(name) %>%
  do(broom::tidy(lm(value ~ FU2_CD4, data = .))) %>%
  filter(!term %in% "(Intercept)") %>%
  mutate(
    sig_test = paste0("Est = ", round(estimate, 3), "\nP = ", round(p.value, 2))
  )

panelA <- ggplot(adivs_bl, aes(x = FU2_CD4, y = value)) +
  geom_point(shape = 21, fill = "lightblue", size = 3) +
  facet_wrap(~name, scales = "free_y", nrow = 1) +
  geom_text(
    data = stats,
    aes(x = 350, y = Inf, label = sig_test),
    vjust = 1.25,
    size = 3,
    inherit.aes = F
  ) +
  theme(legend.position = "top") +
  geom_smooth(method = "lm", fill = "lightblue", color = "darkblue") +
  labs(x = "CD4 at Month 12", y = "Alpha Diversity Value")

## start panel B
amp_bl <- amp_subset_samples(amp, Project_Timepoint %in% c("Month 12"))

bray <- amp_ordinate(
  amp_bl,
  filter_species = 0,
  type = "PCOA",
  transform = "none",
  distmeasure = "bray",
  sample_color_by = "FU2_CD4",
  detailed_output = T
)

can <- amp_ordinate(
  amp_bl,
  filter_species = 0,
  type = "PCOA",
  transform = "none",
  distmeasure = "canberra",
  sample_color_by = "FU2_CD4",
  detailed_output = T
)

unifrac <- amp_ordinate(
  amp_bl,
  filter_species = 0,
  type = "PCOA",
  transform = "none",
  distmeasure = "unifrac",
  sample_color_by = "FU2_CD4",
  detailed_output = T
)

wunifrac <- amp_ordinate(
  amp_bl,
  filter_species = 0,
  type = "PCOA",
  transform = "none",
  distmeasure = "wunifrac",
  sample_color_by = "FU2_CD4",
  detailed_output = T
)

bray_perm <- vegan::adonis2(
  bray$distmatrix ~ amp_bl$metadata$FU2_CD4,
  permutations = 1999
)
can_perm <- vegan::adonis2(
  can$distmatrix ~ amp_bl$metadata$FU2_CD4,
  permutations = 1999
)
unif_perm <- vegan::adonis2(
  unifrac$distmatrix ~ amp_bl$metadata$FU2_CD4,
  permutations = 1999
)
wuf_perm <- vegan::adonis2(
  wunifrac$distmatrix ~ amp_bl$metadata$FU2_CD4,
  permutations = 1999
)

plt_func <- function(ordination, perm_result, title) {
  ggplot(ordination$plot@data, aes(x = PCo1, y = PCo2, fill = FU2_CD4)) +
    geom_point(size = 3, shape = 21) +
    scale_fill_viridis_c(option = "magma") +
    labs(
      x = ordination$plot@labels$x,
      y = ordination$plot@labels$y,
      fill = "CD4 at Month 12"
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

ggsave(
  paste0(fig_loc, saved_title, ".pdf"),
  full_plot,
  height = 9,
  width = 7,
  dpi = 300
)
