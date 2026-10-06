# Figure 1: Any_IRIS
source("~/OneDrive/projects/Irini_HIV_Pandora/figs/code/setup.R")

script_name <- "SuppFigure_M12_MB_NNRTI"
saved_title <- get_figure_label(script_name)

## Panel A: alpha diversity

adivs_bl <- subset(adivs, Project_Timepoint %in% c("Month 12")) %>%
  pivot_longer(cols = diversity$metric) %>%
  mutate(name = factor(name, levels = diversity$metric, labels = diversity$name))

stats <- adivs_bl %>%
  group_by(name) %>%
  do(model1 = broom::tidy(lm(value ~ NNRTI_yr1, data = .)),
     model2 = broom::tidy(lm(value ~ NNRTI_yr1 + mycobacterial_infection + SystemicSteroids_yr1, data = .))) %>%
  mutate(model1 = filter(model1, grepl("NNRTI", model1$term)),
         model2 = filter(model2, grepl("NNRTI", model2$term))) %>%
  mutate(
    sig_test = paste0("Est = ", round(model1$estimate, 2), "; P = ", round(model1$p.value, 2), "\n",
                      "Adj: Est = ", round(model2$estimate, 2), "; P= ", round(model2$p.value, 2))
  )

panelA <- ggplot(adivs_bl, aes(x = NNRTI_yr1, y = value, fill = NNRTI_yr1, shape = NNRTI_yr1)) +
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
  scale_fill_manual(values = colors$nnrti) +
  scale_shape_manual(values= shapes$iris) + 
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) + 
  labs(x = "NNRTI in Year 1", y = "Alpha Diversity Value", 
       fill = "NNRTI in Year 1", shape = "NNRTI in Year 1")

## Panel B: beta diversity

amp_12mo <- amp_subset_samples(amp, Project_Timepoint %in% c("Month 12"))

bray <- amp_ordinate(
  amp_12mo,
  filter_species = 0,
  type = "PCOA",
  transform = "none",
  distmeasure = "bray",
  sample_color_by = "NNRTI_yr1",
  detailed_output = T
)

can <- amp_ordinate(
  amp_12mo,
  filter_species = 0,
  type = "PCOA",
  transform = "none",
  distmeasure = "canberra",
  sample_color_by = "NNRTI_yr1",
  detailed_output = T
)

unifrac <- amp_ordinate(
  amp_12mo,
  filter_species = 0,
  type = "PCOA",
  transform = "none",
  distmeasure = "unifrac",
  sample_color_by = "NNRTI_yr1",
  detailed_output = T
)

wunifrac <- amp_ordinate(
  amp_12mo,
  filter_species = 0,
  type = "PCOA",
  transform = "none",
  distmeasure = "wunifrac",
  sample_color_by = "NNRTI_yr1",
  detailed_output = T
)

bray_perm <- vegan::adonis2(
  bray$distmatrix ~ amp_12mo$metadata$NNRTI_yr1,
  permutations = 1999
)
can_perm <- vegan::adonis2(
  can$distmatrix ~ amp_12mo$metadata$NNRTI_yr1,
  permutations = 1999
)
unif_perm <- vegan::adonis2(
  unifrac$distmatrix ~ amp_12mo$metadata$NNRTI_yr1,
  permutations = 1999
)
wuf_perm <- vegan::adonis2(
  wunifrac$distmatrix ~ amp_12mo$metadata$NNRTI_yr1,
  permutations = 1999
)

bray_perm_adj <- vegan::adonis2(
  bray$distmatrix ~ amp_12mo$metadata$NNRTI_yr1 + amp_12mo$metadata$mycobacterial_infection + amp_12mo$metadata$SystemicSteroids_yr1,
  permutations = 1999, by="margin"
)
can_perm_adj <- vegan::adonis2(
  can$distmatrix ~ amp_12mo$metadata$NNRTI_yr1 + amp_12mo$metadata$mycobacterial_infection + amp_12mo$metadata$SystemicSteroids_yr1,
  permutations = 1999, by="margin"
)
unif_perm_adj <- vegan::adonis2(
  unifrac$distmatrix ~ amp_12mo$metadata$NNRTI_yr1 + amp_12mo$metadata$mycobacterial_infection + amp_12mo$metadata$SystemicSteroids_yr1,
  permutations = 1999, by="margin"
)
wuf_perm_adj <- vegan::adonis2(
  wunifrac$distmatrix ~ amp_12mo$metadata$NNRTI_yr1 + amp_12mo$metadata$mycobacterial_infection + amp_12mo$metadata$SystemicSteroids_yr1,
  permutations = 1999, by="margin"
)

plt_func <- function(ordination, perm_result, perm_adj, title) {
  ggplot(ordination$plot@data, aes(x=PCo1, y = PCo2, fill = NNRTI_yr1, shape = NNRTI_yr1)) +
    geom_point(size = 3) +
    stat_ellipse(aes(color = NNRTI_yr1), level = 0.2) +
    scale_shape_manual(values = shapes$iris) +
    scale_fill_manual(values = colors$nnrti) +
    scale_color_manual(values = colors$nnrti) +
    labs(x = ordination$plot@labels$x,
         y = ordination$plot@labels$y,
         fill = "NNRTI in Year 1",
         shape = "NNRTI in Year 1",
         color = "NNRTI in Year 1") +
    ggtitle(title, 
            subtitle=paste0("R^2^ = ", round(perm_result$R2[1],3), "; P = ", perm_result$`Pr(>F)`[1], "<br>",
                       "Adj: R^2^ = ", round(perm_adj$R2[1],3), "; P = ", perm_adj$`Pr(>F)`[1])) +
    theme(plot.subtitle=element_markdown())
}

panelB <- ggpubr::ggarrange(
  plt_func(bray, bray_perm, bray_perm_adj, "Bray Curtis"),
  plt_func(can, can_perm, can_perm_adj, "Canberra"),
  plt_func(unifrac, unif_perm, unif_perm_adj, "Unweighted UniFrac"),
  plt_func(wunifrac, wuf_perm, wuf_perm_adj, "Weighted UniFrac"),
  ncol = 2,
  nrow = 2,
  common.legend = T
)


full_plot <- ggpubr::ggarrange(panelA, 
                               panelB, 
                               ncol = 1, 
                               heights = c(1, 1.3),
                               labels = c("A","B"),
                               common.legend = TRUE)
ggsave(paste0(fig_loc, saved_title, ".pdf"), full_plot, height = 9, width = 7, dpi = 300)
