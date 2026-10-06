# Figure 1: Any_IRIS
source("~/OneDrive/projects/Irini_HIV_Pandora/figs/code/setup.R")

script_name <- "Figure1"
saved_title <- get_figure_label(script_name)

## Panel A: alpha diversity

adivs_bl <- subset(adivs, Project_Timepoint %in% c("Month 12")) %>%
  pivot_longer(cols = diversity2$metric) %>%
  mutate(Any_IRIS = factor(Any_IRIS)) %>% 
  mutate(name = factor(name, levels = diversity2$metric, labels = diversity2$name))

stats <- adivs_bl %>%
  group_by(name) %>%
  do(broom::tidy(lm(value ~ Any_IRIS, data = .))) %>%
  filter(!term %in% "(Intercept)") %>%
  mutate(
    sig_test = paste0("Est = ", round(estimate, 2), "\nP = ", round(p.value, 2))
  )

panelA <- ggplot(adivs_bl, aes(x = Any_IRIS, y = value, fill = Any_IRIS, shape = Any_IRIS)) +
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
  scale_fill_manual(values = colors$iris) +
  scale_shape_manual(values= shapes$iris) + 
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) + 
  labs(x = "Any IRIS", y = "Alpha Diversity Value", fill = "Any IRIS", shape = "Any IRIS")

## Panel B: beta diversity (keep all matrices or select?)

amp_bl <- amp_subset_samples(amp, Project_Timepoint %in% c("Month 12"))

bray <- amp_ordinate(
  amp_bl,
  filter_species = 0,
  type = "PCOA",
  transform = "none",
  distmeasure = "bray",
  sample_color_by = "Any_IRIS",
  detailed_output = T
)

can <- amp_ordinate(
  amp_bl,
  filter_species = 0,
  type = "PCOA",
  transform = "none",
  distmeasure = "canberra",
  sample_color_by = "Any_IRIS",
  detailed_output = T
)

unifrac <- amp_ordinate(
  amp_bl,
  filter_species = 0,
  type = "PCOA",
  transform = "none",
  distmeasure = "unifrac",
  sample_color_by = "Any_IRIS",
  detailed_output = T
)

wunifrac <- amp_ordinate(
  amp_bl,
  filter_species = 0,
  type = "PCOA",
  transform = "none",
  distmeasure = "wunifrac",
  sample_color_by = "Any_IRIS",
  detailed_output = T
)

bray_perm <- vegan::adonis2(
  bray$distmatrix ~ amp_bl$metadata$Any_IRIS,
  permutations = 1999
)
can_perm <- vegan::adonis2(
  can$distmatrix ~ amp_bl$metadata$Any_IRIS,
  permutations = 1999
)
unif_perm <- vegan::adonis2(
  unifrac$distmatrix ~ amp_bl$metadata$Any_IRIS,
  permutations = 1999
)
wuf_perm <- vegan::adonis2(
  wunifrac$distmatrix ~ amp_bl$metadata$Any_IRIS,
  permutations = 1999
)

plt_func <- function(ordination, perm_result, title) {
  ggplot(ordination$plot@data, aes(x=PCo1, y = PCo2, fill = Any_IRIS, shape = Any_IRIS)) +
    geom_point(size = 3) +
    stat_ellipse(aes(color = Any_IRIS), level = 0.2) +
    scale_shape_manual(values = shapes$iris) +
    scale_fill_manual(values = colors$iris) +
    scale_color_manual(values = colors$iris) +
    labs(x = ordination$plot@labels$x,
         y = ordination$plot@labels$y,
         fill = "Any IRIS",
         shape = "Any IRIS",
         color = "Any IRIS") +
    ggtitle(title, 
            subtitle=paste0("R^2^ = ", round(perm_result$R2[1],3), "; P = ", perm_result$`Pr(>F)`[1])) +
    theme(plot.subtitle=element_markdown())
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


full_plot <- ggpubr::ggarrange(panelA, 
                               panelB, 
                               ncol = 1, 
                               heights = c(1, 1.3),
                               labels = c("A","B"),
                               common.legend = TRUE)
ggsave(paste0(fig_loc, saved_title, ".pdf"), full_plot, height = 9, width = 8.5)
