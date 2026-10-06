## Baseline to FU1

source("~/OneDrive/projects/Irini_HIV_Pandora/figs/code/setup.R")

script_name <- "SuppFigure_Baseline_2mo_Shift"
saved_title <- get_figure_label(script_name)

bl_fu1 <- subset(adivs, Project_Timepoint %in% c("Baseline", "Month 2")) %>%
  pivot_longer(cols = diversity$metric) %>% 
  mutate(name = factor(name, levels = diversity$metric, labels = diversity$name))

stats <- bl_fu1 %>%
  group_by(name) %>%
  do(broom.mixed::tidy(
    lmer(value ~ Project_Timepoint + (1 | subjectID), data = .),
    effect = "fixed"
  )) %>%
  filter(!term %in% "(Intercept)") %>%
  mutate(
    sig_test = paste0("Est = ", round(estimate, 2), "\nP = ", round(p.value, 2))
  )

panelA <- ggplot(
  bl_fu1,
  aes(x = Project_Timepoint, y = value, fill = Project_Timepoint, shape = Project_Timepoint)
) +
  geom_boxplot(outliers = F) +
  geom_line(aes(group = subjectID)) +
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
  scale_fill_manual(values = colors$time) +
  scale_shape_manual(values= shapes$time) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

## Panel B: beta diversity (keep all matrices or select?)

amp_bl_fu1 <- amp_subset_samples(
  amp,
  Project_Timepoint %in% c("Baseline", "Month 2")
)

bray <- amp_ordinate(
  amp_bl_fu1,
  filter_species = 0,
  type = "PCOA",
  transform = "none",
  distmeasure = "bray",
  sample_color_by = "Project_Timepoint",
  detailed_output = T
)

can <- amp_ordinate(
  amp_bl_fu1,
  filter_species = 0,
  type = "PCOA",
  transform = "none",
  distmeasure = "canberra",
  sample_color_by = "Project_Timepoint",
  detailed_output = T
)

unifrac <- amp_ordinate(
  amp_bl_fu1,
  filter_species = 0,
  type = "PCOA",
  transform = "none",
  distmeasure = "unifrac",
  sample_color_by = "Project_Timepoint",
  detailed_output = T
)

wunifrac <- amp_ordinate(
  amp_bl_fu1,
  filter_species = 0,
  type = "PCOA",
  transform = "none",
  distmeasure = "wunifrac",
  sample_color_by = "Project_Timepoint",
  detailed_output = T
)

bray_perm <- vegan::adonis2(
  bray$distmatrix ~ amp_bl_fu1$metadata$Project_Timepoint,
  permutations = 1999
)
can_perm <- vegan::adonis2(
  can$distmatrix ~ amp_bl_fu1$metadata$Project_Timepoint,
  permutations = 1999
)
unif_perm <- vegan::adonis2(
  unifrac$distmatrix ~ amp_bl_fu1$metadata$Project_Timepoint,
  permutations = 1999
)
wuf_perm <- vegan::adonis2(
  wunifrac$distmatrix ~ amp_bl_fu1$metadata$Project_Timepoint,
  permutations = 1999
)

plt_func <- function(ordination, perm_result, title) {
  ggplot(ordination$plot@data, aes(x=PCo1, y = PCo2, fill = Project_Timepoint, shape = Project_Timepoint)) +
    geom_line(aes(group = subjectID), color="grey50") +
    geom_point(size = 3) +
    stat_ellipse(aes(color = Project_Timepoint), level = 0.2) +
    scale_shape_manual(values = shapes$time) +
    scale_fill_manual(values = colors$time) +
    scale_color_manual(values = colors$time) +
    labs(x = ordination$plot@labels$x,
         y = ordination$plot@labels$y,
         fill = "Project Timepoint",
         shape = "Project Timepoint",
         color = "Project Timepoint") +
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
ggsave(paste0(fig_loc, saved_title, ".pdf"), full_plot, height = 9, width = 7, dpi = 300)
