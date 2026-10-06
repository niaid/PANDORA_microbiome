# Younger vs older adiv

source("~/OneDrive/projects/Irini_HIV_Pandora/figs/code/setup.R")

script_name <- "SuppFigure_ADiv_by_Age"
saved_title <- get_figure_label(script_name)

bl_fu1 <- subset(adivs, Project_Timepoint %in% c("Baseline", "Month 2")) %>%
  pivot_longer(cols = diversity$metric) %>%
  mutate(
    name = factor(name, levels = diversity$metric, labels = diversity$name)
  )

stats <- bl_fu1 %>%
  group_by(name, median_Age) %>%
  do(broom.mixed::tidy(
    lmer(value ~ Project_Timepoint + (1 | subjectID), data = .),
    effect = "fixed"
  )) %>%
  filter(!term %in% "(Intercept)") %>%
  mutate(
    sig_test = paste0("Est = ", round(estimate, 2), "\nP = ", round(p.value, 2))
  )

int_stats <- bl_fu1 %>%
  group_by(name) %>%
  do(broom.mixed::tidy(
    lmer(value ~ Project_Timepoint * median_Age + (1 | subjectID), data = .),
    effect = "fixed"
  )) %>%
  filter(grepl(":", term)) %>%
  mutate(
    sig_test = paste0("P(int) = ", round(p.value, 2))
  )


panelA <- ggplot(
  bl_fu1,
  aes(
    x = median_Age,
    y = value,
    fill = Project_Timepoint,
    shape = Project_Timepoint
  )
) +
  geom_boxplot(outliers = F) +
  geom_point(position = position_jitterdodge()) +
  facet_wrap(~name, scales = "free_y", nrow = 1) +
  geom_text(
    data = stats,
    aes(x = median_Age, y = Inf, label = sig_test),
    vjust = 2,
    size = 3,
    inherit.aes = F
  ) +
  geom_text(
    data = int_stats,
    aes(x = 1.5, y = Inf, label = sig_test),
    vjust = 1.75,
    size = 3,
    inherit.aes = F
  ) +
  theme(legend.position = "top") +
  scale_fill_manual(values = colors$time) +
  scale_shape_manual(values = shapes$time) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(
    x = "Median Age",
    y = "Alpha Diversity Value",
    fill = "Project Timepoint",
    shape = "Project Timepoint"
  )

ggsave(
  panelA,
  filename = paste0(fig_loc, saved_title, ".pdf"),
  width = 8,
  height = 5.5,
  dpi = 300
)
