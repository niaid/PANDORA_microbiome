## Supp Figure 1

source("~/OneDrive/projects/Irini_HIV_Pandora/figs/code/setup.R")

script_name <- "SuppFigure_Heatmap"
saved_title <- get_figure_label(script_name)

## Heatmap of all Samples

library(ComplexHeatmap)

## modify species name for this figure so that when there are multiple species it reads spp instead of a long string:
amp$tax$Species2 <- ifelse(
  grepl("/", amp$tax$Species),
  "spp.",
  ifelse(amp$tax$Species %in% "none", "spp.", amp$tax$Species)
)

htmp_dat <- as.matrix(data.frame(amp$abund))
rownames(htmp_dat) <- paste0(
  "*",
  amp$tax$Genus,
  " ",
  amp$tax$Species2,
  "* (",
  rownames(htmp_dat),
  ")"
)
topSums <- names(rowSums(htmp_dat)[order(rowSums(htmp_dat), decreasing = T)[
  1:30
]])
topHtmp <- htmp_dat[topSums, ]
cols <- viridis::viridis(30)
twelveIle <- quantile(c(log(topHtmp + 1)), c(1:30 / 30))
myCol <- circlize::colorRamp2(twelveIle, cols)

pdf(file = paste0(fig_loc, saved_title, ".pdf"), height = 8, width = 10)
ht <- Heatmap(
  log(topHtmp + 1),
  name = "Top Taxa (log2)",
  col = myCol,
  row_labels = gt_render(rownames(topHtmp)),
  column_title = NULL,
  cluster_column_slices = F,
  show_column_names = F,
  row_names_max_width = max_text_width(
    rownames(topHtmp),
    gp = gpar(fontsize = 10)
  )
)
draw(ht)
dev.off()
