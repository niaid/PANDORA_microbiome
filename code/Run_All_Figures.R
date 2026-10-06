rm(list = ls())

home <- c("~/OneDrive/Projects/Irini_HIV_Pandora/")
fig_loc <- paste0(home, "figs/")

## this file (figure_map) is really only helpful when running all figures at once
## so I am moving the figure_map and get_figure_label function here
# so that it can be used by all figure scripts when I am running everything together.
figure_map <- read.csv(paste0(fig_loc, "figure_code_linker.csv"))

get_figure_label <- function(script_name) {
  label <- figure_map$figureID[figure_map$code %in% script_name]
  if (length(label) == 0) {
    stop(paste0("No figure label found for script: ", script_name))
  }
  return(label)
}

diff_abund_table_names <- c(
  "ASV",
  "Coef",
  "Std.Err.",
  "P-value",
  "Q-value",
  "N",
  "non-Zero N",
  "Kingdom",
  "Phylum",
  "Class",
  "Order",
  "Family",
  "Genus",
  "Species"
)

for (script in unique(figure_map$code)) {
  source(paste0(fig_loc, "code/", script, ".R"))
}

writeLines(
  c(
    capture.output(sessioninfo::session_info()),
    paste("Bioconductor version:", as.character(BiocManager::version()))
  ),
  paste0(fig_loc, "sessionInfo.txt")
)
