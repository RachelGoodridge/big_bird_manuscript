# Load required libraries
library(tidyverse)
library(ggpubr)
library(gridExtra)
library(grid)

# Load the records and admixture data for k=5
records <- read.table("/SSD/reg259/Daphne_data_analysis/older_Daphne_data/Daphne_neutral_aims_pruned_Daphne_AIMs.csv", sep=",", header=TRUE)
data <- read.table("finches.11.Q")
fam <- read.table("finches.fam")

# Assign individual IDs and column names
data$ind <- fam$V2
data <- data %>% rename(new_anc_fort = V9, new_anc_scan = V3, new_anc_fuli = V4, new_anc_mag = V7, new_anc_mag2 = V6)
data <- subset(data, select = -c(V1, V2, V5, V8, V10, V11))

# Match and add species info
data$Species <- records$Species[match(data$ind, records$sample)]
data <- data[!is.na(data$Species), ]

# Add old ancestry values
ancestries <- c("fort", "scan", "mag", "fuli", "mag2")
for (anc in ancestries) {
  data[[paste0("old_anc_", anc)]] <- records[[paste0("anc_", anc)]][match(data$ind, records$sample)]
}
data <- data[!is.na(data$old_anc_mag), ]  # filter complete cases

# Generate plots
for (sp in c("fortis", "scandens")) {
  subset_data <- data[data$Species == sp, ]
  plots <- list()

  for (anc in ancestries) {
    old_col <- paste0("old_anc_", anc)
    new_col <- paste0("new_anc_", anc)

    p <- ggplot(subset_data, aes(x = .data[[old_col]], y = .data[[new_col]])) + labs(x = paste("Old", anc, "anc prop"), y = paste("New", anc, "anc prop")) +
      theme_minimal(base_size = 15) + geom_abline(slope = 1, intercept = 0, color = "red") + geom_point() + stat_cor(method = "pearson") +
      scale_x_continuous(limits = c(0, 1)) + scale_y_continuous(limits = c(0, 1)) + theme(plot.margin = margin(0.3, 0.3, 0.3, 0.3, unit = "cm"))
    plots[[length(plots) + 1]] <- p
  }

  # Arrange in 2 rows x 3 columns layout
  png(file = paste0("ancestry_comparison_", sp, ".png"), width = 2400, height = 1600, res = 200)
  grid.arrange(grobs = plots, ncol = 3, nrow = 2, top = textGrob(paste("G.", sp), gp = gpar(fontsize = 20, fontface = "italic")))
  dev.off()
}
