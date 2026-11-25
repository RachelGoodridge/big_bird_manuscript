# load libraries
library(tidyverse)

# move to the directory of this script
args <- commandArgs(trailingOnly = FALSE)
script_path <- gsub("--file=", "", args[grep("--file=", args)])
if (length(script_path) > 0) setwd(dirname(script_path))

# same as the box plot in whole_data_old/con_anc_dist.R, but for the new K=9 projection
k <- 9
con_col <- "V2"
print(paste("plotting K =", k))

# get all the Daphne fortis/scandens records from the master spreadsheet
# (one row per bird, with the newest lineage classifications, so A line / B line / Big Birds aren't counted as fortis, the species corrections,
# and the first / last years alive; see ../../master.R)
source("../../master.R")
records <- load_master() %>% filter(Island == "Daphne", Species %in% c("fortis", "scandens")) %>%
  select(ind = meaningful.unique, Species, start = year, end = last.year)

# add a column for conirostris ancestry
data <- read.table(paste0("finches.", k, ".Q"))
data$ind <- read.table("finches.fam")$V2
records$con_anc <- data[[con_col]][match(records$ind, data$ind)]
records <- records[!is.na(records$con_anc) & !is.na(records$start),]

# color mapping (species in the legend written with the genus, in italics)
species_labels <- function(x) parse(text = paste0("italic('G. ", x, "')"))
species_colors <- c(fortis = "dodgerblue2", scandens = "seagreen3")

# one row per individual per year alive
df_expanded <- records %>% rowwise() %>% mutate(year = list(start:end)) %>% unnest(year) %>% filter(!year %in% unsampled_years())  # no fieldwork
box_data <- df_expanded %>% rename(value = con_anc) %>% filter(year >= 1983)

########################## BOX PLOT ##########################
# make a bunch of skinny box plots
outliers <- box_data %>% group_by(year, Species) %>% mutate(Q1 = quantile(value, 0.25), Q3 = quantile(value, 0.75), IQR = Q3 - Q1,
    is_outlier = value < Q1 - 1.5 * IQR | value > Q3 + 1.5 * IQR) %>% filter(is_outlier)

boxplot_gg <- ggplot(box_data, aes(y = value, x = factor(year))) +
    geom_hline(yintercept = pretty(box_data$value, n = 6), color = "gray90", linewidth = 0.3) +
    geom_boxplot(aes(fill = Species, color = Species), position = position_dodge(width = 0.6), width = 0.5, outlier.shape = NA, alpha = 0.4, size = 0.3) +
    geom_point(data = outliers, aes(y = value, x = factor(year), color = Species), position = position_dodge(width = 0.6), size = 1, alpha = 0.6, shape = 16) +
    scale_fill_manual(values = species_colors, name = NULL, labels = species_labels) + scale_color_manual(values = species_colors, name = NULL, labels = species_labels) + scale_x_discrete() +
    labs(y = expression(italic("G. conirostris") ~ "Ancestry"), x = NULL, title = NULL) + theme(legend.position = "right", legend.text = element_text(hjust = 0), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
    size = 11, margin = margin(t = -10), color = "black"), axis.ticks.x = element_blank(), axis.text.y = element_text(size=11, color = "black"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.title.x = element_text(size = 12))

ggsave(paste0("con_anc_boxplot_over_time_k", k, ".png"), boxplot_gg, width = 9, height = 6, dpi = 300)
