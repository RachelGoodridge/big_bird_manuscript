# load libraries
library(tidyverse)

# same as the box plot in whole_data_old/con_anc_dist.R, but for the new K=9 projection
project_dir <- "/SSD/reg259/big_bird_admixture/whole_data_new/admixture"
k <- 9
con_col <- "V2"  # conirostris ancestry cluster for K=9 (see species_ancestry in 02_plot_admixture.R)

# get all the Daphne fortis/scandens records from the master spreadsheet
records <- read.csv("/local/storage/all_birds_highcov_lowcov/Finch_Database_V2.csv", header = TRUE, fill = TRUE, na.strings = c("", "NA", "XX"), stringsAsFactors = FALSE)
records <- records[!duplicated(records$meaningful.unique),]  # remove duplicates (probably they are all re-captures)

# update the newest lineage classifications (as in 02_plot_admixture.R), so A line / B line / Big Birds aren't counted as fortis
bbs_lines <- read.csv("/local/storage/BigBirds/BigBirds.csv")
bbs_lines <- bbs_lines[!duplicated(bbs_lines$meaningful.unique),]
records$Species <- coalesce(bbs_lines$lineage[match(records$meaningful.unique, bbs_lines$meaningful.unique)], records$Species)
records <- records %>% filter(Island == "Daphne", Species %in% c("fortis", "scandens"))

# fill in missing years with the minimum year, year born, or last year, in that order of preference (as in 02_plot_admixture.R)
records <- records %>% mutate(start = coalesce(First.year, First.year.min, Year.born, Last.year), end = coalesce(Last.year, start)) %>%
  select(ind = meaningful.unique, Species, start, end)

# add a column for conirostris ancestry
data <- read.table(file.path(project_dir, paste0("finches.", k, ".Q")))
data$ind <- read.table(file.path(project_dir, "finches.fam"))$V2
records$con_anc <- data[[con_col]][match(records$ind, data$ind)]
records <- records[!is.na(records$con_anc) & !is.na(records$start),]

# color mapping
species_colors <- c(fortis = "dodgerblue2", scandens = "seagreen3")

# one row per individual per year alive
df_expanded <- records %>% rowwise() %>% mutate(year = list(start:end)) %>% unnest(year)
box_data <- df_expanded %>% rename(value = con_anc) %>% filter(year >= 1983)

########################## BOX PLOT ##########################
# make a bunch of skinny box plots
outliers <- box_data %>% group_by(year, Species) %>% mutate(Q1 = quantile(value, 0.25), Q3 = quantile(value, 0.75), IQR = Q3 - Q1,
    is_outlier = value < Q1 - 1.5 * IQR | value > Q3 + 1.5 * IQR) %>% filter(is_outlier)

boxplot_gg <- ggplot(box_data, aes(y = value, x = factor(year))) +
    geom_hline(yintercept = pretty(box_data$value, n = 6), color = "gray90", linewidth = 0.3) +
    geom_boxplot(aes(fill = Species, color = Species), position = position_dodge(width = 0.6), width = 0.5, outlier.shape = NA, alpha = 0.4, size = 0.3) +
    geom_point(data = outliers, aes(y = value, x = factor(year), color = Species), position = position_dodge(width = 0.6), size = 1, alpha = 0.6, shape = 16) +
    scale_fill_manual(values = species_colors, name = NULL) + scale_color_manual(values = species_colors, name = NULL) + scale_x_discrete() +
    labs(y = "Conirostris Ancestry", x = NULL, title = NULL) + theme(legend.position = "right", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
    size = 11, margin = margin(t = -10), color = "black"), axis.ticks.x = element_blank(), axis.text.y = element_text(size=11, color = "black"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.title.x = element_text(size = 12))

ggsave(paste0("con_anc_boxplot_over_time_k", k, ".png"), boxplot_gg, width = 9, height = 6, dpi = 300)
