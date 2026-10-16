# load libraries
library(tidyverse)
library(ggridges)
library(patchwork)

# get all the records from the master spreadsheet
records <- read.csv("/local/storage/all_birds_highcov_lowcov/Finch_Database_V2.csv", header = TRUE, fill = TRUE, na.strings = c("", "NA", "XX"), stringsAsFactors = FALSE)
records <- records %>% filter(Island == "Daphne")
records <- records[,c("meaningful.unique", "Species", "First.year", "Last.year")]
records <- records %>% filter(Species %in% c("fortis", "scandens"))
records <- records %>% rename(ind=meaningful.unique, start=First.year, end=Last.year)
records <- records %>% mutate(end = coalesce(end, start))  # if no end year, use start year
records <- records[!duplicated(records$ind),]  # remove duplicates (probably they are all re-captures)

# add a column for conirostris ancestry
data <- read.table("finches.11.Q")
data$ind <- read.table("finches.fam")$V2
names(data)[8] <- "con_anc"
records$con_anc <- data$con_anc[match(records$ind, data$ind)]
records <- records[!is.na(records$con_anc),]

# color and label mapping
species_colors <- c(fortis = "dodgerblue2", scandens = "seagreen3")
trait_names <- c("con_anc")
trait_labels <- c(con_anc = "Conirostris Ancestry")

# long format for ridge plotting
records <- records[!is.na(records$start), ]  # remove rows with no start year
df_expanded <- records %>% rowwise() %>% mutate(year = list(start:end)) %>% unnest(year)
ridge_data <- df_expanded %>% pivot_longer(cols = c(con_anc), names_to = "trait", values_to = "value") %>% filter(!is.na(value), year >= 1983)

# define a theme
ridge_theme <- theme_minimal(base_size = 14) + theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 10),
    strip.text = element_text(size = 12, face = "bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    axis.text.x = element_text(color = "black"))

# create the ridge plot function
plot_ridge_trait <- function(trait_name) {
    ridge_data %>% filter(trait == trait_name) %>% ggplot(aes(x = value, y = factor(year), fill = Species, color = Species)) +
    geom_density_ridges(scale = 2.5, rel_min_height = 0.01, linewidth = 0.3, alpha = 0.6, show.legend = FALSE) + scale_fill_manual(values = species_colors) +
    scale_color_manual(values = species_colors) + scale_y_discrete(limits = rev) + labs(x = trait_labels[[trait_name]], title = NULL) + ridge_theme}

# create the plots
ridge_plots <- map(trait_names, plot_ridge_trait)
ridge_grid <- wrap_plots(ridge_plots, nrow = 1, ncol = 1) + plot_layout(guides = "collect") & theme(legend.position = "right")
ggsave("con_anc_dist_over_time.png", ridge_grid, width = 8, height = 6, dpi = 300)

# now create a box plot of only 2022
fig <- ggplot(records %>% filter(start == 2022), aes(x = con_anc, fill = Species)) + geom_histogram(position = "identity", alpha = 0.5, binwidth = 0.01) +
    labs(x = "Conirostris Ancestry in 2022") + ridge_theme + scale_fill_manual(values = species_colors) + theme(axis.line = element_line(color = "gray"),
    axis.text.y = element_text(color = "black")) + scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0)))
ggsave("con_anc_hist_2022.png", fig, width = 8, height = 6, dpi = 300)

########################## BOX PLOT ##########################
# make a bunch of skinny box plots instead of the ridge plots
outliers <- ridge_data %>% group_by(year, Species) %>% mutate(Q1 = quantile(value, 0.25), Q3 = quantile(value, 0.75), IQR = Q3 - Q1,
    is_outlier = value < Q1 - 1.5 * IQR | value > Q3 + 1.5 * IQR) %>% filter(is_outlier)

boxplot_gg <- ggplot(ridge_data, aes(y = value, x = factor(year))) +
    geom_hline(yintercept = pretty(ridge_data$value, n = 6), color = "gray90", linewidth = 0.3) +
    geom_boxplot(aes(fill = Species, color = Species), position = position_dodge(width = 0.6), width = 0.5, outlier.shape = NA, alpha = 0.4, size = 0.3) +
    geom_point(data = outliers, aes(y = value, x = factor(year), color = Species), position = position_dodge(width = 0.6), size = 1, alpha = 0.6, shape = 16) +
    scale_fill_manual(values = species_colors, name = NULL) + scale_color_manual(values = species_colors, name = NULL) + scale_x_discrete() +
    labs(y = "Conirostris Ancestry", x = NULL, title = NULL) + theme(legend.position = "right", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
    size = 11, margin = margin(t = -10), color = "black"), axis.ticks.x = element_blank(), axis.text.y = element_text(size=11, color = "black"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.title.x = element_text(size = 12))

ggsave("con_anc_boxplot_over_time.png", boxplot_gg, width = 9, height = 6, dpi = 300)