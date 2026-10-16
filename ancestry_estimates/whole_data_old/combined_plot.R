# load libraries
library(tidyverse)
library(patchwork)

## ----------------------------- ##
## 1. ANCESTRY OVER TIME PLOT    ##
## ----------------------------- ##

# get all the records from the master spreadsheet
records <- read.csv("/local/storage/all_birds_highcov_lowcov/Finch_Database_V2.csv", header = TRUE, fill = TRUE, na.strings = c("", "NA", "XX"), stringsAsFactors = FALSE)
records <- records %>% filter(Island == "Daphne")
records <- records[,c("meaningful.unique", "Species", "First.year", "Last.year")]
records <- records %>% filter(Species %in% c("fortis", "scandens"))
records <- records %>% rename(ind=meaningful.unique, start=First.year, end=Last.year)
records <- records %>% mutate(end = coalesce(end, start))  # if no end year, use start year
records <- records[!duplicated(records$ind),]  # remove duplicates (probably they are all re-captures)

# add a column for conirostris ancestry
data <- read.table("/SSD/reg259/expanded_big_bird_analysis/Daphne_project/finches.11.Q")
data$ind <- read.table("/SSD/reg259/expanded_big_bird_analysis/Daphne_project/finches.fam")$V2
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

# make a bunch of skinny box plots instead of the ridge plots
outliers <- ridge_data %>% group_by(year, Species) %>% mutate(Q1 = quantile(value, 0.25), Q3 = quantile(value, 0.75), IQR = Q3 - Q1,
    is_outlier = value < Q1 - 1.5 * IQR | value > Q3 + 1.5 * IQR) %>% filter(is_outlier)

# make sure all the years are integers
ridge_data <- ridge_data %>% mutate(year = as.integer(year))
outliers <- outliers %>% mutate(year = as.integer(year))

# ancestry boxplot over time
boxplot_gg <- ggplot(ridge_data, aes(y = value, x = year, group = interaction(year, Species))) +
    geom_hline(yintercept = 0.5, color="gray60", linewidth=0.4) + geom_hline(yintercept = 0.0, color="gray60", linewidth=0.4) +
    geom_hline(yintercept = pretty(ridge_data$value, n = 6), color = "gray90", linewidth = 0.3) +
    geom_boxplot(aes(fill = Species, color = Species), position = position_dodge(width = 0.6), width = 0.5, outlier.shape = NA, alpha = 0.4, size = 0.3) +
    geom_point(data = outliers, aes(y = value, x = year, color = Species), position = position_dodge(width = 0.6), size = 1, alpha = 0.6, shape = 16) +
    scale_fill_manual(values = species_colors, name = NULL) + scale_color_manual(values = species_colors, name = NULL) + 
    scale_x_continuous(limits = c(1982, 2023), breaks = seq(1980, 2025, 2), oob = scales::squish) + theme_minimal(base_size = 14) +
    scale_y_continuous(limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
    labs(y = "Mean Proportion of Conirostris Ancestry", x = NULL, title = NULL) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.text.y = element_text(size=12, color = "black"),  panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
    panel.background = element_blank(), axis.title.y = element_text(size = 14), legend.position = c(1, 0.95), legend.justification = c("right", "top"),
    legend.background = element_rect(fill = "white", color = NA), legend.key.size = unit(0.7, "lines"), legend.text = element_text(face = "italic", size = 14))

## ----------------------------- ##
## 2. IBD PLOT                   ##
## ----------------------------- ##

# get the IBD data and filter for species/years
pivot_seg <- read.csv("/SSD/reg259/IBIS/plots/pivot_seg.csv")
pivot_seg <- pivot_seg %>% filter(Species %in% c("fortis", "scandens"), Year >= 1983)

# compute negative control max error (years < 1987)
neg_con <- pivot_seg %>% filter(Year < 1987)
neg_max_error <- max(neg_con$fraction_snps_ibd, na.rm = TRUE)

# compute yearly means for each species
avg_df <- pivot_seg %>% group_by(Year, Species) %>% summarise(mean_ibd = mean(fraction_snps_ibd, na.rm = TRUE), .groups = "drop")

# scatter plot with trend lines
ibd_plot <- ggplot(pivot_seg, aes(x = Year, y = fraction_snps_ibd)) + geom_hline(yintercept = neg_max_error, color="gray40", linetype="dashed", linewidth=0.7) +
    geom_vline(xintercept = 1987, color="gray40", linetype="dashed", linewidth=0.7) + geom_hline(yintercept = 1.0, color="gray60", linewidth=0.4) +
    geom_hline(yintercept = 0.0, color="gray60", linewidth=0.4) + geom_point(aes(color = Species), size = 1.8, alpha = 0.7) +
    geom_line(data = avg_df, aes(y = mean_ibd, color = Species), linewidth = 1) + scale_color_manual(values = species_colors) +
    labs(x = NULL, y = "Fraction of SNPs in IBD with the Founder") + theme_minimal(base_size = 14) +
    scale_x_continuous(limits = c(1982, 2023), breaks = seq(1980, 2025, 2), oob = scales::squish) + scale_y_continuous(limits = c(0, 1.0), breaks = seq(0, 1.0, 0.2)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12, margin = margin(t = -10), color = "black"), axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14), axis.text = element_text(color = "black"), axis.text.y = element_text(size=12, color = "black"), legend.title = element_blank(),
    panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), legend.position = c(1, 0.95), legend.justification = c("right", "top"),
    legend.background = element_rect(fill = "white", color = NA), legend.key.size = unit(0.7, "lines"), legend.text = element_text(face = "italic", size = 14))

## ----------------------------- ##
## 3. COMBINE THE TWO FIGURES    ##
## ----------------------------- ##

ggsave("ind_avg_ibd_plot.png", ibd_plot, width = 10, height = 6, dpi = 300)

boxplot_gg <- boxplot_gg + theme(plot.margin = margin(t=2,r=5,b=-1,l=5))
ibd_plot <- ibd_plot + theme(plot.margin = margin(t=-1,r=5,b=5,l=5))
combined_plot <- boxplot_gg / ibd_plot + plot_layout(heights = c(1, 1), byrow = TRUE, design = NULL) & theme(plot.margin = margin(1, 5, 1, 5))
ggsave("combined_ancestry_ibd_plot.png", combined_plot, width = 10, height = 12, dpi = 300)