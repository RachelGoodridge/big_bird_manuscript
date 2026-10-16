# switch to working in R and define inputs
library(tidyverse)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(grid)

# get records from the newer Daphne dataset
records1 <- read.table("/SSD/reg259/finch_genome_imputation/final_project/Daphne_finches_band_records_2022_CLEANED.csv", sep=",", header=TRUE)
records1 <- records1[,c("meaningful.unique", "Species", "First.year")]
records1 <- records1 %>% rename(ind=meaningful.unique, start=First.year)
records1$end <- records1$start
records1 <- records1[!duplicated(records1$ind), ]  # remove duplicates

# get records from the older Daphne dataset
records2 <- read.csv("/SSD/reg259/all_birds_highcov_lowcov/Finch_Database_V2.csv", header = TRUE, fill = TRUE, na.strings = c("", "NA", "XX"), stringsAsFactors = FALSE)
records2 <- records2 %>% filter(Island == "Daphne")
records2 <- records2[,c("meaningful.unique", "Species", "First.year", "Last.year")]
records2 <- records2 %>% filter(Species %in% c("fortis", "scandens"))
records2 <- records2 %>% rename(ind=meaningful.unique, start=First.year, end=Last.year)
records2 <- records2[!duplicated(records2$ind), ]  # remove duplicates

# merge them all and label the duplicated one
records <- rbind(records1, records2)
records[records$ind=="03Dap19726",][1,"ind"]="2:03Dap19726"
records <- records[!duplicated(records$ind), ]  # remove the last duplicate

# only include fortis and scandens
species_order <- c("fortis", "scandens")
records <- records %>% filter(Species %in% species_order)

# choose k=10 and read the admixture data
data <- read.table("finches.11.Q")
fam <- read.table("finches.fam")
data$ind <- fam$V2
data$Species <- records$Species[match(data$ind, records$ind)]
data <- data[!is.na(data$Species), ]  # remove individuals not in the records

# add the years and set last year as first year (if NA)
data$start <- records$start[match(data$ind, records$ind)]
data$end <- records$end[match(data$ind, records$ind)]
data$end[is.na(data$end)] <- data$start[is.na(data$end)]
data <- data[!is.na(data$start),]

# rename the columns based on ancestry
data <- data %>% rename(anc_fort = V9, anc_scan = V3, anc_ful = V4, anc_mag = V7, anc_mag2 = V6, anc_con = V8)
data <- data[, c("ind", "Species", "start", "end", "anc_fort", "anc_scan", "anc_ful", "anc_mag", "anc_mag2", "anc_con")]

# Function to summarize ancestry for each species
df_expanded <- data %>% rowwise() %>% mutate(year = list(start:end)) %>% unnest(year)
summarize_species <- function(species_name) {
  df_expanded %>% filter(Species == species_name) %>% group_by(year) %>% summarise(n_alive = n(), across(starts_with("anc_"),
    list(mean = ~mean(.x, na.rm = TRUE), sd = ~sd(.x, na.rm = TRUE)), .names = "{.col}_{.fn}")) %>% pivot_longer(cols = starts_with("anc_"),
    names_to = c("ancestry", ".value"), names_pattern = "(anc_.*)_(mean|sd)") %>% mutate(ancestry = recode(ancestry, anc_fort = "fortis",
    anc_scan = "scandens", anc_ful = "fuliginosa", anc_mag = "magnirostris", anc_mag2 = "magnirostris 2", anc_con = "conirostris"), Species = species_name)}
summary_all <- bind_rows(lapply(species_order, summarize_species))

# Plotting function
ancestry_colors <- c("fortis" = "dodgerblue2", "scandens" = "seagreen3", "fuliginosa" = "slategray2",
  "magnirostris" = "goldenrod2", "magnirostris 2" = "mediumpurple1", "conirostris" = "chocolate3")
plot_species <- function(species_name, show_y, show_legend) {
  df <- summary_all %>% filter(Species == species_name)
  p <- ggplot(df, aes(x = year, y = mean, color = ancestry)) +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd, color = ancestry, group = ancestry), width = 0.5, alpha = 0.8, linewidth = 0.5) +
    geom_line(linewidth = 1) + geom_point(aes(size = n_alive), show.legend = FALSE) + scale_color_manual(values = ancestry_colors) +
    scale_size_continuous(name = "Number Alive") + theme_minimal(base_size = 15) + scale_y_continuous() + coord_cartesian(ylim = c(0, 1)) +
    labs(title = paste("G.", species_name), x = NULL, color = "Ancestry") +
    theme(plot.title = element_text(hjust = 0.5, face = "italic"), legend.position = if (show_legend) "right" else "none")
  if (show_y) {
    p <- p + labs(y = "Mean Ancestry Proportion") + theme(plot.margin = margin(0.3, 0, 0, 0.3, unit = "cm"))}
  else {
    p <- p + labs(y = NULL) + theme(plot.margin = margin(0.3, 0, 0, 0, unit = "cm"), axis.text.y = element_blank(), axis.ticks.y = element_blank())}
  return(p)}
plots <- lapply(seq_along(species_order), function(i) {plot_species(species_name = species_order[i], show_y = (i == 1), show_legend = (i == length(species_order)))})

# Save to PNG
png(file = "ancestry_time_k11.png", width = 9000, height = 4000, res = 600)
wrap_plots(plots, nrow = 1, guides = "collect") & theme(legend.position = "right")
dev.off()