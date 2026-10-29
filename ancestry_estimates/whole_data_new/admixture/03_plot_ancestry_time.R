# switch to working in R and define inputs
library(tidyverse)
library(patchwork)

# move to the directory of this script
args <- commandArgs(trailingOnly = FALSE)
script_path <- gsub("--file=", "", args[grep("--file=", args)])
if (length(script_path) > 0) setwd(dirname(script_path))

# same as whole_data_old/ancestry_time.R, but for the new K=9 projection and only fortis + scandens ancestry
k <- 9
print(paste("plotting K =", k))

# get all the Daphne fortis/scandens records from the master spreadsheet
species_order <- c("fortis", "scandens")
records <- read.csv("/local/storage/all_birds_highcov_lowcov/Finch_Database_V2.csv", header = TRUE, fill = TRUE, na.strings = c("", "NA", "XX"), stringsAsFactors = FALSE)
records <- records[!duplicated(records$meaningful.unique), ]  # remove duplicates (probably they are all re-captures)

# update the newest lineage classifications (as in 02_plot_admixture.R), so A line / B line / Big Birds aren't counted as fortis
bbs_lines <- read.csv("/local/storage/BigBirds/BigBirds.csv")
bbs_lines <- bbs_lines[!duplicated(bbs_lines$meaningful.unique),]
records$Species <- coalesce(bbs_lines$lineage[match(records$meaningful.unique, bbs_lines$meaningful.unique)], records$Species)

# correct the species for an individual (before keeping only fortis / scandens)
records <- records |> mutate(Species = if_else(meaningful.unique == "01Dap21269", "magnirostris", Species))
records <- records %>% filter(Island == "Daphne", Species %in% species_order)

# first year alive: the earliest of year born, first year seen, and earliest year recorded (last year seen if none are known), as in 02_plot_admixture.R
records <- records %>% mutate(start = coalesce(pmin(Year.born, First.year, First.year.min, na.rm = TRUE), Last.year), end = coalesce(Last.year, start)) %>%
  select(ind = meaningful.unique, Species, start, end)

# read the admixture data and add the species + years
data <- read.table(paste0("finches.", k, ".Q"))
data$ind <- read.table("finches.fam")$V2
data <- inner_join(data, records, by = "ind") %>% filter(!is.na(start))

# rename the columns based on ancestry (see species_ancestry in 02_plot_admixture.R)
data <- data %>% rename(anc_fort = V1, anc_scan = V6) %>% select(ind, Species, start, end, anc_fort, anc_scan)

# Function to summarize ancestry for each species
df_expanded <- data %>% rowwise() %>% mutate(year = list(start:end)) %>% unnest(year)
summarize_species <- function(species_name) {
  df_expanded %>% filter(Species == species_name) %>% group_by(year) %>% summarise(n_alive = n(), across(starts_with("anc_"),
    list(mean = ~mean(.x, na.rm = TRUE), sd = ~sd(.x, na.rm = TRUE)), .names = "{.col}_{.fn}")) %>% pivot_longer(cols = starts_with("anc_"),
    names_to = c("ancestry", ".value"), names_pattern = "(anc_.*)_(mean|sd)") %>% mutate(ancestry = recode(ancestry, anc_fort = "fortis",
    anc_scan = "scandens"), Species = species_name)}
summary_all <- bind_rows(lapply(species_order, summarize_species))

# Plotting function
ancestry_colors <- c("fortis" = "dodgerblue2", "scandens" = "seagreen3")
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
png(file = paste0("ancestry_time_k", k, ".png"), width = 9000, height = 4000, res = 600)
print(wrap_plots(plots, nrow = 1, guides = "collect") & theme(legend.position = "right"))
dev.off()
