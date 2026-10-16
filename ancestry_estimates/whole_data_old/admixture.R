# switch to working in R and define inputs
library(tidyverse)
library(tidyr)
library(patchwork)
library(grid)

# get all the records from the master spreadsheet
records <- read.csv("/local/storage/all_birds_highcov_lowcov/Finch_Database_V2.csv", header = TRUE, fill = TRUE, na.strings = c("", "NA", "XX"), stringsAsFactors = FALSE)
records <- records[!duplicated(records$meaningful.unique),]  # remove duplicates (probably they are all re-captures)
records <- records %>% filter(Island == "Daphne")
records <- records[, c("meaningful.unique", "Species")]

# update the newest lineage classifications
bbs_lines <- read.csv("/local/storage/BigBirds/BigBirds.csv")
bbs_lines <- bbs_lines[!duplicated(bbs_lines$meaningful.unique),]
records$lineage <- bbs_lines$lineage[match(records$meaningful.unique, bbs_lines$meaningful.unique)]
records <- records %>% mutate(lineage = coalesce(lineage, Species))
records <- records %>% select(-Species) %>% rename(Species = lineage)

# function to create each admixture subplot
make_admixture_plot <- function(k, show_x = FALSE, spacing = 10) {
  # read the data and add species information
  data <- read.table(paste0("finches.", k, ".Q"))
  fam <- read.table("finches.fam")
  data$ind <- fam$V2
  data$Species <- records$Species[match(data$ind, records$meaningful.unique)]

  # first create an order by species (and filter out non-Daphne species)
  populations <- c("fortis", "scandens", "A line", "B line", "Big Birds")
  data <- data %>% filter(Species %in% populations)
  data <- data %>% mutate(Species = factor(Species, levels = populations))

  # figure out the highest mean ancestry for each species
  ancestry_cols <- grep("^V\\d+", names(data), value = TRUE)
  species_means <- data %>% group_by(Species) %>% summarise(across(all_of(ancestry_cols), \(x) mean(x, na.rm = TRUE)), .groups = "drop")
  species_dominant_ancestry <- species_means %>% pivot_longer(cols = all_of(ancestry_cols), names_to = "Ancestry", values_to = "MeanProp") %>%
    group_by(Species) %>% slice_max(order_by = MeanProp, n = 1, with_ties = FALSE) %>% select(Species, DominantAncestry = Ancestry)
  
  # # write the averages to a file
  # output_file <- "ancestry_proportions.txt"
  # if (!file.exists(output_file)) {
  #   cat("Average ancestry proportions for K\n", file = output_file)}
  # cat("\nAverage ancestry proportions for K =", k, "\n", file = output_file, append = TRUE)
  # write.table(species_means %>% column_to_rownames("Species") %>% round(3),
  #   file = output_file, append = TRUE, sep = "\t", row.names = TRUE, col.names = !file.exists(output_file))

  # add the dominant ancestry column and sort
  data <- left_join(data, species_dominant_ancestry, by = "Species")
  data <- data %>% rowwise() %>% mutate(DomValue = get(DominantAncestry)) %>% ungroup() %>% arrange(Species, desc(DomValue))

  # create space in between species
  blank_row <- as_tibble(matrix(NA_real_, nrow = spacing, ncol = length(ancestry_cols)), .name_repair = "minimal")
  colnames(blank_row) <- ancestry_cols

  # combine data with blank rows and re-format for plotting
  df_with_gaps <- data %>% group_split(Species) %>% purrr::map_dfr(~bind_rows(.x, blank_row)) %>%
    mutate(RowID = row_number()) %>% tidyr::fill(Species, .direction = "downup")
  df_long <- df_with_gaps %>% pivot_longer(cols = all_of(ancestry_cols), names_to = "Ancestry", values_to = "Proportion") %>% filter(!is.na(Proportion))

  # make the plot
  colors <- c("dodgerblue2", "seagreen3", "slategray1", "goldenrod2", "mediumpurple1", "chocolate3",
              "firebrick1", "darkorange", "forestgreen", "darkorchid4", "lightpink1", "maroon3")
  p <- ggplot(df_long, aes(x = RowID, y = Proportion, fill = Ancestry)) + geom_bar(stat = "identity", width = 1) + scale_fill_manual(values = colors[1:k]) +
    theme_minimal(base_size = 24) + labs(y = paste0("K = ", k)) + theme(axis.ticks.x = element_blank(), axis.text.y = element_blank(), panel.grid = element_blank(),
    axis.ticks.y = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(size = 24), legend.position = "none")

  # add species labels only to the bottom plot
  if (show_x) {
    species_positions <- df_long %>% group_by(Species) %>% summarise(mid = mean(RowID, na.rm = TRUE), .groups = "drop") %>% filter(Species %in% populations)
    p <- p + theme(axis.text.x = element_text(size = 24), plot.margin = margin(-0.2, 0, 0, 0, unit = "cm"))
    p <- p + scale_x_continuous(breaks = species_positions$mid, labels = species_positions$Species, expand = expansion(add = 1))}
  else {
    p <- p + theme(axis.text.x = element_blank(), plot.margin = margin(-0.2, 0, -0.2, 0, unit = "cm")) + scale_x_continuous(expand = expansion(add = 1))}

  return(p)
}

# # create and save the plot for all the k values
# png(file = "all_Daphne_admix_proj.png", width = 3000, height = 1500, res = 200)
# p1 <- make_admixture_plot(10)
# p2 <- make_admixture_plot(11, TRUE)
# grid.newpage()
# grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2)))
# dev.off()

# create and save the plot for all the k values
png(file = "all_Daphne_admix_proj_MS.png", width = 3000, height = 750, res = 350)
make_admixture_plot(11, TRUE)
dev.off()