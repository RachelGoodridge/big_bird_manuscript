# switch to working in R and define inputs
library(tidyverse)
library(tidyr)
library(patchwork)
library(grid)
library(ggnewscale)

# ── set these once you've picked K and run 01_run_projection.sh <K> ─────────
project_dir <- "/SSD/reg259/big_bird_admixture/whole_data_new/admixture"
k <- 9

# the ancestry cluster each named population is dominated by, for this K.
species_ancestry <- c(fortis = "V1", scandens = "V6", conirostris = "V2")
highlight_colors <- c(fortis = "#1C86EE", scandens = "#43CD80", conirostris = "#ff3030dd")

# get all the records from the master spreadsheet
all_records <- read.csv("/local/storage/all_birds_highcov_lowcov/Finch_Database_V2.csv", header = TRUE, fill = TRUE, na.strings = c("", "NA", "XX"), stringsAsFactors = FALSE)
all_records <- all_records[!duplicated(all_records$meaningful.unique),]  # remove duplicates (probably they are all re-captures)
all_records <- all_records[, c("meaningful.unique", "Species", "Island", "First.year", "First.year.min", "Year.born", "Last.year")]

# fill in missing years with the minimum year, year born, or last year, in that order of preference
all_records <- all_records %>% mutate(First.year = coalesce(First.year, First.year.min, Year.born, Last.year),
                                       Last.year = coalesce(Last.year, First.year))

# update the newest lineage classifications
bbs_lines <- read.csv("/local/storage/BigBirds/BigBirds.csv")
bbs_lines <- bbs_lines[!duplicated(bbs_lines$meaningful.unique),]
all_records$lineage <- bbs_lines$lineage[match(all_records$meaningful.unique, bbs_lines$meaningful.unique)]
all_records <- all_records %>% mutate(lineage = coalesce(lineage, Species))
all_records <- all_records %>% select(-Species) %>% rename(Species = lineage)

# 01_whole_data_bed.sh now projects every sample in the raw VCF
# Restrict this plot to exactly the populations of interest.
allowed_combos <- tribble(
  ~Species, ~Island,
  "fortis", "Daphne",
  "scandens", "Daphne",
  "A line", "Daphne",
  "B line", "Daphne",
  "Big Birds", "Daphne",
  "founder", "Daphne",
  "conirostris", "Espanola")
all_records <- semi_join(all_records, allowed_combos, by = c("Species", "Island")) %>% select(-Island)

# abbreviated x-axis labels; founder gets no label at all (it's a single individual, technically a conirostris)
display_labels <- c(scandens = "scandens", fortis = "fortis", conirostris = "conirostris",
                     founder = NA, "A line" = "A line", "B line" = "B line", "Big Birds" = "Big Birds")
italic_size_factor <- 0.9

# function to create each admixture subplot
make_admixture_plot <- function(k, show_x = FALSE, spacing = 2, year_cutoff = 1987, label_size = 20) {
  # read the whole-data individuals and add species + first-year information
  data <- read.table(file.path(project_dir, paste0("finches.", k, ".Q")))
  fam <- read.table(file.path(project_dir, "finches.fam"))
  data$ind <- fam$V2
  data <- left_join(data, all_records, by = c("ind" = "meaningful.unique"))

  # order populations left to right
  populations <- c("scandens", "fortis", "conirostris", "founder", "A line", "B line", "Big Birds")
  data <- data %>% filter(Species %in% populations)
  data <- data %>% mutate(Species = factor(Species, levels = populations))

  # restrict fortis and scandens to birds first seen before year_cutoff
  data <- data %>% filter(!(Species %in% c("fortis", "scandens")) | First.year < year_cutoff)

  # the gap after each population's panel; conirostris sits almost flush against founder
  # since founder is itself a conirostris hybrid, everything else uses the standard spacing
  gap_after <- setNames(rep(spacing, length(populations)), populations)
  gap_after[["conirostris"]] <- 1

  # figure out the highest mean ancestry for each species (used only to order individuals within each panel)
  ancestry_cols <- grep("^V\\d+", names(data), value = TRUE)
  species_means <- data %>% group_by(Species) %>% summarise(across(all_of(ancestry_cols), \(x) mean(x, na.rm = TRUE)), .groups = "drop")
  species_dominant_ancestry <- species_means %>% pivot_longer(cols = all_of(ancestry_cols), names_to = "Ancestry", values_to = "MeanProp") %>%
    group_by(Species) %>% slice_max(order_by = MeanProp, n = 1, with_ties = FALSE) %>% select(Species, DominantAncestry = Ancestry)

  # add the dominant ancestry column and sort
  data <- left_join(data, species_dominant_ancestry, by = "Species")
  data <- data %>% rowwise() %>% mutate(DomValue = get(DominantAncestry)) %>% ungroup() %>% arrange(Species, desc(DomValue))

  # create space in between species (variable-size gaps, per gap_after)
  df_with_gaps <- data %>% group_split(Species) %>% purrr::map_dfr(function(grp) {
    n_gap <- gap_after[[as.character(grp$Species[1])]]
    blank <- as_tibble(matrix(NA_real_, nrow = n_gap, ncol = length(ancestry_cols)), .name_repair = "minimal")
    colnames(blank) <- ancestry_cols
    bind_rows(grp, blank)
  }) %>% mutate(RowID = row_number()) %>% tidyr::fill(Species, .direction = "downup")

  # rank each ancestry cluster within each species by its mean proportion, so every panel can be
  # stacked bottom-to-top from that species' own dominant ancestry down to its rarest one; every
  # non-highlighted ("gray") ancestry cluster is pushed to the top of the stack regardless of its size
  species_rank <- species_means %>% pivot_longer(cols = all_of(ancestry_cols), names_to = "Ancestry", values_to = "MeanProp") %>%
    group_by(Species) %>% mutate(is_gray = !(Ancestry %in% unname(species_ancestry)), Rank = rank(is_gray * 1e6 - MeanProp, ties.method = "first")) %>%
    ungroup() %>% select(Species, Ancestry, Rank)

  df_long <- df_with_gaps %>% pivot_longer(cols = all_of(ancestry_cols), names_to = "Ancestry", values_to = "Proportion") %>% filter(!is.na(Proportion))
  df_long <- df_long %>% left_join(species_rank, by = c("Species", "Ancestry"))
  df_long <- df_long %>% arrange(RowID, Rank) %>% group_by(RowID) %>% mutate(ymax = cumsum(Proportion), ymin = lag(ymax, default = 0)) %>% ungroup()

  # only fortis/scandens/conirostris ancestry get a distinct color; every other ancestry cluster is lumped into gray
  colors <- setNames(rep("grey80", length(ancestry_cols)), ancestry_cols)
  colors[unname(species_ancestry)] <- highlight_colors[names(species_ancestry)]

  p <- ggplot(df_long, aes(xmin = RowID - 0.5, xmax = RowID + 0.5, ymin = ymin, ymax = ymax, fill = Ancestry)) + geom_rect() + scale_fill_manual(values = colors) +
    theme_minimal(base_size = 24) + theme(axis.ticks.x = element_blank(), axis.text.y = element_blank(), panel.grid = element_blank(),
    axis.ticks.y = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "none")

  # flag the single founder individual with an asterisk
  founder_mid <- mean(df_long$RowID[df_long$Species == "founder"], na.rm = TRUE)
  p <- p + annotate("text", x = founder_mid, y = 1.05, label = "*", size = label_size / 2.845276, fontface = "bold") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12)))

  # add species labels only to the bottom plot
  if (show_x) {
    species_positions <- df_long %>% group_by(Species) %>% summarise(mid = mean(RowID, na.rm = TRUE), .groups = "drop") %>%
      filter(Species %in% populations) %>% mutate(label = display_labels[as.character(Species)]) %>% filter(!is.na(label)) %>%
      mutate(is_species_name = as.character(Species) %in% names(species_ancestry),
             fontface = ifelse(is_species_name, "italic", "plain"),
             text_size = ifelse(is_species_name, label_size * italic_size_factor, label_size))
    p <- p + theme(axis.text.x = element_blank(), plot.margin = margin(-0.2, 0.3, 0.7, 0.3, unit = "cm")) +
      coord_cartesian(clip = "off") +
      scale_x_continuous(breaks = NULL, expand = expansion(add = c(0, 5))) +
      geom_text(data = species_positions, aes(x = mid, y = -0.12, label = label, fontface = fontface),
                size = species_positions$text_size / 2.845276, colour = "#4D4D4DFF", inherit.aes = FALSE)}
  else {
    p <- p + theme(axis.text.x = element_blank(), plot.margin = margin(-0.2, 0.3, -0.2, 0.3, unit = "cm")) + scale_x_continuous(expand = expansion(add = c(0, 5)))}

  return(p)
}

# build a small standalone plot whose only purpose is to generate a horizontal legend
legend_text_size <- calc_element("legend.text", theme_void(base_size = 24))$size
species_df <- data.frame(label = factor(names(species_ancestry), levels = names(species_ancestry)), x = 1, y = 1)
other_df <- data.frame(label = factor("other ancestries"), x = 1, y = 1)
legend_plot <- ggplot() +
  geom_col(data = species_df, aes(x = x, y = y, fill = label), alpha = 0) +
  scale_fill_manual(values = highlight_colors, name = NULL, labels = parse(text = paste0("italic(", names(species_ancestry), ")")),
                     guide = guide_legend(order = 1, override.aes = list(alpha = 1),
                                           theme = theme(legend.text = element_text(size = legend_text_size * italic_size_factor)))) +
  ggnewscale::new_scale_colour() +
  geom_point(data = other_df, aes(x = x, y = y, colour = label), alpha = 0) +
  scale_colour_manual(values = c("other ancestries" = "grey80"), name = NULL,
                       guide = guide_legend(order = 2, override.aes = list(alpha = 1, shape = 15, size = 6),
                                             theme = theme(legend.text = element_text(size = legend_text_size)))) +
  theme_void(base_size = 24) + theme(legend.position = "bottom", legend.direction = "horizontal", legend.box = "horizontal", legend.justification = "left")

# create and save the plot
stopifnot("set k and species_ancestry above before running" = !is.na(k), !any(grepl("V\\?", species_ancestry)))
png(file = "admix_proj_whole_data.png", width = 9000, height = 850, res = 350)
p <- make_admixture_plot(k, TRUE)
print(legend_plot / p + plot_layout(heights = c(1, 7.5)))
dev.off()
