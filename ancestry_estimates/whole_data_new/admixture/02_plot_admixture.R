# switch to working in R and define inputs
library(tidyverse)
library(tidyr)
library(patchwork)
library(grid)
library(ggnewscale)

# ── set these once you've picked K and run 01_run_projection.sh <K> ─────────
# whole_data_new/plink/bird_selection.txt is training's keep_list.txt union every
# Daphne bird, so this projection's Q file already contains founder/conirostris/
# fuliginosa/magnirostris/etc. directly -- unlike the old admixture_with_training.R,
# there's no need to separately pull them in from a training-set Q file anymore.
project_dir <- "/SSD/reg259/big_bird_admixture/whole_data_new/admixture"
k <- 9  # <- set to the K you chose from training_new/admixture/finches.cv.error + subset_admix.png

# the ancestry cluster each named population is dominated by, for this K. to find these,
# run once with any placeholder mapping, then in R:
#   data <- read.table(file.path(project_dir, paste0("finches.", k, ".Q")))
#   fam <- read.table(file.path(project_dir, "finches.fam")); data$ind <- fam$V2
#   data$Species <- all_records$Species[match(data$ind, all_records$meaningful.unique)]
#   data %>% group_by(Species) %>% summarise(across(starts_with("V"), \(x) mean(x, na.rm = TRUE))) %>% print(width = Inf)
# and read off which V-column has the highest mean for each named population.
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

# 01_whole_data_bed.sh now projects every sample in the raw VCF (not just the curated
# 4330-sample bird_selection.txt), which includes e.g. fortis from Santa Cruz and
# conirostris from Genovesa. Restrict this plot to exactly the populations of interest.
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

# italic type reads visually larger/heavier than upright type at the same point size (a standard
# optical effect, not a font/size bug), so species-name labels are drawn a bit smaller to compensate
italic_size_factor <- 0.9

# function to create each admixture subplot
make_admixture_plot <- function(k, show_x = FALSE, spacing = 2, year_cutoff = 1987, label_size = 20) {
  # read the whole-data individuals and add species + first-year information
  data <- read.table(file.path(project_dir, paste0("finches.", k, ".Q")))
  fam <- read.table(file.path(project_dir, "finches.fam"))
  data$ind <- fam$V2
  data <- left_join(data, all_records, by = c("ind" = "meaningful.unique"))

  # order populations left to right. difficilis/magnirostris/fuliginosa reference birds
  # are also present in this Q file if you want to add them to this list
  populations <- c("scandens", "fortis", "conirostris", "founder", "A line", "B line", "Big Birds")
  data <- data %>% filter(Species %in% populations)
  data <- data %>% mutate(Species = factor(Species, levels = populations))

  # fortis and scandens have far more individuals than the other populations, so for legibility
  # restrict them to birds first seen before year_cutoff; the rest are left as-is
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

  # flag the single founder individual with an asterisk, since it's a conirostris hybrid
  # sitting flush against the conirostris panel and doesn't get its own x-axis label
  founder_mid <- mean(df_long$RowID[df_long$Species == "founder"], na.rm = TRUE)
  p <- p + annotate("text", x = founder_mid, y = 1.05, label = "*", size = label_size / 2.845276, fontface = "bold") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12)))

  # add species labels only to the bottom plot. drawn as a manual geom_text layer (rather than the
  # usual axis.text.x theme) so italic species names can be sized down relative to the upright
  # population labels -- ggplot2's theme system only allows one font size per axis, not per label
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
# with exactly 4 entries (the 3 highlighted ancestries + "other") -- this avoids showing
# every one of the K ancestry columns in the real legend, since only 3 of them have a
# distinct (non-gray) color. its bars are invisible (alpha = 0); only the legend is kept,
# left-justified via legend.justification so it doesn't sit over the founder's asterisk.
#
# the species entries and the "other ancestries" entry are drawn via two independent
# guide_legend()s (fill for species, colour for "other") combined into one row via
# legend.box = "horizontal" -- this is the only way to give the species entries their own
# smaller guide theme (italic_size_factor, to compensate for italic's optical size effect)
# while "other ancestries" keeps the normal legend text size; ggplot2's theme system only
# supports one legend.text size per guide, not per key
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
