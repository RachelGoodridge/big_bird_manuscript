# switch to working in R and define inputs
library(tidyverse)
library(tidyr)
library(patchwork)
library(grid)

# get species + island for every training individual from the master spreadsheet
master <- read.csv("/local/storage/all_birds_highcov_lowcov/Finch_Database_V2.csv", header = TRUE, fill = TRUE, na.strings = c("", "NA", "XX"), stringsAsFactors = FALSE)
master <- master[!duplicated(master$meaningful.unique),]
master <- master[, c("meaningful.unique", "Species", "Island")]

# the founder is its own population (it's the individual that seeded the Big Bird lineage)
bbs_lines <- read.csv("/local/storage/BigBirds/BigBirds.csv")
bbs_lines <- bbs_lines[!duplicated(bbs_lines$meaningful.unique),]
master$lineage <- bbs_lines$lineage[match(master$meaningful.unique, bbs_lines$meaningful.unique)]
master <- master %>% mutate(Species = coalesce(lineage, Species)) %>% select(-lineage)

# clean up island names: Culpepper is the old name for Darwin island, and Wenman for Wolf
master <- master %>% mutate(Island = recode(Island, "Culpepper(Darwin)" = "Darwin", "Wenman(Wolf)" = "Wolf"))

# the master sheet lumps some island populations under a species name that's no longer
# used for that island: conirostris on Genovesa is actually propinqua, and difficilis on
# Genovesa/Darwin/Wolf are actually acutirostris/septentrionalis respectively
master <- master %>% mutate(Species = case_when(
  Species == "conirostris" & Island == "Genovesa" ~ "propinqua",
  Species == "difficilis" & Island == "Genovesa" ~ "acutirostris",
  Species == "difficilis" & Island %in% c("Darwin", "Wolf") ~ "septentrionalis",
  TRUE ~ Species
))
master$spp_isl <- paste(master$Species, master$Island, sep = "_")

# verified population list for the 156-individual training panel (species x island
# combinations actually present in keep_list.txt), left to right on the plot. the
# lone "bigbirds"-labeled individual is 04Dap5110, which the lineage coalesce above
# already relabels to "founder" -- so there's no separate bigbirds population here.
populations <- c("fortis_Daphne", "fortis_SantaCruz", "scandens_Daphne", "fuliginosa_Daphne", "fuliginosa_SantaCruz", "fuliginosa_Santiago",
  "magnirostris_Daphne", "magnirostris_Genovesa", "conirostris_Espanola", "founder_Daphne", "propinqua_Genovesa",
  "acutirostris_Genovesa", "septentrionalis_Darwin", "septentrionalis_Wolf",
  "difficilis_Fernandina", "difficilis_Pinta", "difficilis_Santiago")

# function to create each admixture subplot
make_admixture_plot <- function(k, show_x = FALSE, spacing = 1) {
  # read the data and add species information
  data <- read.table(paste0("finches.", k, ".Q"))
  fam <- read.table("finches.fam")
  data$ind <- fam$V2
  data$Species <- master$spp_isl[match(data$ind, master$meaningful.unique)]
  data <- data %>% filter(Species %in% populations) %>% mutate(Species = factor(Species, levels = populations))

  # figure out the highest mean ancestry for each population and sort individuals by it
  ancestry_cols <- grep("^V\\d+", names(data), value = TRUE)
  species_means <- data %>% group_by(Species) %>% summarise(across(all_of(ancestry_cols), \(x) mean(x, na.rm = TRUE)), .groups = "drop")
  species_dominant_ancestry <- species_means %>% pivot_longer(cols = all_of(ancestry_cols), names_to = "Ancestry", values_to = "MeanProp") %>%
    group_by(Species) %>% slice_max(order_by = MeanProp, n = 1, with_ties = FALSE) %>% select(Species, DominantAncestry = Ancestry)
  data <- left_join(data, species_dominant_ancestry, by = "Species")
  data <- data %>% rowwise() %>% mutate(DomValue = get(DominantAncestry)) %>% ungroup() %>% arrange(Species, desc(DomValue))

  # create space in between populations
  blank_row <- as_tibble(matrix(NA_real_, nrow = spacing, ncol = length(ancestry_cols)), .name_repair = "minimal")
  colnames(blank_row) <- ancestry_cols
  df_with_gaps <- data %>% group_split(Species) %>% purrr::map_dfr(~bind_rows(.x, blank_row)) %>%
    mutate(RowID = row_number()) %>% tidyr::fill(Species, .direction = "downup")
  df_long <- df_with_gaps %>% pivot_longer(cols = all_of(ancestry_cols), names_to = "Ancestry", values_to = "Proportion") %>% filter(!is.na(Proportion))

  # make the plot
  colors <- c("dodgerblue2", "seagreen3", "slategray1", "goldenrod2", "mediumpurple1", "chocolate3",
              "firebrick1", "darkorange", "forestgreen", "darkorchid4", "lightpink1", "maroon3")
  p <- ggplot(df_long, aes(x = RowID, y = Proportion, fill = Ancestry)) + geom_bar(stat = "identity", width = 1) + scale_fill_manual(values = colors[1:k]) +
    theme_minimal(base_size = 12) + labs(y = paste0("K = ", k)) + theme(axis.ticks.x = element_blank(), axis.text.y = element_blank(), panel.grid = element_blank(),
    axis.ticks.y = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(size = 12), legend.position = "none")

  if (show_x) {
    species_positions <- df_long %>% group_by(Species) %>% summarise(mid = mean(RowID, na.rm = TRUE), .groups = "drop") %>% filter(Species %in% populations)
    p <- p + theme(axis.text.x = element_text(size = 9, angle = 90, vjust = 0.5, hjust = 1), plot.margin = margin(-0.2, 0, 0, 0, unit = "cm"))
    p <- p + scale_x_continuous(breaks = species_positions$mid, labels = species_positions$Species, expand = expansion(add = 1))
  } else {
    p <- p + theme(axis.text.x = element_blank(), plot.margin = margin(-0.2, 0, -0.2, 0, unit = "cm")) + scale_x_continuous(expand = expansion(add = 1))
  }
  return(p)
}

# CV-error curve, for reference alongside the visual comparison below.
# NOTE: with this many known reference populations, the CV-minimum K often
# under-resolves visually distinct groups (e.g. last time CV bottomed out around
# K ~ 6 but K = 10-11 was needed to keep every reference population visually
# distinct) -- use this plot together with subset_admix.png, not instead of it.
cv <- read.table("finches.cv.error", col.names = c("K", "CV_error"))
png(file = "cv_error.png", width = 1200, height = 900, res = 200)
print(ggplot(cv, aes(K, CV_error)) + geom_line() + geom_point() + theme_light() +
      scale_x_continuous(breaks = cv$K) + labs(x = "k populations", y = "cv error"))
dev.off()

# multi-panel comparison across k = 2..12
png(file = "subset_admix.png", width = 2000, height = 1800, res = 200)
plots <- map(2:12, ~make_admixture_plot(.x, show_x = (.x == 12)))
grid.newpage()
grid.draw(do.call(rbind, map(plots, ggplotGrob)))
dev.off()
