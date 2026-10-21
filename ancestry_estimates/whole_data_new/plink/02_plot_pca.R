# switch to working in R and define inputs
library(tidyverse)

# get species for every whole-data individual from the master spreadsheet
master <- read.csv("/local/storage/all_birds_highcov_lowcov/Finch_Database_V2.csv", header = TRUE, fill = TRUE, na.strings = c("", "NA", "XX"), stringsAsFactors = FALSE)
master <- master[!duplicated(master$meaningful.unique),]
master <- master[, c("meaningful.unique", "Species", "Island")]

# the Big Bird pedigree (A line / B line / Big Birds / founder) overrides raw Species where tracked
bbs_lines <- read.csv("/local/storage/BigBirds/BigBirds.csv")
bbs_lines <- bbs_lines[!duplicated(bbs_lines$meaningful.unique),]
master$lineage <- bbs_lines$lineage[match(master$meaningful.unique, bbs_lines$meaningful.unique)]
master <- master %>% mutate(Species = coalesce(lineage, Species)) %>% select(-lineage)

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
master <- semi_join(master, allowed_combos, by = c("Species", "Island")) %>% select(-Island)

# read the PCA output from 01_whole_data_bed.sh (plink2 always writes FID then IID first)
pca <- read_table("finches.eigenvec", col_names = TRUE)
names(pca)[1] <- "IID"
eigenval <- scan("finches.eigenval")
pve <- data.frame(PC = seq_along(eigenval), pve = eigenval / sum(eigenval) * 100)

pca$Species <- master$Species[match(pca$IID, master$meaningful.unique)]
pca <- pca %>% filter(!is.na(Species))

# make the pca plot -- restricted (via allowed_combos above) to fortis/scandens/A line/
# B line/Big Birds/founder on Daphne, and conirostris on Espanola. colors match the
# standard project-wide palette (color_map_hex)
ancestry_colors <- c(fortis = "#1C86EE", scandens = "#43CD80", conirostris = "#8B4513", founder = "#FF3030",
                     "A line" = "#FF8C00", "B line" = "#FFAEB9", "Big Birds" = "#800080")
pca <- pca[order(pca$Species == "founder"),] # draw founder last so its on top
p <- ggplot(pca, aes(PC1, PC2, col = Species)) + geom_point(size = 1.5) + theme_light() + scale_color_manual(values = ancestry_colors) +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))
ggsave("whole_data_pca.png", plot = p, height = 5, width = 7, units = "in", dpi = 300)