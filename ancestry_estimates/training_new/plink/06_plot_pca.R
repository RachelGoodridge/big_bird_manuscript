# switch to working in R and define inputs
library(tidyverse)

# get species + island for every training individual from the master spreadsheet
master <- read.csv("/local/storage/all_birds_highcov_lowcov/Finch_Database_V2.csv", header = TRUE, fill = TRUE, na.strings = c("", "NA", "XX"), stringsAsFactors = FALSE)
master <- master[!duplicated(master$meaningful.unique),]
master <- master[, c("meaningful.unique", "Species", "Island")]

# the founder is its own population (it's the individual that seeded the Big Bird lineage)
bbs_lines <- read.csv("/local/storage/BigBirds/BigBirds.csv")
bbs_lines <- bbs_lines[!duplicated(bbs_lines$meaningful.unique),]
master$lineage <- bbs_lines$lineage[match(master$meaningful.unique, bbs_lines$meaningful.unique)]
master <- master %>% mutate(Species = coalesce(lineage, Species)) %>% select(-lineage)

# read the PCA output from 05_training_bed.sh (plink2 writes just IID here, no separate FID column)
pca <- read_table("finches.eigenvec", col_names = TRUE)
names(pca)[1] <- "IID"
eigenval <- scan("finches.eigenval")
pve <- data.frame(PC = seq_along(eigenval), pve = eigenval / sum(eigenval) * 100)

pca$Species <- master$Species[match(pca$IID, master$meaningful.unique)]
pca$Island <- master$Island[match(pca$IID, master$meaningful.unique)]

# conirostris on Genovesa is actually propinqua, not in the palette so it drops out below
pca$Species <- ifelse(pca$Species == "conirostris" & pca$Island == "Genovesa", "propinqua", pca$Species)

# make the pca plot -- standard project-wide palette; anything not in it (just
# difficilis, in this training panel) falls back to "unknown" gray
color_map_hex <- c(fortis = "#1C86EE", scandens = "#43CD80", fuliginosa = "#C6E2FF", hybrid = "#AB82FF", "A line" = "#FF8C00", founder = "#FF3030",
                   magnirostris = "#EEB422", conirostris = "#8B4513", "B line" = "#FFAEB9", unknown = "#808080", "Big Birds" = "#800080")
pca$Species <- ifelse(pca$Species %in% names(color_map_hex), pca$Species, "unknown")
pca <- pca[pca$Species != "unknown",]
pca <- pca[order(pca$Species == "founder"),] # draw founder last so its on top
p <- ggplot(pca, aes(PC1, PC2, col = Species)) + geom_point(size = 1.5) + theme_light() + scale_color_manual(values = color_map_hex) +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))
ggsave("training_pca.png", plot = p, height = 5, width = 7, units = "in", dpi = 300)