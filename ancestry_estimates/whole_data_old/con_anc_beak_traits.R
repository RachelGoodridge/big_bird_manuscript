# switch to working in R and define inputs
library(tidyverse)
library(tidyr)
library(patchwork)
library(grid)

# get all the records from the master spreadsheet
records <- read.csv("/local/storage/all_birds_highcov_lowcov/Finch_Database_V2.csv", header = TRUE, fill = TRUE, na.strings = c("", "NA", "XX"), stringsAsFactors = FALSE)
records <- records %>% filter(Island == "Daphne")
records <- records[,c("meaningful.unique", "Species", "bill.depth", "bill.length", "bill.width", "First.year", "Last.year")]
records <- records %>% filter(Species %in% c("fortis", "scandens", "bigbirds"))
records <- records %>% rename(ind=meaningful.unique, depth=bill.depth, length=bill.length, width=bill.width, start=First.year, end=Last.year)
records <- records %>% mutate(end = coalesce(end, start))  # if no end year, use start year
records <- records[!duplicated(records$ind),]  # remove duplicates (probably they are all re-captures)

# --- Define a function to read and prepare ancestry data ---
prepare_ancestry_data <- function(k_value, con, fort, scan) {
  # read the anccestry data
  data1 <- read.table(paste0("finches.", k_value, ".Q"))
  fam <- read.table("finches.fam")
  data1$ind <- fam$V2
  
  # get the big bird data
  data2 <- read.table(paste0("/SSD/reg259/expanded_big_bird_analysis/admixture_project/finches.", k_value, ".Q"))
  fam2 <- read.table("/SSD/reg259/expanded_big_bird_analysis/admixture_project/finches.fam")
  data2$ind <- fam2$V2

  # merge the two datasets
  bbs = records[records$Species == "bigbirds", "ind"]
  data2 <- data2[data2$ind %in% bbs, ]
  data <- rbind(data1, data2)

  # add species and beak information
  data$Species <- records$Species[match(data$ind, records$ind)]
  data <- data %>% filter(Species %in% c("fortis", "scandens", "bigbirds"))
  data$depth <- records$depth[match(data$ind, records$ind)]
  data$length <- records$length[match(data$ind, records$ind)]
  data$width <- records$width[match(data$ind, records$ind)]
  data$sum <- data$depth + data$length + data$width
  data$year <- records$start[match(data$ind, records$ind)]
  data$pointy <- data$length / data$depth

  # label the founder
  data[data$ind == "04Dap5110", "Species"] <- "founder"
  
  # calculate conirostris ancestry as a function of "main" ancestry
  data <- data %>% rename(anc_con = all_of(con), anc_fort = all_of(fort), anc_scan = all_of(scan))
  data <- data %>% mutate(anc_main = ifelse(Species == "scandens", anc_scan, anc_fort))
  data$con_modified <- data$anc_con / (1 - data$anc_main)
  data <- data[, c("ind", "Species", "depth", "length", "width", "sum", "anc_con", "con_modified", "anc_fort", "year", "pointy")]

  return(data)
}

# --- Prepare summaries for k=10 and k=11 ---
summary_k10 <- prepare_ancestry_data(k_value = 10, con = "V9", fort = "V10", scan = "V5")
summary_k11 <- prepare_ancestry_data(k_value = 11, con = "V8", fort = "V9", scan = "V3")

# --- Calculate correlation coefficients ---
k11_fort <- summary_k11 %>% filter(Species == "fortis")
k11_scan <- summary_k11 %>% filter(Species == "scandens")
cor.test(k11_fort$con_modified, k11_fort$sum, method = "pearson")
cor.test(k11_scan$con_modified, k11_scan$sum, method = "pearson")
cor.test(k11_fort$anc_con, k11_fort$sum, method = "pearson")
cor.test(k11_scan$anc_con, k11_scan$sum, method = "pearson")

# --- Define a function to plot the data ---
plot_con_ancestry <- function(summary_df, title_text, show_y, show_legend) {
  species_colors <- c("scandens" = "seagreen3", "fortis" = "dodgerblue2", "bigbirds" = "lightpink1", "founder" = "gray")
  
  p <- ggplot(summary_df, aes(x = anc_con, y = sum, color = Species)) + scale_color_manual(values = species_colors, labels = scales::parse_format()) +
    theme_minimal(base_size = 15) + labs(title = title_text, x = "Conirostris Ancestry", color = "Species") + geom_point() +
    theme(plot.title = element_text(hjust = 0.5, face = "plain"), legend.position = if (show_legend) "right" else "none", axis.text.y = if (!show_y) element_blank() else NULL,
      axis.ticks.y = if (!show_y) element_blank() else NULL, axis.title.y = if (!show_y) element_blank() else element_text(margin = margin(r = 10)))
  
  if (show_y) {p <- p + labs(y = "Sum of Beak Dimensions")}
  return(p)
}

# --- Create plots with control over legend and y-axis ---
plot_k10 <- plot_con_ancestry(summary_k10, title_text = "K=10", show_y = TRUE, show_legend = FALSE)
plot_k11 <- plot_con_ancestry(summary_k11, title_text = "K=11", show_y = FALSE, show_legend = TRUE)

# --- Save combined plot ---
png("con_anc_beak_traits.png", width = 5800, height = 2400, res = 400)
plot_k10 + plot_k11 + plot_layout(axes = "collect")
dev.off()

# --- individuals with >20% conirostris ancestry ---
png("year_con_anc.png", width = 2600, height = 1800, res = 400)
species_colors <- c("scandens" = "seagreen3", "fortis" = "dodgerblue2", "bigbirds" = "lightpink1", "founder" = "gray")
sub1 <- summary_k11[summary_k11$anc_con > 0.1, c("anc_con", "Species", "year")]
sub1 <- sub1 %>% filter(Species %in% c("fortis", "scandens", "bigbirds", "founder"))
ggplot(sub1, aes(x = year, y = anc_con, color = Species)) + scale_color_manual(values = species_colors, labels = scales::parse_format()) + theme_minimal(base_size = 15) +
  labs(title = "K=11", x = "First Year", y = "Conirostris Ancestry", color = "Species") + geom_point() + theme(plot.title = element_text(hjust = 0.5), legend.position = "right")
dev.off()


# --- anc_con vs anc_fort, shaded by pointy (scandens only) ---
k11_scan_pointy <- k11_scan %>% filter(!is.na(pointy))
p_k11_anc_scatter <- ggplot(k11_scan_pointy, aes(x = anc_con, y = anc_fort, color = pointy)) + geom_point(size = 2) + scale_color_viridis_c(option = "plasma") +
  theme_minimal(base_size = 15) + labs(title = "beak pointedness of scandens by ancestry", x = "conirostris ancestry", y = "fortis ancestry",
  color = "pointedness\n(length / depth)") + theme(plot.title = element_text(hjust = 0.5))

# save the plot
png("scandens_anc_con_vs_fort_pointy.png", width = 3200, height = 2000, res = 400)
p_k11_anc_scatter
dev.off()

# --- anc_con vs anc_fort, shaded by beak sum (scandens only) ---
k11_scan_sum <- k11_scan %>% filter(!is.na(sum))
p_k11_anc_scatter <- ggplot(k11_scan_sum, aes(x = anc_con, y = anc_fort, color = sum)) + geom_point(size = 2) + scale_color_viridis_c(option = "plasma") +
  theme_minimal(base_size = 15) + labs(title = "beak size of scandens by ancestry", x = "conirostris ancestry", y = "fortis ancestry",
  color = "beak size") + theme(plot.title = element_text(hjust = 0.5))

# save the plot
png("scandens_anc_con_vs_fort_size.png", width = 3200, height = 2000, res = 400)
p_k11_anc_scatter
dev.off()