# switch to working in R and define inputs
library(tidyverse)
library(tidyr)
library(patchwork)
library(grid)

# get all the records from the master spreadsheet
records <- read.csv("/SSD/reg259/all_birds_highcov_lowcov/Finch_Database_V2.csv", header = TRUE, fill = TRUE, na.strings = c("", "NA", "XX"), stringsAsFactors = FALSE)
records <- records %>% filter(Island == "Daphne")
records <- records[,c("meaningful.unique", "Species", "bill.depth", "bill.length", "bill.width", "First.year", "Last.year")]
records <- records %>% filter(Species %in% c("fortis", "scandens"))
records <- records %>% rename(ind=meaningful.unique, depth=bill.depth, length=bill.length, width=bill.width, start=First.year, end=Last.year)
records <- records %>% mutate(end = coalesce(end, start))  # if no end year, use start year
records <- records[!duplicated(records$ind),]  # remove duplicates (probably they are all re-captures)

# --- Define a function to read and prepare ancestry data ---
prepare_ancestry_data <- function(k_value, con, fort, scan) {
  # read the anccestry data
  data <- read.table(paste0("finches.", k_value, ".Q"))
  fam <- read.table("finches.fam")
  data$ind <- fam$V2
  
  # add species and year information
  data$Species <- records$Species[match(data$ind, records$ind)]
  data <- data %>% filter(Species %in% c("fortis", "scandens"))
  data$start <- records$start[match(data$ind, records$ind)]
  data$end <- records$end[match(data$ind, records$ind)]
  data <- data[!is.na(data$start),]
  data$end[is.na(data$end)] <- data$start[is.na(data$end)]  # end year is start year if not specified

  # remove individuals from the A line
  a_line_inds <- c("01Dap19800", "01Dap19210", "01Dap19565", "01Dap19644")
  data <- data[!(data$ind %in% a_line_inds),]
  
  # calculate conirostris ancestry as a function of "main" ancestry
  data <- data %>% rename(anc_con = all_of(con), anc_fort = all_of(fort), anc_scan = all_of(scan))
  data <- data %>% mutate(anc_main = ifelse(Species == "fortis", anc_fort, anc_scan))
  data$con_modified <- data$anc_con / (1 - data$anc_main)
  data <- data[, c("ind", "Species", "start", "end", "anc_con", "con_modified")]
  
  # expand the data and summarize
  df_expanded <- data %>% rowwise() %>% mutate(year = list(start:end)) %>% unnest(year)
  df_summary <- df_expanded %>% group_by(Species, year) %>% summarise(n_alive = n(), mean = mean(anc_con, na.rm = TRUE), max = max(anc_con, na.rm = TRUE),
    sd = sd(anc_con, na.rm = TRUE), .groups = "drop") %>% mutate(Species = recode(Species, fortis = "italic('fortis')", scandens = "italic('scandens')"))
  
  return(df_summary)
}

# --- Define a function to plot the data ---
plot_con_ancestry <- function(summary_df, title_text, show_y = TRUE, show_legend = TRUE) {
  species_colors <- c("italic('scandens')" = "seagreen3", "italic('fortis')" = "dodgerblue2")
  
  p <- ggplot(summary_df, aes(x = year, y = max, color = Species)) + geom_line(linewidth = 1) + geom_point(aes(size = n_alive), show.legend = FALSE) +
    #geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.3, alpha = 0.7) +
    scale_color_manual(values = species_colors, labels = scales::parse_format()) +
    scale_size_continuous(name = "Number Alive") + theme_minimal(base_size = 15) + labs(title = title_text, x = "Year", color = "Species") +
    theme(plot.title = element_text(hjust = 0.5, face = "plain"), legend.position = if (show_legend) "right" else "none", axis.text.y = if (!show_y) element_blank() else NULL,
      axis.ticks.y = if (!show_y) element_blank() else NULL, axis.title.y = if (!show_y) element_blank() else element_text(margin = margin(r = 10)))
  
  if (show_y) {p <- p + labs(y = "Maximum Conirostris Ancestry")}
  return(p)
}

# --- Prepare summaries for k=10 and k=11 ---
#summary_k10 <- prepare_ancestry_data(k_value = 10, con = "V9", fort = "V10", scan = "V5")
summary_k11 <- prepare_ancestry_data(k_value = 11, con = "V8", fort = "V9", scan = "V3")

# --- Create plots with control over legend and y-axis ---
#plot_k10 <- plot_con_ancestry(summary_k10, title_text = "K=10", show_y = TRUE, show_legend = FALSE)
plot_k11 <- plot_con_ancestry(summary_k11, title_text = "K=11", show_y = TRUE, show_legend = TRUE)

# --- Save combined plot ---
png("max_con_anc_time_k11.png", width = 3000, height = 2200, res = 400)
#plot_k10 + plot_k11
plot_k11
dev.off()

# for identifying top individuals
#inds <- c()
#for (year in 1976:2022) {
#  top <- df_expanded[(df_expanded$Species == "fortis") & (df_expanded$year == year),] %>% arrange(desc(anc_con)) %>% head(1)
#  inds <- c(inds, top$ind)
#  }
#print(inds)