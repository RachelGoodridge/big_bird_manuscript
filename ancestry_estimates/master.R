# the master database for the ancestry plots: exactly the table of load_master(island = None) in tune_trim_ibis/scripts/common.py (one row per
# bird with its entries combined, the newest lineage classifications, SPECIES_CORRECTIONS, measurements averaged with NESTLING_MEASUREMENTS
# removed, first / last years alive), produced by that function itself so the R and Python analyses can never disagree.
# Use: source("../../master.R") from a script two folders down, then load_master()
REPO <- normalizePath(file.path(dirname(sys.frame(1)$ofile), ".."))

load_master <- function() {
  table <- tempfile(fileext = ".csv")
  on.exit(unlink(table))
  code <- paste0("import sys; sys.path.insert(0, '", REPO, "/tune_trim_ibis/scripts'); from common import load_master; ",
                 "load_master(island=None).to_csv(sys.argv[1], index=False)")
  status <- system2("python", c("-c", shQuote(code), table))
  if (status != 0) stop("load_master() in tune_trim_ibis/scripts/common.py failed")
  read.csv(table, na.strings = c("", "NA"), stringsAsFactors = FALSE, check.names = FALSE)
}

# the years without fieldwork on Daphne (UNSAMPLED_YEARS in the same file), left out wherever birds are counted in every year they were alive
unsampled_years <- function() {
  code <- paste0("import sys; sys.path.insert(0, '", REPO, "/tune_trim_ibis/scripts'); from common import UNSAMPLED_YEARS; ",
                 "print(*UNSAMPLED_YEARS)")
  as.integer(strsplit(system2("python", c("-c", shQuote(code)), stdout = TRUE), " ")[[1]])
}
