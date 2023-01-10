source("Functions.R")

myfilename <- "EddyCenter_CoAs"

# Import Skyline file --------------------------------------------------
filenames <- RemoveCsv(list.files(path = 'data_raw', pattern = '*.csv'))

for (datafile in filenames) {
  filepath <- file.path('data_raw', paste(datafile, ".csv", sep = ""))
  assign(datafile, read.csv(filepath, stringsAsFactors = FALSE))
}

CoAs <- get(datafile)

# Adjust file for analysis--------------------------------------------------

CoAs <- CoAs %>%
  select(-Protein.Name, -Start.Time, -End.Time)

# Export Skyline file--------------------------------------------------
currentDate <- Sys.Date()
csvFileName <- paste("data_processed/Skyline_", myfilename, "_", currentDate, ".csv", sep = "")

write.csv(CoAs, csvFileName, row.names = FALSE)

rm(list = ls())







