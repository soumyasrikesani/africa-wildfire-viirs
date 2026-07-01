library(raster)
library(tidyverse)
library(sf)

# Set download and temp file paths relative to current working directory
curdir <- "./downloads"
tempdir <- "./temp"

# Read in file with list of African country names
africa_data <- read.csv("Africa_country_table.csv")

# List downloaded data zipfiles
zipfiles <- list.files(path = curdir,
                       pattern = "*.zip",
                       full.names = TRUE)

# Loop through years of data
years <- 2012:2024
for(i in 1:length(zipfiles)) {
  curzip <- zipfiles[i]     # Get zipfile for the current year
  curyear <- years[i]       # Get the current year
  # Generate a list of csv names for all African countries in the current zipfile
  africa_files <- paste0("viirs-snpp/", curyear, "/viirs-snpp_", 
                         curyear, "_", africa_data[,1], ".csv")
  # Extract the csv files for African countries to temp directory
  unzip(curzip, 
        files = africa_files, 
        exdir = tempdir,
        junkpaths = TRUE)
  # Generate a list of the csv file names for African countries in the temp directory
  csvfile <- list.files(path = tempdir,
                        pattern = "*.csv",
                        full.names = TRUE)
  # Combine the individual country csv files into a single data frame
  for(j in 1:length(csvfile)) {
    print(csvfile[j])
    incsv <- read.csv(csvfile[j])
    if(j == 1) {
      outcsv <- incsv
    } else {
      outcsv <- bind_rows(outcsv, incsv)
    }
  }
  # Write the combined CSV files for Africa into a single csv
  outname <- paste0("africa_fires_", curyear, ".csv")
  write.csv(outcsv, file = outname, row.names=F)
  # Delete the temporary files
  unlink(csvfile)
}