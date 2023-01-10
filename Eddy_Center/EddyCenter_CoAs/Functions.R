library(plyr)
library(reshape2)
library(stringr)
library(tidyverse)
options(scipen=999)


RemoveCsv <- function(full.filepaths) {
  # Gathers all files in given directory and drops the csv extension.
  #
  # Args
  #   full.filepaths: list of files in a directory matching given patterns.
  #
  # Returns
  #   no.path: list of files, with filepath and csv extension removed.
  #
  no.path <- substr(full.filepaths, 1, nchar(full.filepaths)-4)
  no.ID <-   gsub("\\_.*","", no.path)
  
  return(no.path)
}

StandardizeMetabolites <- function(df) {
  df.standardized <- df %>%
    mutate(Metabolite.Name = ifelse(str_detect(Metabolite.Name, "Ingalls_"), sapply(strsplit(Metabolite.Name, "_"), `[`, 2), Metabolite.Name)) 
  
  df.standardized$Replicate.Name <- gsub("^.{0,1}", "", df.standardized$Replicate.Name)
  
  return(df.standardized)
}
