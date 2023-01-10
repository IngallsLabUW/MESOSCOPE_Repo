library(ggplot2)
library(rlist)
library(stringr)
library(tidyverse)
library(tidyr)
options(scipen=999)


ChangeClasses <- function(df) {
  # Change specified columns from factors to numeric.
  #
  # Args
  #   df: MSDial dataframe containing sample columns.
  #
  # Returns
  #   df: MSDial dataframe, with specified columns changed to a numeric class. 
  for (i in c(10:ncol(df))) {
    df[, i] <- as.numeric(as.character(df[, i]))
  }
  return(df)
}

IdentifyDuplicates <- function(df) {
  # Determine which compounds are detected in both positive and negative HILIC runs.
  # 
  # Args
  #   df: MSDial dataframe, containing all required parameters (MZ, SN, Area, etc),
  #       and modified to long form instead of wide.
  # 
  # Returns
  #   duplicates: Simple dataframe of listed compounds that have been identified as duplicates.
  #
  duplicates <- df %>%
    group_by(Metabolite.name, Replicate.Name) %>%
    mutate(number = 1) %>%
    mutate(ticker = cumsum(number)) %>%
    filter(ticker == 2) %>%
    ungroup() %>%
    select(Metabolite.name) %>%
    unique()
  return(duplicates)
}

IdentifyRunTypes <- function(msdial.file) {
  # Identify run typfes and return each unique value present in the Skyline output.
  #
  # Args
  #   msdial.file: Raw output file from Skyline.
  #
  # Returns
  #   run.type: list of labels identifying the run types, isolated from Replicate.Name.
  #   Options conssist of samples (smp), pooled (poo), standards (std), and blanks (blk).
  #
  run.type <- tolower(str_extract(msdial.file$Replicate.Name, "(?<=_)[^_]+(?=_)"))
  print(paste("Your runtypes are:", toString(unique(run.type))))
}

RearrangeDatasets <- function(df, parameter) {
  df <- df %>%
    tidyr::gather(
      key = "Replicate.Name",
      value = "parameter",
      starts_with("X")) %>%
    select(Replicate.Name, parameter, everything())
  
  names(df)[2] <- parameter
  
  return(df)
}

RemoveCsv <- function(full.filepaths) {
  # Remove a .csv file extension and obtain basename from a given list of filepaths.
  #
  # Args
  #   Character strings of filepaths in a directory.
  #
  # Returns
  #   Character strings of file basenames, without a csv extension.
  #
  no.path <- substr(full.filepaths, 1, nchar(full.filepaths)-4)
  no.ID <-   gsub("\\_.*","", no.path)
  
  return(no.path)
}

SetHeader <- function(df) {
  # Remove empty or unnecessary lines from machine output, and make column names headers.
  #
  # Args
  #   df: Raw output file from MSDial.
  #
  # Returns
  #   df: modified dataframe with correct headers and no empty lines.
  #
  df <- df[!(is.na(df[1]) | df[1]==""), ]
  colnames(df) <- make.names(as.character(unlist(df[1,])))
  df <- df[-1, ]
  
  return(df)
}

TrimWhitespace <- function (x) gsub("^\\s+|\\s+$", "", x)








StandardizeMetabolites <- function(df) {
  df.standardized <- df %>%
    mutate(Metabolite.name = ifelse(str_detect(Metabolite.name, "Ingalls_"), sapply(strsplit(Metabolite.name, "_"), `[`, 2), Metabolite.name)) 
  
  df.standardized$Replicate.Name <- gsub("^.{0,1}", "", df.standardized$Replicate.Name)
  
  return(df.standardized)
}

CheckStandards <- function (df) {
  # Mutates a new column identifying standard run types, then prints number of unique run types.
  #
  # Args
  #   df: Dataset of containing a Replicate.Name column, pre-filtered to include only standard runs.
  #
  # Returns
  #   df.checked: Dataset with a new column describing run types, and a printed message stating how many 
  #               unique types there are.
  #
  df.checked <- df %>%
    mutate(Type = paste(Env = ifelse(str_detect(Replicate.Name, "StdMix"), "Standards", "Water"),
                        Matrix = ifelse(str_detect(Replicate.Name, "Matrix"), "Matrix", "Water"), sep = "_"))
  
  print(paste("Number of standard run types:", length(unique(df.checked$Type))))
  print(unique(df.checked$Type))
  
  return(df.checked)
}