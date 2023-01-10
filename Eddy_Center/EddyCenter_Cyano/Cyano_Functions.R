library(ggplot2)
library(plyr)
library(reshape2)
library(rlist)
library(stringr)
library(tidyverse)
library(tidyr)
options(scipen=999)
options(digits=3)


CheckFragments <- function(skyline.file, runtype) { 
  # Modifies transformed skyline output to prepare for standard ion ratio detection by running several tests.
  # 1. Isolate standards.
  # 2. Isolate unique Product.Mz fragments per compound.
  # 3. Check if each compound has two unique fragments.
  # 4. If it does, identify which fragment is quantitative and which is secondary by comparing them to the master compound list.
  # 5. Find 5% of the quantitative fragment.
  # 6. Determine if the secondary trace is > the 5% value from step 5.
  #
  # Args:
  #   skyline.file: Output from skyline that has had its variables modified to numeric values.
  #
  # Returns:
  #   fragments.checked: Modified data frame with added columns reflecting the above tests.
  #
  fragment.check <- skyline.output %>%
    filter(str_detect(Replicate.Name, runtype)) %>%
    select(Replicate.Name, Precursor.Ion.Name, Area, Precursor.Mz, Product.Mz)
  
  fragment.unique <-unique(fragment.check %>% select(Precursor.Ion.Name, Precursor.Mz, Product.Mz))
  
  fragment.multi.unique <- fragment.unique %>%
    count(Precursor.Ion.Name) %>%
    mutate(Two.Fragments = ifelse((n==1), FALSE, TRUE)) %>%
    select(-n)
  
  fragments.checked <- fragment.check %>%
    left_join(fragment.multi.unique, by = "Precursor.Ion.Name") %>%
    #left_join(standard.types, by = "Replicate.Name") %>%
    merge(y = master.file,
          by.x = c("Precursor.Ion.Name", "Product.Mz"),
          by.y = c("Compound.Name", "Daughter"),
          all.x = TRUE) %>%
    select(Replicate.Name, Precursor.Ion.Name, Area, Precursor.Mz, Product.Mz, Two.Fragments, Quan.Trace, Second.Trace) %>%
    mutate(Second.Trace = ifelse(Second.Trace == "", FALSE, TRUE)) %>%
    mutate(Quan.Trace = ifelse(Quan.Trace == "no", FALSE, TRUE)) %>%
    mutate(QT.Five.Percent = ifelse((Two.Fragments == TRUE & Quan.Trace == TRUE), 0.05 * Product.Mz, NA)) %>%
    mutate(Significant.Size = QT.Five.Percent < Product.Mz) %>%
    group_by(Replicate.Name, Precursor.Ion.Name, Product.Mz) %>%
    summarise_all(first) %>%
    arrange(Precursor.Ion.Name) %>%
    select(Replicate.Name, Precursor.Ion.Name, Precursor.Mz, Product.Mz, Area, Two.Fragments:Significant.Size)
  
  return(fragments.checked)
}

IdentifyRunTypes <- function(skyline.file) {
  # Identify run typfes and return each unique value present in the Skyline output.
  #
  # Args
  #   skyline.file: Raw output file from Skyline.
  #
  # Returns
  #   run.type: list of labels identifying the run types, isolated from Replicate.Name.
  #   Options conssist of samples (smp), pooled (poo), standards (std), and blanks (blk).
  #
  run.type <- tolower(str_extract(skyline.file$Replicate.Name, "(?<=_)[^_]+(?=_)"))
  print(paste("Your runtypes are:", toString(unique(run.type))))
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





CheckSmpFragments <- function(areas.transformed) {
  # Modifies transformed skyline output to prepare for standard ion ratio detection by running several tests.
  # 1. Isolate samples and pooled runs.
  # 2. Isolate unique Product.Mz fragments per compound.
  # 3. Check if each compound has two unique fragments.
  # 4. If it does, identify which fragment is quantitative and which is secondary by comparing them to the master compound list.
  # 5. Find 5% of the quantitative fragment.
  # 6. Determine if the secondary trace is > the 5% value from step 5.
  # 7. Find ion ratio by dividing the area of the quantitative trace by the area of the secondary trace.
  #
  # Args:
  #   areas.transformed: Output from skyline that has had its variables modified to numeric values.
  #
  # Returns:
  #   all.samples.IR: Modified data frame with added columns reflecting the above tests.
  unique.smp.frags <- unique(all.samples %>% select(Precursor.Ion.Name, Precursor.Mz, Product.Mz))
  
  unique.smp.frags2 <- unique.smp.frags %>%
    count(Precursor.Ion.Name) %>%
    mutate(Two.Fragments = ifelse((n==1), FALSE, TRUE)) %>%
    select(-n)
  
  all.samples.IR <- all.samples %>%
    left_join(unique.smp.frags2, by = "Precursor.Ion.Name" ) %>%
    merge(y = master.file,
          by.x = c("Precursor.Ion.Name", "Product.Mz"),
          by.y = c("Compound.Name", "Daughter"),
          all.x = TRUE) %>%
    select(Replicate.Name, Precursor.Ion.Name, Protein.Name, Precursor.Mz, Product.Mz, Area, Two.Fragments, Quan.Trace, Second.Trace) %>%
    mutate(Second.Trace = ifelse(Second.Trace == "", FALSE, TRUE)) %>%
    mutate(Quan.Trace = ifelse(Quan.Trace == "no", FALSE, TRUE)) %>%
    mutate(QT.Five.Percent = ifelse((Two.Fragments == TRUE & Quan.Trace == TRUE), 0.05 * Product.Mz, NA)) %>%
    mutate(Significant.Size = QT.Five.Percent < Product.Mz) %>%
    group_by(Precursor.Ion.Name) %>%
    mutate(IR.Ratio = ifelse(TRUE %in% Significant.Size, (Area[Quan.Trace == TRUE] / Area[Second.Trace == TRUE]), NA))
  return(all.samples.IR)
}


TrimWhitespace <- function (x) gsub("^\\s+|\\s+$", "", x)

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