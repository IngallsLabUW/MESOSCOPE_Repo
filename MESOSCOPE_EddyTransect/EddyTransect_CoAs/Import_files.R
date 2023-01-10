source("Functions.R")

# User input --------------------------------------------------

matching.variable <- "CoA"

columns.to.drop <- c('Average.Rt.min.', 'Formula', 'Ontology', 'INCHIKEY', 
                     'SMILES', 'Isotope.tracking.parent.ID', 'Isotope.tracking.weight.number',
                     'MS1.isotopic.spectrum', 'MS.MS.spectrum', 'Average.Mz', 'Post.curation.result', 
                     'Fill..', 'Annotation.tag..VS1.0.', 'RT.matched',
                     'm.z.matched', 'MS.MS.matched', 'Manually.modified', 'Total.score', 
                     'RT.similarity', 'Dot.product', 'Reverse.dot.product', 'Fragment.presence..')

# Import all MSDial files --------------------------------------------------
filenames <- RemoveCsv(list.files(path = 'data_raw', pattern = '*.csv'))

for (i in filenames) {
  filepath <- file.path('data_raw', paste(i, ".csv", sep = ""))
  assign(make.names(i), read.csv(filepath, stringsAsFactors = FALSE))
}


# Set header, filter unknowns ---------------------------------------
runs <- grep(matching.variable, names(.GlobalEnv), value = TRUE, ignore.case = TRUE)
runlist <- do.call("list", mget(runs))

headers.set <- lapply(names(runlist), function(x) SetHeader(runlist[[x]]))
names(headers.set) <- runs


for (df in seq_along(headers.set)) { 
  headers.set[[df]] <- headers.set[[df]] %>% filter(!Metabolite.name == "Unknown")
  headers.set[[df]] <- headers.set[[df]] %>% select(-one_of(columns.to.drop))
}


# Change variable classes -------------------------------------------------
classes.changed <- lapply(names(headers.set), function(x) ChangeClasses(headers.set[[x]]))
names(classes.changed) <- runs

list2env(classes.changed, globalenv())


# Rearrange datasets ------------------------------------------------------
Area <- RearrangeDatasets(Area_CoA_EddyTransect, parameter = "Area.Value")
Mz   <- RearrangeDatasets(Mz_CoA_EddyTransect, parameter = "Mz.Value")
RT   <- RearrangeDatasets(RT_CoA_EddyTransect, parameter = "RT.Value")
SN   <- RearrangeDatasets(SN_CoA_EdddyTransect, parameter = "SN.Value")


# Combine to one dataset --------------------------------------------------
combined <- Area %>%
  left_join(Mz) %>%
  left_join(SN) %>%
  left_join(RT) %>%
  select(Replicate.Name, Area.Value, Mz.Value, RT.Value, SN.Value, everything()) %>%
  mutate(Metabolite.name = ifelse(str_detect(Metabolite.name, "Ingalls_"), sapply(strsplit(Metabolite.name, "_"), `[`, 2), Metabolite.name)) 

combined$Replicate.Name <- gsub("^.{0,1}", "", combined$Replicate.Name)

currentDate <- Sys.Date()
csvFileName <- paste("data_processed/MSDial_combined_", currentDate, ".csv", sep = "")

write.csv(combined, csvFileName, row.names = FALSE)

rm(list = ls())