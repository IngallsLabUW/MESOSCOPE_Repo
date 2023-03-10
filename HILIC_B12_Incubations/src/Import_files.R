source("src/Functions.R")

# Import all MSDial files --------------------------------------------------

filenames <- RemoveCsv(list.files(path = 'data_raw', pattern = '*.csv'))

for (i in filenames) {
  filepath <- file.path('data_raw', paste(i, ".csv", sep = ""))
  assign(make.names(i), read.csv(filepath, stringsAsFactors = FALSE))
}

matching.variable <- "hilic"


columns.to.drop <- c('Average.Rt.min.', 'Formula', 'Ontology', 'INCHIKEY', 
                     'SMILES', 'Isotope.tracking.parent.ID', 'Isotope.tracking.weight.number',
                     'MS1.isotopic.spectrum', 'MS.MS.spectrum', 'Average.Mz', 'Post.curation.result', 
                     'Fill..', 'Annotation.tag..VS1.0.', 'RT.matched',
                     'm.z.matched', 'MS.MS.matched', 'Manually.modified', 'Total.score', 
                     'RT.similarity', 'Dot.product', 'Reverse.dot.product', 'Fragment.presence..')

# Set header, filter unknowns ---------------------------------------
runs <- grep(matching.variable, names(.GlobalEnv), value = TRUE, ignore.case = TRUE)
runlist <- do.call("list", mget(runs))

headers.set <- lapply(names(runlist), function(x) SetHeader(runlist[[x]]))
names(headers.set) <- runs

for (df in seq_along(headers.set)) { 
  headers.set[[df]] <- headers.set[[df]] %>% filter(!Metabolite.Name == "Unknown")
  headers.set[[df]] <- headers.set[[df]] %>% select(-one_of(columns.to.drop))
}

# Change variable classes -------------------------------------------------
classes.changed <- lapply(names(headers.set), function(x) ChangeClasses(df = headers.set[[x]], start.column = 10))
names(classes.changed) <- runs

list2env(classes.changed, globalenv())


# Rearrange datasets ------------------------------------------------------

# Positive
Area.pos <- RearrangeDatasets(Area_HILIC.POS_B12.Inc, parameter = "Area.Value")
Mz.pos   <- RearrangeDatasets(Mz_HILIC.POS_B12.Inc, parameter = "Mz.Value")
RT.pos   <- RearrangeDatasets(RT_HILIC.POS_B12.Inc, parameter = "RT.Value")
SN.pos   <- RearrangeDatasets(SN_HILIC.POS_B12.Inc, parameter = "SN.Value")


# Negative
Area.neg <- RearrangeDatasets(Area_HILIC.NEG_B12.Inc, parameter = "Area.Value")
Mz.neg   <- RearrangeDatasets(Mz_HILIC.NEG_B12.Inc, parameter = "Mz.Value")
RT.neg   <- RearrangeDatasets(RT_HILIC.NEG_B12.Inc, parameter = "RT.Value")
SN.neg   <- RearrangeDatasets(SN_HILIC.NEG_B12.Inc, parameter = "SN.Value")


# Combine to one dataset --------------------------------------------------
combined.pos <- Area.pos %>%
  left_join(Mz.pos) %>%
  left_join(SN.pos) %>%
  left_join(RT.pos) %>%
  mutate(Column = "HILICPos") %>%
  select(Replicate.Name, Column, Area.Value, Mz.Value, RT.Value, SN.Value, everything())

combined.neg <- Area.neg %>%
  left_join(Mz.neg) %>%
  left_join(SN.neg) %>%
  left_join(RT.neg) %>%
  mutate(Column = "HILICNeg") %>%
  select(Replicate.Name, Column, Area.Value, Mz.Value, RT.Value, SN.Value, everything())

combined <- rbind(combined.pos, combined.neg) %>%
  mutate(Metabolite.Name = gsub(pattern = "Ingalls_", replacement = "", x = .$Metabolite.Name)) 

currentDate <- Sys.Date()
csvFileName <- paste("data_processed/MSDial_combined_", currentDate, ".csv", sep = "")

#############################
options(scipen = 999)
test <- combined %>%
  select(Metabolite.Name, Replicate.Name, Area.Value) %>%
  filter(!str_detect(Replicate.Name, "Sept29QC|TruePooWeek1|TruePooWeek2|TruePooWeek3|TruePooWeek4|DSW700m|Process|Std")) %>%
  mutate(Replicate.Name = recode(Replicate.Name, 
                                 "171002_Smp_IT0_1" ="171002_Smp_IL1IT0_1", 
                                 "171002_Smp_IT0_2" = "171002_Smp_IL1IT0_2",
                                 "171002_Smp_IT0_3" = "171002_Smp_IL1IT0_3",
                                 "171009_Smp_IT05um_1" = "171009_Smp_IL1IT05um_1",
                                 "171009_Smp_IT05um_2" = "171009_Smp_IL1IT05um_2",
                                 "171009_Smp_IT05um_3" = "171009_Smp_IL1IT05um_3",
                                 "171016_Smp_IT0_1" = "171016_Smp_IL2IT0_1",
                                 "171016_Smp_IT0_2" = "171016_Smp_IL2IT0_2",
                                 "171016_Smp_IT0_3" = "171016_Smp_IL2IT0_3",
                                 "171023_Smp_IT05um_1" = "171023_Smp_IL2IT05um_1",
                                 "171023_Smp_IT05um_2" = "171023_Smp_IL2IT05um_2",
                                 "171023_Smp_IT05um_3" = "171023_Smp_IL2IT05um_3")) %>%
  separate(Replicate.Name, into = c("one", "two", "SampID", "four"), fill = "right", remove = FALSE) %>%
  select(Metabolite.Name, SampID, Area.Value) %>%
  drop_na() %>%
  group_by(Metabolite.Name, SampID) %>%
  mutate(Averages = mean(Area.Value, na.rm = TRUE)) %>%
  select(Metabolite.Name, SampID, Averages) %>%
  unique() %>%
  mutate(Five_um = ifelse(str_detect(SampID, "5um"), TRUE, FALSE)) 

ggplot(data = test, aes(x = Metabolite.Name, y = Averages, fill = Five_um)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(plot.subtitle=element_text(size=10, color="black"),
        axis.text.x = element_text(angle = 90))

write.csv(combined, csvFileName, row.names = FALSE)

