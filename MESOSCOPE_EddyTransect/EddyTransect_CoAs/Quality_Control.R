# Quality control script
## CoA, MSDial, QE

source("Functions.R")

area.min   <- 1000
RT.flex    <- 0.4
blk.thresh <- 0.3
SN.min     <- 4

pattern = "combined"


# Import QC'd files and clean parameter data ----------------------------
filename <- RemoveCsv(list.files(path = 'data_processed/', pattern = pattern))
filepath <- file.path('data_processed', paste(filename, ".csv", sep = ""))

combined <- assign(make.names(filename), read.csv(filepath, stringsAsFactors = FALSE, header = TRUE)) %>%
  select(Replicate.Name:Alignment.ID, Metabolite.name) %>%
  mutate(Run.Type = (tolower(str_extract(Replicate.Name, "(?<=_)[^_]+(?=_)")))) 

msdial.runtypes <- IdentifyRunTypes(combined)

## Quick graph
coA.graphs <- combined %>%
  separate(Replicate.Name, c("runDate", "type", "SampID", "replicate"), "_") %>%
  select(SampID, Metabolite.name, Area.Value) %>%
  group_by(Metabolite.name, SampID) %>%
  mutate(Mean.Area = mean(Area.Value, na.rm = TRUE)) %>%
  filter(!str_detect(SampID, "StdMix")) %>%
  arrange(Metabolite.name)

ggplot(coA.graphs, aes(x = Metabolite.name, y = Mean.Area)) +
  geom_bar(stat = "identity", position = "dodge") 

# Create reference tables for Retention Time (RT) and Blanks (blk) --------------------------------

RT.table <- combined %>%
  filter(Run.Type == "std") %>%
  arrange(Metabolite.name) %>%
  group_by(Metabolite.name) %>%
  mutate(RT.min = min(RT.Value, na.rm = TRUE)) %>%
  mutate(RT.max = max(RT.Value, na.rm = TRUE)) %>%
  select(Metabolite.name:RT.max) %>%
  unique()

blank.table <- combined %>%
  filter(Run.Type == "blk") %>%
  mutate(Blk.Area = Area.Value) %>%
  arrange(Metabolite.name) %>%
  group_by(Metabolite.name) %>%
  mutate(Blk.min = min(Area.Value)) %>%
  mutate(Blk.max = max(Area.Value)) %>%
  select(Metabolite.name:Blk.max) %>%
  select(-Blk.Area) %>%
  unique()


# Create datasets for different flag types --------------------------------
SN.Area.Flags <- combined %>%
  arrange(Metabolite.name) %>%
  mutate(SN.Flag       = ifelse(((SN.Value) < SN.min), "SN.Flag", NA)) %>%
  mutate(Area.Min.Flag = ifelse((Area.Value < area.min), "Area.Min.Flag", NA))

# Joining datasets---------------------------------------
add.RT.Flag <- SN.Area.Flags %>%
  group_by(Metabolite.name) %>%
  left_join(RT.table, by = c("Metabolite.name", "Run.Type")) %>%
  mutate(RT.Flag = ifelse((RT.Value >= (RT.max + RT.flex) | RT.Value <= (RT.min - RT.flex)), "RT.Flag", NA)) %>%
  select(-c("RT.max", "RT.min"))

add.blk.Flag <- add.RT.Flag %>%
  left_join(blank.table, by = c("Metabolite.name", "Run.Type")) %>%
  mutate(Blank.Flag = ifelse((Area.Value / Blk.max) < blk.thresh, "Blank.Flag", NA)) %>%
  select(-c("Blk.min", "Blk.max"))


# Combine all the flags ---------------------------------------------------
final.table <- add.blk.Flag %>%
  mutate(all.Flags      = paste(SN.Flag, Area.Min.Flag, RT.Flag, Blank.Flag, sep = ", ")) %>%
  mutate(all.Flags      = as.character(all.Flags %>% str_remove_all("NA, ") %>% str_remove_all("NA"))) %>%
  mutate(all.Flags      = ifelse(all.Flags == "", NA, all.Flags)) %>%
  mutate(Area.with.QC   = ifelse(is.na(Area.Min.Flag), Area.Value, NA)) %>%
  select(Replicate.Name:Area.Value, Area.with.QC, everything()) %>%
  ungroup(Metabolite.name) %>%
  mutate(Metabolite.name = as.character(Metabolite.name)) 

#### Keep only Acetyl CoA and its matched internal standard ##
final.table <- final.table %>%
  filter(Metabolite.name %in% c("Acetyl CoA, 13C2", "Acetyl CoA"))


# Print to file with comments and a new name ------------------------------
Description <- c("Hello! Welcome to the world of MSDIAL QE Quality Control! ",
                 "Minimum area for a real peak: ",
                 "RT flexibility: ",
                 "Blank can be this fraction of a sample: ",
                 "S/N ratio: " ,
                 "Processed on: ")
Value <- c(NA, area.min, RT.flex, blk.thresh, SN.min, Sys.time())

df <- data.frame(Description, Value)
final.table <- bind_rows(df, final.table)

currentDate <- Sys.Date()
csvFileName <- paste("data_processed/QC_HILIC_Output_", currentDate, ".csv", sep = "")

write.csv(final.table, csvFileName, row.names = FALSE)

rm(list = ls())