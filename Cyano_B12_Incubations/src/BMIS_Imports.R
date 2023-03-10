# MSDial BMIS Script

# Imports -----------------------------------------------------------------
# Sample Key
filename <- RemoveCsv(list.files(path = "data_extras/", pattern = file.pattern))
filepath <- file.path("data_extras", paste(filename, ".csv", sep = ""))

SampKey.all <- assign(make.names(filename), read.csv(filepath, stringsAsFactors = FALSE, header = TRUE)) %>%
  rename(Replicate.Name = Sample.Name) %>%
  mutate(Replicate.Name = Replicate.Name %>%
           str_replace("-",".")) 

# Internal Standards
Internal.Standards <- read.csv("https://raw.githubusercontent.com/IngallsLabUW/Ingalls_Standards/master/Ingalls_Lab_Standards_NEW.csv",
                               stringsAsFactors = FALSE, header = TRUE) %>%
  filter(Column == Column.Type) %>%
  filter(Compound.Type == "Internal Standard")
Internal.Standards$Compound.Name <- TrimWhitespace(Internal.Standards$Compound.Name)

# Internal.Standards <- assign(make.names(filename), read.csv(filepath, stringsAsFactors = FALSE, header = TRUE)) %>%
#   filter(Column == Column.Type) %>%
#   filter(Compound.Type == "Internal Standard")
# Internal.Standards$Compound.Name <- TrimWhitespace(Internal.Standards$Compound.Name)

# QC'd output
filename <- RemoveCsv(list.files(path = "data_processed/", pattern = software.pattern))
filepath <- file.path("data_processed", paste(filename, ".csv", sep = ""))

QCd.data <- assign(make.names(filename), read.csv(filepath, stringsAsFactors = FALSE, header = TRUE)) %>%
  slice(-1:-6) %>%
  select(-c(Description, Value)) %>%
  filter(!str_detect(Replicate.Name, "Blk|Std")) %>%
  mutate(Replicate.Name = Replicate.Name %>%
           str_replace("-",".")) 

