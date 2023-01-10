source("Functions.R")

# This code retrieves mol/L from peak areas of targeted compounds.
# This run is for the MESOSCOPE Eddy Transect coA run.
# MSDial, QE, CoA

# Volume of seawater filtered for all TRANSECT samples was 5 L.

# User data ---------------------------------------------------------------
# Column.Type = "CoA"
Compound.ID = "Co-enzyme A"
QC.pattern = "QC"
# BMIS.pattern = "BMIS"
Dilution.Factor = 2
Injection.Volume = 400 # nanomoles
Volume.Filtered = 10 # liters

# Import standards and filter NAs ---------------------------------------------------------------
Ingalls.Standards <- read.csv("data_extras/Ingalls_Lab_Standards.csv", stringsAsFactors = FALSE) %>%
  filter(Compound.Type == Compound.ID) %>%
  rename(Metabolite.name = Compound.Name) %>%
  select(Metabolite.name, Compound.Type, QE.RF.ratio, Conc..uM, Emperical.Formula) %>%
  mutate(QE.RF.ratio = as.numeric(QE.RF.ratio)) %>%
  filter(!is.na(Conc..uM))

# Import sample key ---------------------------------------------------------------
Sample.Key <- read.csv("data_extras/Sample.Key.EddyTransect.CoA.csv", stringsAsFactors = F)

# Import QC'd files and remove parameter data ------------------------------
filename <- RemoveCsv(list.files(path = 'data_processed/', pattern = QC.pattern))
filepath <- file.path('data_processed', paste(filename, ".csv", sep = ""))

CoA.transect.raw <- assign(make.names(filename), read.csv(filepath, stringsAsFactors = FALSE, header = TRUE)) %>%
  slice(-1:-6) %>%
  select(-c(Description, Value)) %>%
  select(Replicate.Name, Metabolite.name, Area.with.QC, Area.Value, Run.Type)


# Apply appropriate filters and isolate standards ---------------------------------------------------------------
CoA.transect <- CoA.transect.raw %>%
  filter(Metabolite.name %in% Ingalls.Standards$Metabolite.name) %>%
  filter(str_detect(Replicate.Name, "Std")) %>%
  left_join(Ingalls.Standards, by = "Metabolite.name") %>%
  select(Replicate.Name, Metabolite.name, Compound.Type, everything()) %>%
  unique()


# Check standard run types ----------------------------------------------------
CoA.transect <- CheckStandards(CoA.transect)

# Get response factors for transect compounds ----------------------------------
CoA.RF.transect <- CoA.transect %>%
  mutate(RF = Area.with.QC/Conc..uM) %>%
  filter(!Compound.Type == "Internal Standard") %>%
  mutate(Replicate.Name = substr(Replicate.Name, 1, nchar(Replicate.Name)-2))


# Calculate RF max and min using only standards in water.
CoA.RF.dimensions <- CoA.RF.transect %>%
  filter(Type == "Standards_Water") %>%
  group_by(Metabolite.name) %>%
  mutate(RF.max = max(RF, na.rm = TRUE),
         RF.min = min(RF, na.rm = TRUE)) %>%
  rowwise() %>% mutate(RF.ave = mean(c(RF.max, RF.min), na.rm = T))

CoA.RF.dimensions <- CoA.RF.dimensions %>%
  mutate(RF.diff = RF.max/RF.min) %>%
  unique()

# Calculate response factor ratios ----------------------------------------
# Usually, calculate the response factor ratios using (Standards in Matrix - Water in Matrix) / Standards in water for each replicate.
# Because there is poor data for this run (no data for Water in Matrix), this time the data is drawn from Ingalls Standards.

# Quantify standards using response factor ratio ---------------------------------
CoA.standards.quantified <- CoA.RF.transect %>%
  left_join(CoA.RF.dimensions %>% select(Metabolite.name, RF.max:RF.diff)) %>%
  rename(RF.ratio = QE.RF.ratio) %>%
  replace(is.na(.), 0) %>%
  mutate(umol.in.vial.ave = Area.with.QC/RF.ave/RF.ratio,
         umol.in.vial.max = Area.with.QC/RF.min/RF.ratio,
         umol.in.vial.min = Area.with.QC/RF.max/RF.ratio) %>%
  unique()
CoA.standards.quantified[CoA.standards.quantified == 0] <- NA
  


# Quantify samples for the BMIS'd dataset ---------------------------------
# Sample.data.BMIS$Run.Cmpd <- strsplit(Sample.data.BMIS$Run.Cmpd," ")[[1]][1]
#



IS.data.transect <- CoA.standards.quantified %>%
  mutate(umol.in.vial_IS = NA) %>%
  select(Replicate.Name, Metabolite.name, Area.with.QC, Conc..uM, RF:umol.in.vial_IS) 


# Calculate umol/vial for samples (using the matched internal standard) -----------------
IS.smp.data.transect <- CoA.transect.raw %>%
  left_join(IS.data.transect %>% select(Metabolite.name, Conc..uM)) %>%
  unique() %>%
  filter(!str_detect(Replicate.Name, "Std")) %>%
  mutate(Std.Type = ifelse(str_detect(Metabolite.name, ","), "Internal_std", "Standard")) %>%
  group_by(Replicate.Name) %>%
  mutate(umol.in.vial_IS = (Area.with.QC[Std.Type == "Standard"] / Area.with.QC[Std.Type == "Internal_std"]) * (Conc..uM[Std.Type == "Standard"])) %>%
  select(Replicate.Name, Metabolite.name, Area.with.QC, umol.in.vial_IS) 



# Add transect matched IS_smp info back into main frame ------------------
all.info <- IS.smp.data.transect %>%
  bind_rows(IS.data.transect) %>%
  arrange(desc(Replicate.Name)) %>%
  filter(!Metabolite.name == "Acetyl CoA, 13C2")
  
  
# Add in dilution factor and filtered volume ------------------------------
all.info.quant <- all.info %>%
  mutate(nmol.in.Enviro = (umol.in.vial_IS*10^-6*Injection.Volume/Volume.Filtered*1000*Dilution.Factor)) %>% 
  left_join(Ingalls.Standards %>% select(Metabolite.name, Emperical.Formula)) %>%
  unique()

# Get molecules of carbon and nitrogen ------------------------------------
all.info.molecules <- all.info.quant  %>%
  mutate(C = ifelse(is.na(str_extract(Emperical.Formula, "^C\\d\\d")),
                    str_extract(Emperical.Formula, "^C\\d"),
                    str_extract(Emperical.Formula, "^C\\d\\d"))) %>%
  mutate(C = as.numeric(str_replace_all(C, "C", ""))) %>%
  mutate(N = ifelse(str_detect(Emperical.Formula, "N\\D"),
                    1,
                    str_extract(Emperical.Formula, "N\\d"))) %>%
  mutate(N = as.numeric(str_replace_all(N, "N", ""))) %>%
  mutate(nmol.C.ave = nmol.in.Enviro*C,
         nmol.N.ave = nmol.in.Enviro*N)
  #separate(Sample.Name, into = c("Date","type", "SampID", "replicate"), remove = FALSE) %>%
  #select(Metabolite.name, SampID, Replicate.Name, everything())


# Summarize for each metabolite ------------------------------------
all.info.summed <- all.info.molecules %>%
  group_by(Metabolite.name) %>%
  summarise(nmol.Enviro.med = median(nmol.in.Enviro, na.rm  = T),
            nmol.Enviro.min = min(nmol.in.Enviro, na.rm  = T),
            nmol.Enviro.max = max(nmol.in.Enviro, na.rm  = T),
            nmol.C.med = median(nmol.C.ave, na.rm  = T),
            nmol.C.min = min(nmol.C.ave, na.rm  = T),
            nmol.C.max = max(nmol.C.ave, na.rm  = T)) %>%
  arrange(desc(nmol.Enviro.med))


# Calculate mole fractions of each compound ------------------------------------
quantitative.perID <- all.info.molecules %>%
  separate(Replicate.Name, into = c("Date", "Run", "SampID", "Extras")) %>%
  select(SampID, nmol.C.ave, nmol.N.ave) %>%
  group_by(SampID) %>%
  summarise(totalCmeasured_nM_perID = sum(as.numeric(nmol.C.ave), na.rm = TRUE),
            totalNmeasured_nM_perID = sum(as.numeric(nmol.N.ave), na.rm = TRUE))

quantitative.final <- all.info.molecules %>%
  unique() %>%
  separate(Replicate.Name, into = c("Date", "Run", "SampID", "Extras"), remove = FALSE) %>%
  left_join(quantitative.perID) %>%
  mutate(ratioCN = totalCmeasured_nM_perID / totalNmeasured_nM_perID) %>%
  mutate(molFractionC = nmol.C.ave/totalCmeasured_nM_perID,
         molFractionN = nmol.N.ave/totalNmeasured_nM_perID) %>%
  select(Replicate.Name, Metabolite.name:umol.in.vial_IS, umol.in.vial.ave:nmol.in.Enviro, nmol.C.ave:molFractionN) %>%
  unique()

currentDate <- Sys.Date()
csvFileName.final <- paste("data_processed/Quantified_EddyTransect_CoA_Measurements_", currentDate, ".csv", sep = "")
csvFileName.perID <- paste("data_processed/Quantified_EddyTransect_CoA_perSampID_", currentDate, ".csv", sep = "")

write.csv(all.info.molecules, csvFileName.final, row.names = FALSE)
write.csv(quantitative.perID, csvFileName.perID, row.names = FALSE)

rm(list = ls())