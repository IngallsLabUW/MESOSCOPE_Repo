source("Functions.R")

All.Info.Quantitative <- read.csv("data_processed/Quantified_Measurements_2020-06-29.csv", 
                                  stringsAsFactors = FALSE)

# Isolate Gln/Glu ratios to check Wei's work ------------------------------
Gln_Glu <- All.Info.Quantitative %>%
  select(Metabolite.name, Replicate.Name, umol.in.vial.ave) %>%
  filter(str_detect(Metabolite.name, "Glutamic|Glutamine")) %>%
  filter(!str_detect(Replicate.Name, "Poo")) %>%
  filter(str_detect(Replicate.Name, "MS6|MS12")) %>%
  filter(str_detect(Replicate.Name, "15m")) %>%
  mutate(Eddy = ifelse(str_detect(Replicate.Name, "MS12"), "Cyclonic", "Anticyclonic")) %>%
  separate(Replicate.Name, into = c("Sample.Name", "type", "SampID", "replicate")) %>%
  group_by(Metabolite.name, SampID) %>%
  mutate(Total.umol.Ave = mean(umol.in.vial.ave, na.rm = TRUE)) %>%
  select(Metabolite.name, SampID, Eddy, Total.umol.Ave) %>%
  unique() %>%
  group_by(Eddy) %>%
  mutate(Ratios = (Total.umol.Ave[Metabolite.name == "Glutamine"]) 
                     / (Total.umol.Ave[Metabolite.name == "Glutamic acid"]))

ggplot(data = Gln_Glu, aes(Eddy, Ratios, fill = Eddy)) +
  geom_bar(stat = "identity", position = "dodge") +
  ggtitle("Eddy Transect Cyclonic + Anticyclonic")