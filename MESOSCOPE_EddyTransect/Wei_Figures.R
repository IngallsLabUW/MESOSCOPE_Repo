source("Functions.R")

## Remaking Wei's slide figures

All.Info.Quantitative <- read.csv("data_processed/Quantified_Measurements_2020-06-29.csv", 
                                  stringsAsFactors = FALSE)

# Glutamic acid|Glutamine|Ketoglutaric|Acetyl-Lysine" to be included

makeLineGraph <- function(df, chosen.compound) {
  my.compounds <- df %>%
    select(Metabolite.name, Replicate.Name, umol.in.vial.ave) %>%
    filter(str_detect(Metabolite.name, chosen.compound),
          !str_detect(Replicate.Name, "Poo"),
           str_detect(Replicate.Name, "15m|175m|DCM")) %>%
    separate(Replicate.Name, into = c("Sample.Name", "type", "SampID", "replicate")) %>%
    mutate(Depth = ifelse(str_detect(SampID, "15m"), "15m",
                          ifelse(str_detect(SampID, "175m"), "175m",
                                 ifelse(str_detect(SampID, "DCM"), "DCM", NA)))) %>%
    group_by(SampID) %>%
    mutate(Station = strsplit(SampID, "[C]")[[1]][1]) %>%
    group_by(Metabolite.name, SampID) %>%
    mutate(StdDevs = sd(umol.in.vial.ave, na.rm = TRUE)) %>%
    group_by(Metabolite.name, SampID) %>%
    mutate(Averages = mean(umol.in.vial.ave, na.rm = TRUE)) %>%
    select(-type, -replicate, -Sample.Name, -umol.in.vial.ave) %>%
    unique()

  my.compounds$Station <- factor(my.compounds$Station, 
                                 levels = c("MS4", "MS5", "MS6", "MS7", "MS8", "MS9",
                                            "MS10", "MS11", "MS12", "MS13", "MS14"))

  myplot <- ggplot(my.compounds, aes(x = Station, y = Averages, group = Depth, color = Depth)) + 
    geom_line(size=1) + 
    geom_point()+
    geom_errorbar(aes(ymin=Averages-StdDevs, ymax=Averages+StdDevs), width=.2,
                  position=position_dodge(0.05)) +
    ggtitle(paste(chosen.compound, "EddyTransect", sep = "_"))
  print(myplot)
}

test <- makeLineGraph(All.Info.Quantitative, "Acetyl-Lysine")