# Code 4: Full climate data ----

rm(list=ls())
## Settings ----
source("00_Code/0.1 Settings.R")
source("00_Code/0.2 Packages.R")
source("00_Code/0.3 Functions.R")

# Data path 
data_inp <- "01_Data/Output/"
data_out <- "01_Data/Output/"

## Load data ----
mun <- rio::import(paste0(data_inp, "Data_births_ptb_mun.xlsx"))
cont <- rio::import(paste0(data_inp, "Contamination_2010_2020_series.RData"))
temp <- rio::import(paste0(data_inp, "Temp_NDVI_2010_2020_series.RData")) |> 
  filter(com %in% mun$com)

glimpse(cont); summary(cont)
glimpse(temp); summary(temp)

# Normalize municipality names
length(unique(cont$comuna))
length(unique(temp$name_com))

edit_com <- c(
  "Conchalí" = "Conchali",
  "Estación Central" = "Estacion Central",
  "Maipú" = "Maipu",
  "Ñuñoa" = "Nunoa",
  "Peñalolén" = "Penalolen",
  "San Joaquín" = "San Joaquin",
  "San Ramón" = "San Ramon"
)

cont <- cont |> 
  mutate(
    name_com = recode(comuna, !!!edit_com, .default = comuna)
  ) |> 
  select(-comuna)

setdiff(unique(temp$name_com), unique(cont$name_com))
setdiff(unique(cont$name_com), unique(temp$name_com))

glimpse(temp) # 132495
summary(temp)

glimpse(cont) # 144639
summary(cont)

n_dates <- length(seq(as.Date("2010-01-01"), as.Date("2020-12-31"), by = "days"))
n_dates*33 # 4018 -> 132594

data_clime <- temp |> 
  left_join(cont, by = c("name_com", "date"))

glimpse(data_clime) # 132495
summary(data_clime)

test <- data_clime |> 
  group_by(com) |> 
  summarise(n = n())

# Differecen by 29/02 -> 2012, 2016, 2020
setdiff(
  seq(as.Date("2010-01-01"), as.Date("2020-12-31"), by = "days"), 
  unique(data_clime$date)
)

# Save data 
save(data_clime, file = paste0(data_out, "Contamination_Climate_Data_2010_2020.RData"))

