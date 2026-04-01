# Code 1: Adjust contamination interpolation ----

## Settings ----
source("00_Code/0.1 Settings.R")
source("00_Code/0.2 Packages.R")
source("00_Code/0.3 Functions.R")

# Data path 
data_inp <- "01_Data/Input/Contant_series/"
data_out <- "01_Data/Output/"

## Contamination Data ---- 

# Load data and generate variables
cont <- rio::import(paste0(data_inp, "interpolated_series.RData")) |>
  select("comuna", "date", ends_with("pred"))

glimpse(cont)
summary(cont)
unique(cont$comuna)

# Edit minimum detectable values:
# PM: ATENUACION BETA-MET ONE 1020, daily detection 0.98µg.
# Ozono: FOTOMETRIA UV - THERMO 49i, < 1 ppb
# NO2: 0,01 ppm (tiempo promedio de 60 segundos)
cont <- cont |> 
  mutate(
    across(starts_with("pm25"), ~ ifelse(. <= 0.98, 0.98, .)),
    across(starts_with("o3"), ~ ifelse(. < 1, 1, .)),
    across(starts_with("no2"), ~ ifelse(. < 0.01, 0.01, .))
  )

summary(cont)

# Save data 
save(cont, file = paste0(data_inp, "Contamination_2009_2020_series.RData"))
