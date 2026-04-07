# Code 3.1: Temp data and NDVI ----

rm(list=ls())
## Settings ----
source("00_Code/0.1 Settings.R")
source("00_Code/0.2 Packages.R")
source("00_Code/0.3 Functions.R")

# Data path 
data_inp <- "01_Data/Input/Clime_series/"
data_out <- "01_Data/Output/"

data_dis <- rio::import(paste0(data_inp, "District_data_geo_RM.RData"))  |> 
  mutate(geometry_id = row_number() - 1) |> 
  dplyr::select(geometry_id, codigo_comuna)
glimpse(data_dis)

## NDVI Data ----

ndvi <- "ndvi_daily_district.csv"

green <- rio::import(paste0(data_inp, ndvi)) |> 
  mutate(date = as.Date(date, format = "%Y_%m_%d")) |> 
  mutate(ndvi=ndvi/10000) |> # Adjust scale -1, 1
  dplyr::select(geometry_id, date, ndvi)

summary(green)

ndvi <- data_dis |> 
  left_join(green, by = "geometry_id")

glimpse(ndvi)
summary(ndvi)

# Complete ndvi with mean between two dates
# Impute NA with interpolation and kalman filter 
# Missing daily NDVI values were imputed using a state-space model with Kalman smoothing, 
# which estimates the latent vegetation signal conditional on the observed MODIS composites, accounting for temporal dependence and measurement noise.
# Reference: Harvey, Andrew C. Forecasting, structural time series models and the Kalman filter. Cambridge university press, 1990

ndvi <- ndvi |> 
  group_by(geometry_id, codigo_comuna) |>
  mutate(
    ndvi_kalman = na_kalman(
      ndvi,
      model = "StructTS",
      smooth = TRUE
    )
  ) |>
  ungroup()

glimpse(ndvi)
summary(ndvi)

ndvi <- ndvi |> 
  dplyr::select(-ndvi, -geometry_id) |> 
  rename(ndvi=ndvi_kalman)

glimpse(ndvi)

## Temperature Data ----
# CR2Met 
temp <- "hw_data_1980_2021.RData"
temp <- rio::import(paste0(data_inp, temp)) 
glimpse(temp)

temp <- temp |> 
  dplyr::select(com:tmin, TAD) |> 
  filter(year >= 2010 & year <= 2020) 

glimpse(temp)

## Join and save data ----

temp <- temp |> 
  left_join(ndvi, by = c("com" = "codigo_comuna", "date")) 

glimpse(temp)
summary(temp)

# Save data 
save(temp, file = paste0(data_out, "Temp_NDVI_2010_2020_series.RData"))

