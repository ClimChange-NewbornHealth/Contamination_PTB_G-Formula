# 9.0 DLM results -----
rm(list = ls())

## Settings ----
source("00_Code/0.1 Settings.R")
source("00_Code/0.2 Packages.R")
source("00_Code/0.3 Functions.R")

data_inp <- "01_Data/Output/"
data_out <- "02_Output/Descriptives/"
data_out_model <- "02_Output/Models/"

## Load data ----

births_weeks <- rio::import(paste0(data_inp, "births_2010_2020_exposure_weeks.RData"))
glimpse(births_weeks)

## Descriptive exposition by gestational week ----

g_expo <- births_weeks |> 
  group_by(week_gest_num) |> 
  summarise(
    pm25_mean = mean(pm25_krg, na.rm = TRUE),
    pm25_min = min(pm25_krg, na.rm = TRUE),
    pm25_max = max(pm25_krg, na.rm = TRUE),

    no2_mean = mean(no2_krg, na.rm = TRUE),
    no2_min = min(no2_krg, na.rm = TRUE),
    no2_max = max(no2_krg, na.rm = TRUE),

    o3_mean = mean(o3_krg, na.rm = TRUE),
    o3_min = min(o3_krg, na.rm = TRUE),
    o3_max = max(o3_krg, na.rm = TRUE)
  ) |> 
  ungroup() |> 
  mutate(across(where(is.numeric), ~ format(round(.x, 2), nsmall = 2, decimal.mark = ".")))

write.xlsx(g_expo, paste0(data_out, "Descriptives_exposure_stats_time.xlsx"))

## Weighted DLM ----

