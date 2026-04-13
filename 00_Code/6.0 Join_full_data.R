# 6.0 Join data exposure and births -----
rm(list = ls())

## Settings ----
source("00_Code/0.1 Settings.R")
source("00_Code/0.2 Packages.R")
source("00_Code/0.3 Functions.R")

data_inp <- "01_Data/Output/"
data_out <- "01_Data/Output/"
data_sovi <- "01_Data/Input/SOVI/"

## Open data ----

births <- rio::import(paste0(data_out, "births_2010_2020.RData"))
glimpse(births)

sovi <- rio::import(paste0(data_sovi, "sovi_datasets", ".RData")) |> 
  select(-name_comuna)  |> 
  rename(vulnerability=vulnerablidad) |> 
    mutate(vulnerability = fct_recode(vulnerability,
      "Low" = "Baja",
      "Medium-low" = "Medio-baja",
      "Medium-high" = "Medio-alta"))

sovi <- sovi |> 
  rename(com = cod_com)

glimpse(sovi)

exposure <- rio::import(paste0(data_out, "births_exposure_period_metrics_full_d30_d4_tri.RData"))
glimpse(exposure)

week_exposure <- rio::import(paste0(data_out, "births_exposure_gest_window_means.RData"))
glimpse(week_exposure)

## Join data exposure ----

exposure <- exposure |> 
  select(id, com, pm25_krg_full:ndvi_t3) |> 
  distinct(id, .keep_all = TRUE)

births_exposure <- births |> 
  left_join(sovi, by = "com") |> 
  left_join(exposure, by = c("id", "com")) 

glimpse(births_exposure)

## IQR Exposure vars 

exposure_vars <- colnames(births_exposure)[38:79]
exposure_vars

iqr_vals <- births_exposure |>
  summarise(across(all_of(exposure_vars), ~ IQR(.x, na.rm = TRUE))) |>
  as.list()
iqr_vals # Calculate IQR per variable 

writexl::write_xlsx(iqr_vals |> data.frame(), paste0(data_out, "Data_IQR_ref_values.xlsx"))

births_exposure <- births_exposure |> 
  mutate(across(all_of(exposure_vars), ~ .x / iqr_vals[[cur_column()]], .names = "{.col}_iqr"))

glimpse(births_exposure)

save(births_exposure, file=paste0(data_out, "births_2010_2020_exposure", ".RData"))

## Join data exposure weeks ----

week_exposure <- week_exposure |> 
  select(-c("name_com", "weeks", "date_nac")) 

births_exposure_weeks <- births |> 
  left_join(sovi, by = "com") |> 
  left_join(week_exposure, by=c("id", "com", "date_start_week_gest", "date_ends_week_gest"), multiple = "all") 

glimpse(births_exposure_weeks)
summary(births_exposure_weeks)

## IQR Exposure vars 

exposure_vars <- c("pm25_krg", "o3_krg", "no2_krg", "pm25_idw", "o3_idw", "no2_idw")
exposure_vars

iqr_vals <- births_exposure_weeks |>
  summarise(across(all_of(exposure_vars), ~ IQR(.x, na.rm = TRUE))) |>
  as.list()
iqr_vals # Calculate IQR per variable 

writexl::write_xlsx(iqr_vals |> data.frame(), paste0(data_out, "Data_IQR_ref_values_weeks.xlsx"))

births_exposure_weeks <- births_exposure_weeks |> 
  mutate(across(all_of(exposure_vars), ~ .x / iqr_vals[[cur_column()]], .names = "{.col}_iqr"))

glimpse(births_exposure_weeks)

save(births_exposure_weeks, file=paste0(data_out, "births_2010_2020_exposure_weeks", ".RData"))
