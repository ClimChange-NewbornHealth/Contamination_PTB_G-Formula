# Code 5: Births exposure data ----
rm(list = ls())

## Settings ----
source("00_Code/0.1 Settings.R")
source("00_Code/0.2 Packages.R")
source("00_Code/0.3 Functions.R")

data_inp <- "01_Data/Output/"
data_out <- "01_Data/Output/"

## Load data (exposure and births) ----

births <- rio::import(paste0(data_out, "births_2010_2020.RData"))
glimpse(births)

exposure <- rio::import(paste0(data_out, "Contamination_Climate_Data_2010_2020.RData"))
glimpse(exposure)

## Births long data -----
births_weeks <- births |> 
  rowwise() |> 
  mutate(week_gest = list(seq.Date(date_start_week_gest, date_ends_week_gest, by = "week"))) |>
  unnest(week_gest) |>
  group_by(id) |>
  mutate(week_gest_num = paste0(abs(weeks - row_number())),  
         week_gest_num = (weeks) - as.numeric(week_gest_num), 
         date_start_week = (week_gest - (7 * abs(week_gest_num - row_number()))) - weeks(1), #(abs(week_gest_num - row_number())),
         date_end_week = week_gest - (7 * abs(week_gest_num - row_number()))
         ) |> # ,(abs(week_gest_num - row_number())
  group_by(id) |> 
  distinct(week_gest_num, .keep_all = TRUE) |> 
  arrange(id, week_gest_num) |> 
  ungroup()

glimpse(births_weeks)

# Check results 
t1 <- births_weeks %>%
  group_by(id) %>% 
  summarise(min=min(week_gest_num), 
            max=max(week_gest_num), 
            n=n(), 
            test=if_else(n==max, 1, 0))

table(t1$test)

## Save new births data ----
save(births_weeks, file=paste0(data_out, "births_2010_2020_weeks_long", ".RData"))

bw_data <- births_weeks |>
  select(id, com, name_com, date_nac, weeks, week_gest_num, date_start_week_gest, date_ends_week_gest)

glimpse(bw_data)
rm(births, births_weeks, t1)

## Pollution data -----

table(unique(exposure$name_com) %in% unique(bw_data$name_com))

exposure <- exposure |> 
  rename(
        no2_idw = no2_idw_pred, 
        no2_krg = no2_ok_pred,
        pm25_idw = pm25_idw_pred,
        pm25_krg = pm25_ok_pred,
        o3_idw = o3_idw_pred,
        o3_krg = o3_ok_pred 
        )

glimpse(exposure)
summary(exposure)

cont_data <- data.table::as.data.table(exposure)
cont_tibble <- dplyr::as_tibble(cont_data)

## Exposure data ---- 

calc_exposure_periods <- function(start_date, end_date, cont_data) {
  cont_data |>
    filter(date >= start_date, date <= end_date) |>
    dplyr::summarise(
      pm25_krg    = mean(pm25_krg,    na.rm = TRUE),
      o3_krg      = mean(o3_krg,      na.rm = TRUE),
      no2_krg     = mean(no2_krg,      na.rm = TRUE),
      pm25_idw    = mean(pm25_idw,    na.rm = TRUE),
      o3_idw      = mean(o3_idw,      na.rm = TRUE),
      no2_idw     = mean(no2_idw,      na.rm = TRUE),
      tad         = mean(TAD,      na.rm = TRUE),
      ndvi        = mean(ndvi,      na.rm = TRUE),
    )
}

calc_all_exposures <- function(start_date, end_date, cont_data) {
  # Full / 30 / 4
  full <- calc_exposure_periods(start_date,           end_date, cont_data)
  d30  <- calc_exposure_periods(end_date - days(30),  end_date, cont_data)
  d4   <- calc_exposure_periods(end_date - days(4),   end_date, cont_data)
  
  # Trimester
  t1_end <- min(start_date + weeks(13), end_date)
  t2_end <- min(start_date + weeks(26), end_date)
  
  tri1 <- calc_exposure_periods(start_date,       t1_end, cont_data)
  tri2 <- calc_exposure_periods(t1_end + days(1), t2_end, cont_data)
  tri3 <- calc_exposure_periods(t2_end + days(1), end_date, cont_data)
  
  # combinar todo en un tibble de salida
  tibble(
    # full
    pm25_krg_full    = full$pm25_krg,
    o3_krg_full      = full$o3_krg,
    no2_krg_full      = full$no2_krg,
    pm25_idw_full    = full$pm25_idw,
    o3_idw_full      = full$o3_idw,
    no2_idw_full      = full$no2_idw,
    tad_full          = full$tad,
    ndvi_full      = full$ndvi,
    
    # 30 días
    pm25_krg_30    = d30$pm25_krg,
    o3_krg_30      = d30$o3_krg,
    no2_krg_30      = d30$no2_krg,
    pm25_idw_30    = d30$pm25_idw,
    o3_idw_30      = d30$o3_idw,
    no2_idw_30      = d30$no2_idw,
    tad_30          = d30$tad,
    ndvi_30      = d30$ndvi,
    
    # 4 días
    pm25_krg_4    = d4$pm25_krg,
    o3_krg_4      = d4$o3_krg,
    no2_krg_4      = d4$no2_krg,
    pm25_idw_4    = d4$pm25_idw,
    o3_idw_4      = d4$o3_idw,
    no2_idw_4      = d4$no2_idw,
    tad_4          = d4$tad,
    ndvi_4      = d4$ndvi,
    
    # Trimester (mean)
    pm25_krg_t1 = tri1$pm25_krg,
    pm25_krg_t2 = tri2$pm25_krg,
    pm25_krg_t3 = tri3$pm25_krg,

    o3_krg_t1   = tri1$o3_krg,
    o3_krg_t2   = tri2$o3_krg,
    o3_krg_t3   = tri3$o3_krg,

    no2_krg_t1   = tri1$no2_krg,
    no2_krg_t2   = tri2$no2_krg,
    no2_krg_t3   = tri3$no2_krg,

    pm25_idw_t1 = tri1$pm25_idw,
    pm25_idw_t2 = tri2$pm25_idw,
    pm25_idw_t3 = tri3$pm25_idw,

    o3_idw_t1   = tri1$o3_idw,
    o3_idw_t2   = tri2$o3_idw,
    o3_idw_t3   = tri3$o3_idw,

    no2_idw_t1   = tri1$no2_idw,
    no2_idw_t2   = tri2$no2_idw,
    no2_idw_t3   = tri3$no2_idw,
    
    tad_t1       = tri1$tad,
    tad_t2       = tri2$tad,
    tad_t3       = tri3$tad,

    ndvi_t1      = tri1$ndvi,
    ndvi_t2      = tri2$ndvi,
    ndvi_t3      = tri3$ndvi
  )
}

calculate_period_metrics_row <- function(row, cont_tbl) {
  setDT(row)
  out <- calc_all_exposures(
    start_date = as.Date(row$date_start_week_gest[1]),
    end_date = as.Date(row$date_ends_week_gest[1]),
    cont_data = cont_tbl
  )
  cbind(row, data.table::as.data.table(out))
}

parts_in_dir <- file.path(data_out, "temp_exposure_bw_parts")
parts_out_gest <- file.path(data_out, "temp_exposure_gest_means_proc")
parts_out_period <- file.path(data_out, "temp_exposure_period_metrics_proc")

num_parts <- 50L
parts(bw_data, path = data_out, folder = "temp_exposure_bw_parts", num_parts = num_parts)

parallel::detectCores()
plan(multisession, workers = parallel::detectCores() - 4)

tic()
process_files_parallel(
  input_directory = parts_in_dir,
  output_directory = parts_out_period,
  cont_data = cont_tibble,
  calc_func = calculate_period_metrics_row
)
bw_period_metrics_expo <- combine_processed_parts(parts_out_period)
if ("I" %in% names(bw_period_metrics_expo)) bw_period_metrics_expo[, I := NULL]
toc()
beepr::beep(8)

tic()
process_files_parallel(
  input_directory = parts_in_dir,
  output_directory = parts_out_gest,
  cont_data = cont_data,
  calc_func = calculate_gest_window_means
)
bw_gest_window_expo <- combine_processed_parts(parts_out_gest)
if ("I" %in% names(bw_gest_window_expo)) bw_gest_window_expo[, I := NULL]
toc()
beepr::beep(8)

glimpse(bw_gest_window_expo)
glimpse(bw_period_metrics_expo)
plan(sequential)

## Check results ----

check_tbl <- bw_period_metrics_expo |>
  dplyr::as_tibble() |>
  select(id, name_com, date_start_week_gest, date_ends_week_gest,
         pm25_krg_full,
         pm25_krg_30,
         pm25_krg_4)
check <- if (nrow(check_tbl) <= 3L) check_tbl else dplyr::slice_sample(check_tbl, n = 3L)

mun <- check$name_com[2]
s <- as.Date(check$date_start_week_gest[2])
e <- as.Date(check$date_ends_week_gest[2])
e30 <- e - days(30) 

cont_data |> filter(# municipio==mun & 
  (date >= s & date <= e)
) |> 
  select(date, pm25_krg) |> 
  summarise(mean = mean(pm25_krg))

cont_data |> filter( # municipio==mun & 
  (date >= e30 & date <= e)
) |> 
  select(date, pm25_krg) |> 
  summarise(mean = mean(pm25_krg))

# Results ok

## Save results ----
save(cont_data, file=paste0(data_out, "series_contamination_pm25_o3_kriging_idw", ".RData"))

save(
  bw_gest_window_expo,
  file = paste0(
    data_out,
    "births_exposure_gest_window_means",
    ".RData"
  )
)

save(
  bw_period_metrics_expo,
  file = paste0(
    data_out,
    "births_exposure_period_metrics_full_d30_d4_tri",
    ".RData"
  )
)
