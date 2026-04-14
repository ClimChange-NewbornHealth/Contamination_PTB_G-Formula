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

g_expo_krg <- births_weeks |> 
  group_by(week_gest_num) |> 
  summarise(
    pm25_mean = mean(pm25_krg, na.rm = TRUE),
    pm25_p50 = median(pm25_krg, na.rm = TRUE),
    pm25_sd = sd(pm25_krg, na.rm = TRUE),

    no2_mean = mean(no2_krg, na.rm = TRUE),
    no2_p50 = median(no2_krg, na.rm = TRUE),
    no2_sd = sd(no2_krg, na.rm = TRUE),

    o3_mean = mean(o3_krg, na.rm = TRUE),
    o3_p50 = median(o3_krg, na.rm = TRUE),
    o3_sd = sd(o3_krg, na.rm = TRUE)
  ) |> 
  ungroup() |> 
  mutate(across(where(is.numeric), ~ format(round(.x, 2), nsmall = 2, decimal.mark = ".")))

g_expo_idw <- births_weeks |> 
  group_by(week_gest_num) |> 
  summarise(
    pm25_mean = mean(pm25_idw, na.rm = TRUE),
    pm25_p50 = median(pm25_idw, na.rm = TRUE),
    pm25_sd = sd(pm25_idw, na.rm = TRUE),

    no2_mean = mean(no2_idw, na.rm = TRUE),
    no2_p50 = median(no2_idw, na.rm = TRUE),
    no2_sd = sd(no2_idw, na.rm = TRUE),

    o3_mean = mean(o3_idw, na.rm = TRUE),
    o3_p50 = median(o3_idw, na.rm = TRUE),
    o3_sd = sd(o3_idw, na.rm = TRUE)
  ) |> 
  ungroup() |> 
  mutate(across(where(is.numeric), ~ format(round(.x, 2), nsmall = 2, decimal.mark = ".")))

# Save results
write.xlsx(
  list(
    krg = g_expo_krg, 
    idw = g_expo_idw
  ), 
  paste0(data_out, "Descriptives_exposure_stats_time.xlsx"))

## Weighted DLM ----

# Data are in long format here. For each idbase and contaminant-tipo, we compute
# exposicion_lagged = weighted sum of past exposures, with weights = 1 / (current_week - past_week).
# This section computes lagged weighted exposures for each pollutant.

compute_lagged_exposure <- function(df, idbase, pollulant, week) {
  poll_quo <- rlang::enquo(pollulant)
  week_quo <- rlang::enquo(week)
  
  poll_name <- if (rlang::quo_is_symbol(poll_quo)) {
    candidate <- rlang::as_name(poll_quo)
    if (candidate %in% names(df)) {
      candidate
    } else {
      evaluated <- rlang::eval_tidy(poll_quo)
      if (is.character(evaluated) && length(evaluated) == 1) {
        evaluated
      } else {
        rlang::abort("pollulant must resolve to a valid column name.")
      }
    }
  } else if (is.character(pollulant) && length(pollulant) == 1) {
    pollulant
  } else {
    rlang::abort("pollulant must be a column name or a single string.")
  }
  
  lag_name <- paste0(poll_name, "_lagged")
  
  df |>
    dplyr::mutate(
      .id_tmp = {{ idbase }},
      .poll_tmp = .data[[poll_name]],
      .week_tmp = !!week_quo
    ) |>
    dplyr::arrange(.id_tmp, .week_tmp) |>
    dplyr::group_by(.id_tmp) |>
    dplyr::mutate(
      # Weighted cumulative exposure using previous weeks only
      .lag_tmp = purrr::map_dbl(dplyr::row_number(), function(i) {
        if (is.na(.week_tmp[i]) || .week_tmp[i] == 0) return(NA_real_)
        past_rows <- which(.week_tmp < .week_tmp[i])
        if (length(past_rows) == 0) return(NA_real_)
        weights <- 1 / (.week_tmp[i] - .week_tmp[past_rows])
        exposures <- .poll_tmp[past_rows]
        sum(weights * exposures, na.rm = TRUE)
      })
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-.id_tmp, -.poll_tmp, -.week_tmp) |>
    dplyr::rename(!!lag_name := .lag_tmp)
}

pollutant_vars <- c(
  "pm25_krg", "o3_krg", "no2_krg",
  "pm25_krg_iqr", "o3_krg_iqr", "no2_krg_iqr"
)

data_long <- purrr::reduce(
  pollutant_vars,
  .init = births_weeks,
  .f = function(acc, v) {
    out <- compute_lagged_exposure(
      df = acc,
      idbase = id,
      pollulant = v,
      week = week_gest_num - 1
    )
    lag_col <- paste0(v, "_lagged")
    dplyr::bind_cols(acc, dplyr::select(out, dplyr::all_of(lag_col)))
  }
)

glimpse(data_long)