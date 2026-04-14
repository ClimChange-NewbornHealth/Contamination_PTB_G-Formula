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
  "pm25_idw", "o3_idw", "no2_idw"
)

tic()
data_long <- purrr::reduce(
  pollutant_vars,
  .init = births_weeks,
  .f = function(acc, v) {
    out <- compute_lagged_exposure(
      df = acc,
      idbase = id,
      pollulant = v,
      week = week_gest_num
    )
    lag_col <- paste0(v, "_lagged")
    dplyr::bind_cols(acc, dplyr::select(out, dplyr::all_of(lag_col)))
  }
)
toc() # ~ 6 min 

glimpse(data_long)

save(
  data_long,
  file = paste0(data_inp, "births_2010_2020_exposure_weeks_lagged.RData")
)

rm(births_weeks, g_expo_krg, g_expo_idw)

## Cox models (kriging only) ----

births <- rio::import(paste0(data_inp, "births_exposure_period_metrics_full_d30_d4_tri.RData")) |> 
  select(id, tad_full, ndvi_full) |> 
  distinct(id, .keep_all = TRUE)
glimpse(births)

load(paste0(data_inp, "births_2010_2020_exposure_weeks_lagged.RData"))

data_long <- data_long |> 
  left_join(births, by = "id") |> 
  mutate(
    month_week1 = factor(month_week1), 
    year_week1 = factor(year_week1),
    covid = factor(covid)
  )

control_vars <- c(
  "sex", "age_group_mom", "educ_group_mom", "job_group_mom",
  "age_group_dad", "educ_group_dad", "job_group_dad",
  "month_week1", "year_week1", "covid", "vulnerability", 
  "tad_full", "ndvi_full"
)

dependent_var <- "birth_preterm"
weeks_analysis <- 1:44
lag_weeks <- 2:44

# Build wide inputs for one contaminant with weekly exposure and covariates
build_wide_pollutant <- function(df, pollutant, weeks_keep = 1:37) {
  lag_col <- paste0(pollutant, "_lagged")
  
  wide_one <- df |>
    dplyr::select(id, week_gest_num, all_of(pollutant), all_of(lag_col)) |>
    dplyr::filter(week_gest_num %in% weeks_keep) |>
    tidyr::pivot_wider(
      names_from = week_gest_num,
      values_from = c(all_of(pollutant), all_of(lag_col)),
      names_glue = "{.value}_{week_gest_num}"
    ) |>
    dplyr::rename_with(~ stringr::str_replace(.x, paste0("^", pollutant, "_"), "exposicion_")) |>
    dplyr::rename_with(~ stringr::str_replace(.x, paste0("^", lag_col, "_"), "exposicion_lagged_"))
  
  wide_one
}

# Build modeling base (one row per id)
base_vars <- c("id", "weeks", dependent_var, control_vars)
data_base <- data_long |>
  dplyr::select(all_of(base_vars)) |>
  dplyr::distinct(id, .keep_all = TRUE) 

glimpse(data_base)

### Kriging analysis ----

krg_contaminants <- c("pm25_krg", "o3_krg", "no2_krg") 

# Clean memory 
invisible(gc())                             # clean RAM
while (dev.cur() > 1) dev.off()             # close plots
closeAllConnections()            # close connections    

results_cox <- list()
wide_mat <- list()

tic("Cox adjusted models (kriging)")
for (contam in krg_contaminants) {
  wide_contam <- build_wide_pollutant(data_long, pollutant = contam, weeks_keep = weeks_analysis)
  
  data_model <- data_base |>
    dplyr::left_join(wide_contam, by = "id") |>
    dplyr::mutate(tstart = 28)
  
  results_cox[[contam]] <- data.frame()
  
  for (w in weeks_analysis) {
    exp_var <- paste0("exposicion_", w)
    lag_var <- paste0("exposicion_lagged_", w)
    predictor_terms <- c(
      exp_var,
      lag_var[lag_var %in% paste0("exposicion_lagged_", lag_weeks)]
    )
    predictor <- paste(predictor_terms, collapse = " + ")
    
    tbl_cox <- fit_cox_model(
      dependent = dependent_var,
      predictor = predictor,
      tiempo = paste0("w", w),
      contaminante = contam,
      tipo = "krg",
      data = data_model,
      time_var = "weeks",
      time_start = "tstart"
    )
    
    if (nrow(tbl_cox) > 0) {
      tbl_cox <- tbl_cox |> 
        dplyr::mutate(
          week = w,
          exposure_term = exp_var,
          is_exposure = term == exp_var,
          is_lag = term == lag_var
        )
      results_cox[[contam]] <- dplyr::bind_rows(results_cox[[contam]], tbl_cox)
      wide_mat[[contam]] <- wide_contam 
    }
  }
}
toc() # ~ 12-15 min
beep(8)

dlm_cox_krg_results <- list(
  results_cox = results_cox,
  dependent_var = dependent_var,
  contaminants = krg_contaminants,
  weeks_analysis = weeks_analysis,
  lag_weeks = lag_weeks,
  control_vars = control_vars
)

save(
  dlm_cox_krg_results,
  file = paste0(data_out_model, "DLM_cox_krg_results.RData")
)

# Table with effects

effects_table <- purrr::imap(
  results_cox,
  function(df, contam_name) {
    df |>
      dplyr::filter(is_exposure, week %in% 1:37) |>
      dplyr::transmute(
        week = as.integer(week),
        !!paste0("hr_", contam_name) := sprintf("%.3f", hr),
        !!paste0("ic_left_", contam_name) := sprintf("%.3f", conf.low),
        !!paste0("ic_right_", contam_name) := sprintf("%.3f", conf.high)
      )
  }
) |>
  purrr::reduce(dplyr::full_join, by = "week") |>
  dplyr::arrange(week)

effects_table

write.xlsx(
  effects_table,
  file = paste0(data_out_model, "DLM_cox_krg_effects_table.xlsx")
)

write.xlsx(
  dlm_cox_krg_results,
  file = paste0(data_out_model, "DLM_cox_krg_effects_complete.xlsx")
)

### IDW analysis ----

idw_contaminants <- c("pm25_idw", "o3_idw", "no2_idw") 

invisible(gc())                             # clean RAM
while (dev.cur() > 1) dev.off()             # close plots
closeAllConnections()            # close connections    

results_cox <- list()
wide_mat <- list()

tic("Cox adjusted models (IDW)")
for (contam in idw_contaminants) {
  wide_contam <- build_wide_pollutant(data_long, pollutant = contam, weeks_keep = weeks_analysis)
  
  data_model <- data_base |>
    dplyr::left_join(wide_contam, by = "id") |>
    dplyr::mutate(tstart = 28)
  
  results_cox[[contam]] <- data.frame()
  
  for (w in weeks_analysis) {
    exp_var <- paste0("exposicion_", w)
    lag_var <- paste0("exposicion_lagged_", w)
    predictor_terms <- c(
      exp_var,
      lag_var[lag_var %in% paste0("exposicion_lagged_", lag_weeks)]
    )
    predictor <- paste(predictor_terms, collapse = " + ")
    
    tbl_cox <- fit_cox_model(
      dependent = dependent_var,
      predictor = predictor,
      tiempo = paste0("w", w),
      contaminante = contam,
      tipo = "idw",
      data = data_model,
      time_var = "weeks",
      time_start = "tstart"
    )
    
    if (nrow(tbl_cox) > 0) {
      tbl_cox <- tbl_cox |> 
        dplyr::mutate(
          week = w,
          exposure_term = exp_var,
          is_exposure = term == exp_var,
          is_lag = term == lag_var
        )
      results_cox[[contam]] <- dplyr::bind_rows(results_cox[[contam]], tbl_cox)
      wide_mat[[contam]] <- wide_contam 
    }
  }
}
toc() # ~ 12-15 min
beep(8)

dlm_cox_idw_results <- list(
  results_cox = results_cox,
  dependent_var = dependent_var,
  contaminants = idw_contaminants,
  weeks_analysis = weeks_analysis,
  lag_weeks = lag_weeks,
  control_vars = control_vars
)

save(
  dlm_cox_idw_results,
  file = paste0(data_out_model, "DLM_cox_idw_results.RData")
)

# Table with effects

effects_table <- purrr::imap(
  results_cox,
  function(df, contam_name) {
    df |>
      dplyr::filter(is_exposure, week %in% 1:37) |>
      dplyr::transmute(
        week = as.integer(week),
        !!paste0("hr_", contam_name) := sprintf("%.3f", hr),
        !!paste0("ic_left_", contam_name) := sprintf("%.3f", conf.low),
        !!paste0("ic_right_", contam_name) := sprintf("%.3f", conf.high)
      )
  }
) |>
  purrr::reduce(dplyr::full_join, by = "week") |>
  dplyr::arrange(week)

effects_table

write.xlsx(
  effects_table,
  file = paste0(data_out_model, "DLM_cox_idw_effects_table.xlsx")
)

write.xlsx(
  dlm_cox_idw_results,
  file = paste0(data_out_model, "DLM_cox_idw_effects_complete.xlsx")
)


