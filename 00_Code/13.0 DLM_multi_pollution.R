# 13.0 DLM multicontaminante (ajuste por co-exposición semanal) -----
#
# Extensión de 9.0 DLM_pollution.R: misma estructura, covariables, semanas,
# lag ponderado y métricas. La única diferencia es que cada modelo DLM del
# contaminante principal ajusta además por la exposición semanal (semana w)
# de los otros dos contaminantes:
#   - PM2.5  ~ PM2.5 DLM + NO2_w + O3_w
#   - NO2    ~ NO2 DLM   + PM2.5_w + O3_w
#   - O3     ~ O3 DLM    + PM2.5_w + NO2_w
#
# Requisito previo: 9.0 (genera births_2010_2020_exposure_weeks_lagged.RData)
# Salida: 02_Output/Models/DLM_multi_cox_{krg,idw}_*

rm(list = ls())

## Settings ----
source("00_Code/0.1 Settings.R")
source("00_Code/0.2 Packages.R")
source("00_Code/0.3 Functions.R")

data_inp <- "01_Data/Output/"
data_out_model <- "02_Output/Models/"

lagged_path <- paste0(data_inp, "births_2010_2020_exposure_weeks_lagged.RData")
if (!file.exists(lagged_path)) {
  stop(
    "No se encontró ", lagged_path,
    ". Ejecute primero 00_Code/9.0 DLM_pollution.R para generar los lags DLM."
  )
}

## Load data ----

births <- rio::import(paste0(data_inp, "births_exposure_period_metrics_full_d30_d4_tri.RData")) |>
  dplyr::select(id, ndvi_full) |>
  dplyr::distinct(id, .keep_all = TRUE)

load(lagged_path)

data_long <- data_long |>
  dplyr::left_join(births, by = "id") |>
  dplyr::mutate(
    month_week1 = factor(month_week1),
    year_week1 = factor(year_week1),
    covid = factor(covid)
  )

control_vars <- c(
  "sex", "age_group_mom", "educ_group_mom", "job_group_mom",
  "age_group_dad", "educ_group_dad", "job_group_dad",
  "month_week1", "year_week1", "covid", "vulnerability", "tad",
  "ndvi_full"
)

dependent_var <- "birth_preterm"
weeks_analysis <- 1:44
lag_weeks <- 2:44

pollutant_short <- function(pollutant) {
  sub("_(krg|idw)$", "", pollutant)
}

copollutants_for <- function(primary, all_contaminants) {
  setdiff(all_contaminants, primary)
}

# Wide del contaminante principal (exposición semanal + lag ponderado DLM)
build_wide_pollutant <- function(df, pollutant, weeks_keep = weeks_analysis) {
  lag_col <- paste0(pollutant, "_lagged")

  df |>
    dplyr::select(id, week_gest_num, dplyr::all_of(pollutant), dplyr::all_of(lag_col)) |>
    dplyr::filter(week_gest_num %in% weeks_keep) |>
    tidyr::pivot_wider(
      names_from = week_gest_num,
      values_from = c(dplyr::all_of(pollutant), dplyr::all_of(lag_col)),
      names_glue = "{.value}_{week_gest_num}"
    ) |>
    dplyr::rename_with(~ stringr::str_replace(.x, paste0("^", pollutant, "_"), "exposicion_")) |>
    dplyr::rename_with(~ stringr::str_replace(.x, paste0("^", lag_col, "_"), "exposicion_lagged_"))
}

# Wide semanal de co-contaminantes (solo semana w; sin lag DLM)
build_wide_copollutant <- function(df, pollutant, weeks_keep = weeks_analysis) {
  short <- pollutant_short(pollutant)

  df |>
    dplyr::select(id, week_gest_num, dplyr::all_of(pollutant)) |>
    dplyr::filter(week_gest_num %in% weeks_keep) |>
    tidyr::pivot_wider(
      names_from = week_gest_num,
      values_from = dplyr::all_of(pollutant),
      names_glue = "{.value}_{week_gest_num}"
    ) |>
    dplyr::rename_with(
      ~ stringr::str_replace(.x, paste0("^", pollutant, "_"), paste0("copoll_", short, "_"))
    )
}

build_wide_weekly_var <- function(df, varname, weeks_keep = weeks_analysis) {
  df |>
    dplyr::select(id, week_gest_num, dplyr::all_of(varname)) |>
    dplyr::filter(week_gest_num %in% weeks_keep) |>
    tidyr::pivot_wider(
      names_from = week_gest_num,
      values_from = dplyr::all_of(varname),
      names_glue = "{.value}_{week_gest_num}"
    )
}

base_vars <- c("id", "weeks", dependent_var, control_vars)
data_base <- data_long |>
  dplyr::select(dplyr::all_of(base_vars)) |>
  dplyr::distinct(id, .keep_all = TRUE)

wide_tad <- build_wide_weekly_var(data_long, "tad", weeks_keep = weeks_analysis)
wide_climate <- wide_tad

copollutant_terms_week <- function(primary, all_contaminants, week) {
  copolls <- copollutants_for(primary, all_contaminants)
  paste0("copoll_", pollutant_short(copolls), "_", week)
}

build_multicontaminant_model_data <- function(
    primary,
    all_contaminants,
    data_long,
    data_base,
    wide_climate,
    weeks_keep = weeks_analysis) {

  wide_primary <- build_wide_pollutant(data_long, pollutant = primary, weeks_keep = weeks_keep)
  copolls <- copollutants_for(primary, all_contaminants)

  wide_copolls <- purrr::reduce(
    copolls,
    .init = NULL,
    .f = function(acc, copoll) {
      wide_one <- build_wide_copollutant(data_long, pollutant = copoll, weeks_keep = weeks_keep)
      if (is.null(acc)) {
        wide_one
      } else {
        dplyr::full_join(acc, wide_one, by = "id")
      }
    }
  )

  data_base |>
    dplyr::left_join(wide_primary, by = "id") |>
    dplyr::left_join(wide_copolls, by = "id") |>
    dplyr::left_join(wide_climate, by = "id") |>
    dplyr::mutate(tstart = 28)
}

run_multicontaminant_dlm <- function(
    contaminants,
    exposure_type,
    data_long,
    data_base,
    wide_climate) {

  results_cox <- list()
  wide_mat <- list()

  for (contam in contaminants) {
    data_model <- build_multicontaminant_model_data(
      primary = contam,
      all_contaminants = contaminants,
      data_long = data_long,
      data_base = data_base,
      wide_climate = wide_climate
    )

    results_cox[[contam]] <- data.frame()

    for (w in weeks_analysis) {
      exp_var <- paste0("exposicion_", w)
      lag_var <- paste0("exposicion_lagged_", w)
      tad_var <- paste0("tad_", w)
      co_terms <- copollutant_terms_week(contam, contaminants, w)

      predictor_terms <- c(
        exp_var,
        lag_var[lag_var %in% paste0("exposicion_lagged_", lag_weeks)],
        co_terms,
        tad_var
      )
      predictor <- paste(predictor_terms, collapse = " + ")

      tbl_cox <- fit_cox_model(
        dependent = dependent_var,
        predictor = predictor,
        tiempo = paste0("w", w),
        contaminante = contam,
        tipo = exposure_type,
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
            is_lag = term == lag_var,
            is_copollutant = term %in% co_terms
          )
        results_cox[[contam]] <- dplyr::bind_rows(results_cox[[contam]], tbl_cox)
        wide_mat[[contam]] <- data_model
      }
    }
  }

  list(
    results_cox = results_cox,
    wide_mat = wide_mat
  )
}

save_multicontaminant_results <- function(
    fit_obj,
    contaminants,
    exposure_type,
    file_prefix) {

  out <- list(
    results_cox = fit_obj$results_cox,
    dependent_var = dependent_var,
    contaminants = contaminants,
    weeks_analysis = weeks_analysis,
    lag_weeks = lag_weeks,
    control_vars = control_vars,
    model_type = "multicontaminant_dlm",
    copollutant_adjustment = "weekly concurrent exposure at week w for the other two pollutants"
  )

  rdata_path <- paste0(data_out_model, file_prefix, "_results.RData")
  save(out, file = rdata_path)

  effects_table <- purrr::imap(
    fit_obj$results_cox,
    function(df, contam_name) {
      df |>
        dplyr::filter(.data$is_exposure, week %in% 1:37) |>
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

  write.xlsx(
    effects_table,
    file = paste0(data_out_model, file_prefix, "_effects_table.xlsx")
  )

  write.xlsx(
    out,
    file = paste0(data_out_model, file_prefix, "_effects_complete.xlsx")
  )

  invisible(out)
}

## Cox models multicontaminante (kriging) ----

krg_contaminants <- c("pm25_krg", "o3_krg", "no2_krg")

invisible(gc())
while (dev.cur() > 1) dev.off()
closeAllConnections()

tic("Cox adjusted multicontaminant models (kriging)")
fit_krg <- run_multicontaminant_dlm(
  contaminants = krg_contaminants,
  exposure_type = "krg",
  data_long = data_long,
  data_base = data_base,
  wide_climate = wide_climate
)
toc()
beep(8)

dlm_multi_cox_krg_results <- save_multicontaminant_results(
  fit_obj = fit_krg,
  contaminants = krg_contaminants,
  exposure_type = "krg",
  file_prefix = "DLM_multi_cox_krg"
)

rm(fit_krg)
invisible(gc())

## Cox models multicontaminante (IDW) ----

idw_contaminants <- c("pm25_idw", "o3_idw", "no2_idw")

invisible(gc())
while (dev.cur() > 1) dev.off()
closeAllConnections()

tic("Cox adjusted multicontaminant models (IDW)")
fit_idw <- run_multicontaminant_dlm(
  contaminants = idw_contaminants,
  exposure_type = "idw",
  data_long = data_long,
  data_base = data_base,
  wide_climate = wide_climate
)
toc()
beep(8)

dlm_multi_cox_idw_results <- save_multicontaminant_results(
  fit_obj = fit_idw,
  contaminants = idw_contaminants,
  exposure_type = "idw",
  file_prefix = "DLM_multi_cox_idw"
)

rm(fit_idw, data_long, data_base, wide_climate, wide_tad, births)
invisible(gc())

message("Modelos multicontaminante guardados en ", data_out_model)
