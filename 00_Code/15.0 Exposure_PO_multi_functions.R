# 15.0 Exposure–perinatal outcome models multicontaminante — funciones compartidas (Cox) ----
#
# Extensión de 11.0: misma estructura Cox (overall + trimestres), pero cada modelo
# del contaminante focal ajusta además por las exposiciones overall o trimestrales
# de los otros contaminantes:
#   - Overall (single_multi): PM2.5_full ~ PM2.5_full + NO2_full + O3_full
#   - Trimestres (t1_t2_t3_multi): PM2.5_t1+t2+t3 ~ PM2.5_t1+t2+t3 + NO2_t1+t2+t3 + O3_t1+t2+t3

PO_MULTI_EXPOSURE_POLLUTANTS <- c("pm25", "no2", "o3")
PO_MULTI_EXPOSURE_TIPO <- "krg"
PO_MULTI_EXPOSURE_PERIODS_SINGLE <- "tot"
PO_MULTI_EXPOSURE_SCALE <- c("raw", "iqr")

po_multi_exposure_var <- function(pollutant, period, scale = "raw") {
  period_suffix <- if (period == "tot") "full" else period
  var <- paste0(pollutant, "_krg_", period_suffix)
  if (identical(scale, "iqr")) paste0(var, "_iqr") else var
}

po_multi_copollutants_for <- function(focal, pollutants = PO_MULTI_EXPOSURE_POLLUTANTS) {
  setdiff(pollutants, focal)
}

po_multi_focal_predictors <- function(
    focal,
    model_type,
    exposure_scale = "raw") {
  if (identical(model_type, "single_multi")) {
    po_multi_exposure_var(focal, "tot", exposure_scale)
  } else if (identical(model_type, "t1_t2_t3_multi")) {
    c(
      po_multi_exposure_var(focal, "t1", exposure_scale),
      po_multi_exposure_var(focal, "t2", exposure_scale),
      po_multi_exposure_var(focal, "t3", exposure_scale)
    )
  } else {
    character(0L)
  }
}

po_multi_all_predictors <- function(
    focal,
    model_type,
    exposure_scale = "raw",
    pollutants = PO_MULTI_EXPOSURE_POLLUTANTS) {

  if (identical(model_type, "single_multi")) {
    vapply(
      pollutants,
      function(p) po_multi_exposure_var(p, "tot", exposure_scale),
      character(1L)
    )
  } else if (identical(model_type, "t1_t2_t3_multi")) {
    unlist(lapply(pollutants, function(p) {
      c(
        po_multi_exposure_var(p, "t1", exposure_scale),
        po_multi_exposure_var(p, "t2", exposure_scale),
        po_multi_exposure_var(p, "t3", exposure_scale)
      )
    }), use.names = FALSE)
  } else {
    character(0L)
  }
}

po_multi_predictor_string <- function(
    focal,
    model_type,
    exposure_scale = "raw",
    pollutants = PO_MULTI_EXPOSURE_POLLUTANTS) {
  paste(
    po_multi_all_predictors(focal, model_type, exposure_scale, pollutants),
    collapse = " + "
  )
}

po_multi_na_result_cox <- function(
    predictor,
    dependent,
    tiempo,
    contaminante,
    tipo,
    model_type,
    adjustment,
    exposure_scale,
    n = 0L) {
  data.frame(
    term = predictor,
    estimate = NA_real_,
    conf.low = NA_real_,
    conf.high = NA_real_,
    log_hr = NA_real_,
    log_hr_conf.low = NA_real_,
    log_hr_conf.high = NA_real_,
    std.error = NA_real_,
    statistic = NA_real_,
    p.value = NA_real_,
    dependent_var = dependent,
    predictor = predictor,
    tiempo = tiempo,
    contaminante = contaminante,
    tipo = tipo,
    model_type = model_type,
    adjustment = adjustment,
    exposure_scale = exposure_scale,
    n = n,
    stringsAsFactors = FALSE
  )
}

fit_po_multi_cox_model <- function(
    dependent,
    predictor,
    tiempo,
    contaminante,
    tipo,
    model_type,
    data,
    control_vars,
    time_var = "weeks",
    time_start = "tstart",
    conf.level = 0.95,
    adjustment = "Adjusted",
    exposure_scale = "raw") {

  all_predictors <- po_multi_all_predictors(
    focal = contaminante,
    model_type = model_type,
    exposure_scale = exposure_scale
  )
  focal_predictors <- po_multi_focal_predictors(
    focal = contaminante,
    model_type = model_type,
    exposure_scale = exposure_scale
  )

  missing_predictors <- all_predictors[!all_predictors %in% names(data)]
  if (length(missing_predictors)) {
    return(po_multi_na_result_cox(
      predictor, dependent, tiempo, contaminante, tipo,
      model_type, adjustment, exposure_scale
    ))
  }

  data_subset <- data |>
    dplyr::filter(
      !is.na(.data[[dependent]]),
      !is.na(.data[[time_var]])
    )

  use_delayed <- !is.null(time_start) &&
    is.character(time_start) &&
    nzchar(time_start) &&
    time_start %in% names(data_subset)
  if (use_delayed) {
    data_subset <- data_subset |>
      dplyr::filter(!is.na(.data[[time_start]])) |>
      dplyr::filter(.data[[time_start]] < .data[[time_var]])
  }

  for (pred in all_predictors) {
    data_subset <- data_subset |> dplyr::filter(!is.na(.data[[pred]]))
  }

  if (nrow(data_subset) < 10L) {
    return(po_multi_na_result_cox(
      predictor, dependent, tiempo, contaminante, tipo,
      model_type, adjustment, exposure_scale, nrow(data_subset)
    ))
  }

  if (identical(adjustment, "Adjusted")) {
    available_controls <- control_vars[control_vars %in% names(data_subset)]
    rhs <- if (length(available_controls)) {
      paste(
        paste(all_predictors, collapse = " + "),
        paste("+", paste(available_controls, collapse = " + "))
      )
    } else {
      paste(all_predictors, collapse = " + ")
    }
  } else {
    rhs <- paste(all_predictors, collapse = " + ")
  }

  surv_lhs <- if (use_delayed) {
    paste0("Surv(", time_start, ", ", time_var, ", ", dependent, ")")
  } else {
    paste0("Surv(", time_var, ", ", dependent, ")")
  }
  fml <- stats::as.formula(paste0(surv_lhs, " ~ ", rhs))

  model_fit <- tryCatch(
    survival::coxph(fml, data = data_subset),
    error = function(e) NULL
  )

  if (is.null(model_fit)) {
    return(po_multi_na_result_cox(
      predictor, dependent, tiempo, contaminante, tipo,
      model_type, adjustment, exposure_scale, nrow(data_subset)
    ))
  }

  tbl <- broom::tidy(
    model_fit,
    exponentiate = TRUE,
    conf.int = TRUE,
    conf.level = conf.level
  )
  tbl_exposure <- tbl[tbl$term %in% focal_predictors, , drop = FALSE]

  if (!nrow(tbl_exposure)) {
    n_obs <- model_fit$n
    rm(model_fit)
    return(po_multi_na_result_cox(
      predictor, dependent, tiempo, contaminante, tipo,
      model_type, adjustment, exposure_scale, n_obs
    ))
  }

  out <- tbl_exposure |>
    dplyr::mutate(
      hr = .data$estimate,
      hr_conf.low = .data$conf.low,
      hr_conf.high = .data$conf.high,
      log_hr = log(.data$estimate),
      log_hr_conf.low = log(.data$conf.low),
      log_hr_conf.high = log(.data$conf.high),
      dependent_var = dependent,
      predictor = predictor,
      tiempo = tiempo,
      contaminante = contaminante,
      tipo = tipo,
      model_type = model_type,
      adjustment = adjustment,
      exposure_scale = exposure_scale,
      n = model_fit$n
    ) |>
    dplyr::select(
      term, estimate, conf.low, conf.high, log_hr, log_hr_conf.low, log_hr_conf.high,
      std.error, statistic, p.value, dependent_var, predictor,
      tiempo, contaminante, tipo, model_type, adjustment, exposure_scale, n
    )

  rm(model_fit)
  invisible(gc())
  out
}

build_po_multi_model_grid <- function(
    dependent_vars,
    pollutants = PO_MULTI_EXPOSURE_POLLUTANTS,
    tipo = PO_MULTI_EXPOSURE_TIPO,
    exposure_scale = PO_MULTI_EXPOSURE_SCALE,
    available_predictors) {

  combinations_single <- expand.grid(
    dependent = dependent_vars,
    tiempo = PO_MULTI_EXPOSURE_PERIODS_SINGLE,
    contaminante = pollutants,
    tipo = tipo,
    model_type = "single_multi",
    adjustment = c("Unadjusted", "Adjusted"),
    exposure_scale = exposure_scale,
    stringsAsFactors = FALSE
  )
  combinations_single <- combinations_single |>
    dplyr::rowwise() |>
    dplyr::mutate(
      predictor = po_multi_predictor_string(
        focal = contaminante,
        model_type = model_type,
        exposure_scale = exposure_scale,
        pollutants = pollutants
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::rowwise() |>
    dplyr::mutate(
      predictors_list = list(po_multi_all_predictors(
        focal = contaminante,
        model_type = model_type,
        exposure_scale = exposure_scale,
        pollutants = pollutants
      )),
      all_exist = all(predictors_list %in% available_predictors)
    ) |>
    dplyr::ungroup() |>
    dplyr::filter(all_exist) |>
    dplyr::select(-predictors_list, -all_exist)

  combinations_t1_t2_t3 <- expand.grid(
    dependent = dependent_vars,
    contaminante = pollutants,
    tipo = tipo,
    model_type = "t1_t2_t3_multi",
    adjustment = c("Unadjusted", "Adjusted"),
    exposure_scale = exposure_scale,
    stringsAsFactors = FALSE
  )
  combinations_t1_t2_t3 <- combinations_t1_t2_t3 |>
    dplyr::rowwise() |>
    dplyr::mutate(
      predictor = po_multi_predictor_string(
        focal = contaminante,
        model_type = model_type,
        exposure_scale = exposure_scale,
        pollutants = pollutants
      ),
      tiempo = "t1_t2_t3"
    ) |>
    dplyr::ungroup() |>
    dplyr::rowwise() |>
    dplyr::mutate(
      predictors_list = list(po_multi_all_predictors(
        focal = contaminante,
        model_type = model_type,
        exposure_scale = exposure_scale,
        pollutants = pollutants
      )),
      all_exist = all(predictors_list %in% available_predictors)
    ) |>
    dplyr::ungroup() |>
    dplyr::filter(all_exist) |>
    dplyr::select(-predictors_list, -all_exist)

  dplyr::bind_rows(combinations_single, combinations_t1_t2_t3)
}

prepare_po_multi_plot_table_data <- function(models_df) {
  models_df |>
    dplyr::filter(
      (model_type == "single_multi" & tiempo == "tot") |
        model_type == "t1_t2_t3_multi"
    ) |>
    dplyr::arrange(
      dependent_var, contaminante, tipo, adjustment,
      model_type, exposure_scale, term
    ) |>
    dplyr::group_by(
      dependent_var, contaminante, tipo, adjustment,
      model_type, exposure_scale
    ) |>
    dplyr::mutate(
      exposure = dplyr::case_when(
        model_type == "single_multi" & tiempo == "tot" ~ "Overall",
        model_type == "t1_t2_t3_multi" & dplyr::row_number() == 1L ~ "Trimester 1",
        model_type == "t1_t2_t3_multi" & dplyr::row_number() == 2L ~ "Trimester 2",
        model_type == "t1_t2_t3_multi" & dplyr::row_number() == 3L ~ "Trimester 3",
        TRUE ~ NA_character_
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::filter(!is.na(exposure)) |>
    dplyr::mutate(
      exposure = factor(
        exposure,
        levels = c("Overall", "Trimester 1", "Trimester 2", "Trimester 3")
      ),
      adjustment = factor(adjustment, levels = c("Unadjusted", "Adjusted"))
    )
}

ensure_po_multi_log_hr_columns <- function(plot_data_cox) {
  if (!"log_hr" %in% names(plot_data_cox)) {
    plot_data_cox <- plot_data_cox |>
      dplyr::mutate(
        log_hr = log(.data$estimate),
        log_hr_conf.low = log(.data$conf.low),
        log_hr_conf.high = log(.data$conf.high)
      )
  }
  plot_data_cox
}

format_po_multi_effect_ci <- function(estimate, conf_low, conf_high) {
  sprintf("%.3f (%.3f-%.3f)", estimate, conf_low, conf_high)
}
