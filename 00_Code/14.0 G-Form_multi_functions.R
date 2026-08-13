# 14.0 G-Formula multicontaminante — funciones ----
#
# Extensión de 10.1 G-Form_functions.R: misma lógica g-computation, covariables,
# bootstrap, heatmap y paralelismo. Cada modelo Cox semanal del contaminante
# principal ajusta además por la exposición semanal concurrente (semana w) de
# los otros dos contaminantes (como 13.0 DLM_multi).
#
# Salidas en 02_Output/G-Form-Multi/ (no sobrescribe G-Form).

source("00_Code/10.1 G-Form_functions.R")

## Configuración multicontaminante ----

GFORM_MULTI_DATA_OUT <- "02_Output/G-Form-Multi/"

GFORM_MULTI_KRG_CONTAMINANTS <- c("pm25_krg", "no2_krg", "o3_krg")

GFORM_MULTI_INTERVENTION_ORDER <- c(
  "pm25_krg_pct20",
  "no2_krg_pct20",
  "o3_krg_pct20"
)

GFORM_MULTI_INTERVENTION_REGISTRY <- list(
  pm25_krg_pct20 = c(
    GFORM_INTERVENTION_REGISTRY$pm25_krg_pct20,
    list(
      description = paste0(
        GFORM_INTERVENTION_REGISTRY$pm25_krg_pct20$description,
        " [multicontaminante: ajuste semanal NO2 + O3]"
      )
    )
  ),
  no2_krg_pct20 = c(
    GFORM_INTERVENTION_REGISTRY$no2_krg_pct20,
    list(
      description = paste0(
        GFORM_INTERVENTION_REGISTRY$no2_krg_pct20$description,
        " [multicontaminante: ajuste semanal PM2.5 + O3]"
      )
    )
  ),
  o3_krg_pct20 = c(
    GFORM_INTERVENTION_REGISTRY$o3_krg_pct20,
    list(
      description = paste0(
        GFORM_INTERVENTION_REGISTRY$o3_krg_pct20$description,
        " [multicontaminante: ajuste semanal PM2.5 + NO2]"
      )
    )
  )
)

GFORM_MULTI_ACTIVE <- new.env(parent = emptyenv())

.gform_build_cox_model_frame <- build_cox_model_frame
.gform_fit_cox_one_week <- fit_cox_one_week
.gform_fit_natural_course_models <- fit_natural_course_models
.gform_model_cache_path <- gform_model_cache_path
.gform_natural_course_cache_key <- gform_natural_course_cache_key

gform_multi_pollutant_short <- function(pollutant) {
  sub("_krg$", "", pollutant)
}

gform_multi_copollutants_for <- function(primary) {
  setdiff(GFORM_MULTI_KRG_CONTAMINANTS, primary)
}

gform_multi_copollutant_terms_week <- function(primary, week) {
  copolls <- gform_multi_copollutants_for(primary)
  paste0("copoll_", gform_multi_pollutant_short(copolls), "_", week)
}

build_wide_copollutant_raw <- function(
    df,
    pollutant,
    weeks_keep = GFORM_DEFAULTS$weeks_exposure) {

  short <- gform_multi_pollutant_short(pollutant)
  wide <- build_wide_raw_exposure(df, pollutant, weeks_keep = weeks_keep)
  wide <- as.data.frame(wide)
  names(wide) <- sub(
    paste0("^", pollutant, "_"),
    paste0("copoll_", short, "_"),
    names(wide)
  )
  wide
}

build_wide_copollutants_obs <- function(
    data_long,
    primary,
    weeks_keep = GFORM_DEFAULTS$weeks_exposure) {

  copolls <- gform_multi_copollutants_for(primary)
  wide_copolls <- purrr::reduce(
    copolls,
    .init = NULL,
    .f = function(acc, copoll) {
      wide_one <- build_wide_copollutant_raw(
        data_long,
        pollutant = copoll,
        weeks_keep = weeks_keep
      )
      if (is.null(acc)) {
        wide_one
      } else {
        dplyr::full_join(acc, wide_one, by = "id")
      }
    }
  )
  data.table::as.data.table(wide_copolls)
}

gform_multi_intervention_keys <- function() {
  missing <- setdiff(
    GFORM_MULTI_INTERVENTION_ORDER,
    names(GFORM_MULTI_INTERVENTION_REGISTRY)
  )
  if (length(missing)) {
    stop(
      "GFORM_MULTI_INTERVENTION_ORDER contiene IDs ausentes: ",
      paste(missing, collapse = ", ")
    )
  }
  GFORM_MULTI_INTERVENTION_ORDER
}

print_gform_multi_intervention_menu <- function(selected = NULL) {
  keys <- gform_multi_intervention_keys()
  message("Intervenciones multicontaminante (14.2):")
  for (i in seq_along(keys)) {
    spec <- GFORM_MULTI_INTERVENTION_REGISTRY[[keys[i]]]
    marker <- if (!is.null(selected) && i == selected) "  <-- seleccionada" else ""
    message(sprintf(
      "  %2d. %s — %s%s",
      i, spec$intervention_id, spec$description, marker
    ))
  }
  invisible(keys)
}

resolve_gform_multi_intervention <- function(n) {
  keys <- gform_multi_intervention_keys()
  n <- as.integer(n[[1L]])
  if (length(n) != 1L || is.na(n) || n < 1L || n > length(keys)) {
    stop(
      "intervention_number debe ser un entero entre 1 y ", length(keys),
      ". Ejecute print_gform_multi_intervention_menu()."
    )
  }
  keys[[n]]
}

gform_multi_set_active <- function(primary_pollutant, wide_copollutants) {
  GFORM_MULTI_ACTIVE$primary_pollutant <- primary_pollutant
  GFORM_MULTI_ACTIVE$wide_copollutants <- wide_copollutants
  invisible(NULL)
}

gform_multi_clear_active <- function() {
  GFORM_MULTI_ACTIVE$primary_pollutant <- NULL
  GFORM_MULTI_ACTIVE$wide_copollutants <- NULL
  invisible(NULL)
}

## Frame Cox slim (solo semanas de riesgo con modelo Cox) ----

gform_multi_cox_frame_weeks <- function(
    risk_weeks_vec = GFORM_DEFAULTS$risk_weeks,
    risk_entry_week = GFORM_DEFAULTS$risk_entry_week) {
  risk_weeks_vec[risk_weeks_vec >= risk_entry_week]
}

gform_multi_cox_exposure_columns <- function(weeks, lag_weeks = GFORM_DEFAULTS$lag_weeks) {
  unique(c(
    paste0("exposicion_", weeks),
    paste0("exposicion_lagged_", intersect(weeks, lag_weeks))
  ))
}

gform_multi_cox_copoll_columns <- function(primary, weeks) {
  if (is.null(primary)) return(character(0))
  unlist(
    lapply(weeks, gform_multi_copollutant_terms_week, primary = primary),
    use.names = FALSE
  )
}

gform_multi_cox_tad_columns <- function(weeks) {
  paste0("tad_", weeks)
}

gform_multi_slim_wide_table <- function(dt, keep_cols, id_col = "id") {
  if (!data.table::is.data.table(dt)) {
    dt <- data.table::as.data.table(dt)
  }
  keep <- unique(c(id_col, intersect(keep_cols, names(dt))))
  dt[, ..keep]
}

gform_multi_cox_vars_for_week <- function(
    rw,
    control_vars,
    lag_weeks,
    dependent_var,
    risk_entry_week) {

  if (rw < risk_entry_week) return(character(0))

  exp_var <- paste0("exposicion_", rw)
  lag_var <- paste0("exposicion_lagged_", rw)
  tad_var <- paste0("tad_", rw)
  co_terms <- character(0)
  primary <- GFORM_MULTI_ACTIVE$primary_pollutant
  if (!is.null(primary)) {
    co_terms <- gform_multi_copollutant_terms_week(primary, rw)
  }

  pred_terms <- c(
    exp_var,
    lag_var[lag_var %in% paste0("exposicion_lagged_", lag_weeks)],
    co_terms,
    tad_var
  )

  c("id", "weeks", "tstart", dependent_var, pred_terms, control_vars)
}

## Overrides Cox multicontaminante ----

build_cox_model_frame <- function(
    data_base,
    wide_exposicion,
    wide_tad_obs,
    control_vars,
    dependent_var,
    risk_entry_week) {

  weeks_needed <- gform_multi_cox_frame_weeks(
    risk_weeks_vec = GFORM_DEFAULTS$risk_weeks,
    risk_entry_week = risk_entry_week
  )
  lag_weeks <- GFORM_DEFAULTS$lag_weeks

  wide_exposicion <- gform_multi_slim_wide_table(
    wide_exposicion,
    keep_cols = gform_multi_cox_exposure_columns(weeks_needed, lag_weeks)
  )
  wide_tad_obs <- gform_multi_slim_wide_table(
    wide_tad_obs,
    keep_cols = gform_multi_cox_tad_columns(weeks_needed)
  )

  wide_copoll <- GFORM_MULTI_ACTIVE$wide_copollutants
  if (!is.null(wide_copoll)) {
    primary <- GFORM_MULTI_ACTIVE$primary_pollutant
    wide_copoll <- gform_multi_slim_wide_table(
      wide_copoll,
      keep_cols = gform_multi_cox_copoll_columns(primary, weeks_needed)
    )
    data.table::setkey(wide_copoll, id)
  }

  out <- .gform_build_cox_model_frame(
    data_base = data_base,
    wide_exposicion = wide_exposicion,
    wide_tad_obs = wide_tad_obs,
    control_vars = control_vars,
    dependent_var = dependent_var,
    risk_entry_week = risk_entry_week
  )

  if (!is.null(wide_copoll)) {
    data.table::setkey(out, id)
    out <- wide_copoll[out, on = "id"]
  }
  out
}

fit_cox_one_week <- function(
    rw,
    data_model_df,
    control_vars,
    lag_weeks,
    dependent_var,
    risk_entry_week) {

  vars_needed <- gform_multi_cox_vars_for_week(
    rw, control_vars, lag_weeks, dependent_var, risk_entry_week
  )
  if (!length(vars_needed)) return(NULL)

  pred_terms <- setdiff(vars_needed, c("id", "weeks", "tstart", dependent_var, control_vars))

  rhs <- paste(c(pred_terms, control_vars), collapse = " + ")
  fml <- stats::as.formula(paste0(
    "Surv(tstart, weeks, ", dependent_var, ") ~ ", rhs
  ))

  if (data.table::is.data.table(data_model_df)) {
    model_df <- as.data.frame(data_model_df[, ..vars_needed])
  } else {
    model_df <- data_model_df[, vars_needed, drop = FALSE]
  }
  model_df <- stats::na.omit(model_df)

  if (nrow(model_df) < 50L) {
    warning("Muy pocas filas para la semana ", rw, "; modelo Cox omitido.")
    return(NULL)
  }

  cox_fit <- tryCatch(
    survival::coxph(fml, data = model_df, x = FALSE, y = FALSE),
    error = function(e) NULL
  )
  if (is.null(cox_fit)) {
    warning("Cox no convergió o error en semana ", rw)
    return(NULL)
  }

  slim_coxph_model_store_entry(list(
    model = cox_fit,
    formula = fml,
    id_order = model_df$id,
    pred_terms = pred_terms,
    control_vars = control_vars,
    baseline_hazard = survival::basehaz(cox_fit, centered = FALSE),
    risk_week = rw
  ))
}

fit_natural_course_models <- function(
    data_base,
    wide_exposicion_natural,
    wide_tad_obs,
    risk_weeks_vec,
    control_vars = GFORM_DEFAULTS$control_vars,
    lag_weeks = GFORM_DEFAULTS$lag_weeks,
    dependent_var = GFORM_DEFAULTS$dependent_var,
    risk_entry_week = GFORM_DEFAULTS$risk_entry_week,
    parallel = TRUE) {

  cox_frame_natural <- build_cox_model_frame(
    data_base = data_base,
    wide_exposicion = wide_exposicion_natural,
    wide_tad_obs = wide_tad_obs,
    control_vars = control_vars,
    dependent_var = dependent_var,
    risk_entry_week = risk_entry_week
  )
  if (!data.table::is.data.table(cox_frame_natural)) {
    cox_frame_natural <- data.table::as.data.table(cox_frame_natural)
  }
  data.table::setkey(cox_frame_natural, id)

  weeks_to_fit <- risk_weeks_vec[risk_weeks_vec >= risk_entry_week]

  fit_one <- function(rw) {
    vars_needed <- gform_multi_cox_vars_for_week(
      rw, control_vars, lag_weeks, dependent_var, risk_entry_week
    )
    week_df <- cox_frame_natural[, ..vars_needed]
    fit_cox_one_week(
      rw,
      week_df,
      control_vars,
      lag_weeks,
      dependent_var,
      risk_entry_week
    )
  }

  if (parallel && length(weeks_to_fit) > 1L && requireNamespace("furrr", quietly = TRUE)) {
    gform_setup_parallel(task = "cox")
    fitted <- furrr::future_map(
      weeks_to_fit,
      fit_one,
      .options = gform_furrr_options()
    )
  } else {
    fitted <- lapply(weeks_to_fit, fit_one)
  }

  model_store <- vector("list", length(risk_weeks_vec))
  names(model_store) <- as.character(risk_weeks_vec)
  for (i in seq_along(weeks_to_fit)) {
    model_store[[as.character(weeks_to_fit[i])]] <- fitted[[i]]
  }
  model_store <- slim_model_store(model_store)

  list(
    model_store = model_store,
    cox_frame_natural = cox_frame_natural
  )
}

gform_model_cache_path <- function(pollutant, dir_models) {
  file.path(dir_models, paste0("natural_course_", pollutant, "_multi.rds"))
}

gform_natural_course_cache_key <- function(
    data_base,
    pollutant,
    risk_weeks_vec,
    control_vars,
    lag_weeks,
    dependent_var,
    risk_entry_week,
    sample_frac = NULL) {

  c(
    list(model_type = "multicontaminant_dlm", cox_frame_slim = TRUE),
    .gform_natural_course_cache_key(
      data_base = data_base,
      pollutant = pollutant,
      risk_weeks_vec = risk_weeks_vec,
      control_vars = control_vars,
      lag_weeks = lag_weeks,
      dependent_var = dependent_var,
      risk_entry_week = risk_entry_week,
      sample_frac = sample_frac
    )
  )
}
