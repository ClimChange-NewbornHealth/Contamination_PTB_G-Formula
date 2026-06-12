# 10.1 G-Formula — funciones modulares (DLM + g-computación con Cox) ----
#
# Etapa 1: generar objetos RDS de intervención (historias semanales intervenidas).
# Etapa 2: ajustar Cox por semana de exposición bajo curso natural (referencia 9.0)
#          y predecir hazards bajo natural vs. contrafactual.
#
# Referencia metodológica: mismo predictor que 9.0 (exposición + lag ponderado +
# TAD semanal + covariables fijas), con Surv(tstart, weeks, birth_preterm).

## Parámetros por defecto ----

GFORM_DEFAULTS <- list(
  dependent_var = "birth_preterm",
  control_vars = c(
    "sex", "age_group_mom", "educ_group_mom", "job_group_mom",
    "age_group_dad", "educ_group_dad", "job_group_dad",
    "month_week1", "year_week1", "covid", "vulnerability", "tad",
    "ndvi_full"
  ),
  krg_contaminants = c("pm25_krg", "o3_krg", "no2_krg"),
  max_follow_up = 37L,
  weeks_exposure = 1:44,
  lag_weeks = 2:44,
  risk_weeks = 1:36,
  risk_entry_week = 28L,
  follow_up_weeks = 28:36,
  baseline_scenario = "observed",
  population_week = 36L,
  boot_iter = 200L,
  boot_seed = 2026L,
  figure4_pct = 0.20
)

## Registro completo de intervenciones globales (ver intervention.txt) ----

GFORM_INTERVENTION_REGISTRY <- list(
  pm25_krg_lt20 = list(
    pollutant = "pm25_krg",
    intervention_id = "pm25_krg_lt20",
    output_stub = "pm25_lt20",
    intervention = list(type = "cap", cap = 20),
    description = "PM2.5 < 20 µg/m³ cada semana gestacional"
  ),
  pm25_krg_lt15 = list(
    pollutant = "pm25_krg",
    intervention_id = "pm25_krg_lt15",
    output_stub = "pm25_lt15",
    intervention = list(type = "cap", cap = 15),
    description = "PM2.5 < 15 µg/m³ cada semana gestacional"
  ),
  pm25_krg_lt10 = list(
    pollutant = "pm25_krg",
    intervention_id = "pm25_krg_lt10",
    output_stub = "pm25_lt10",
    intervention = list(type = "cap", cap = 10),
    description = "PM2.5 < 10 µg/m³ cada semana gestacional"
  ),
  pm25_krg_lt5 = list(
    pollutant = "pm25_krg",
    intervention_id = "pm25_krg_lt5",
    output_stub = "pm25_lt5",
    intervention = list(type = "cap", cap = 5),
    description = "PM2.5 < 5 µg/m³ cada semana gestacional"
  ),
  pm25_krg_pct20 = list(
    pollutant = "pm25_krg",
    intervention_id = "pm25_krg_pct20",
    output_stub = "pm25_pct20",
    intervention = list(type = "pct_reduce", pct = 0.20),
    description = "PM2.5 reducción 20% cada semana gestacional"
  ),
  no2_krg_lt20 = list(
    pollutant = "no2_krg",
    intervention_id = "no2_krg_lt20",
    output_stub = "no2_lt20",
    intervention = list(type = "cap", cap = 20),
    description = "NO2 < 20 ppbv cada semana gestacional"
  ),
  no2_krg_lt15 = list(
    pollutant = "no2_krg",
    intervention_id = "no2_krg_lt15",
    output_stub = "no2_lt15",
    intervention = list(type = "cap", cap = 15),
    description = "NO2 < 15 ppbv cada semana gestacional"
  ),
  no2_krg_lt10 = list(
    pollutant = "no2_krg",
    intervention_id = "no2_krg_lt10",
    output_stub = "no2_lt10",
    intervention = list(type = "cap", cap = 10),
    description = "NO2 < 10 ppbv cada semana gestacional"
  ),
  no2_krg_lt5 = list(
    pollutant = "no2_krg",
    intervention_id = "no2_krg_lt5",
    output_stub = "no2_lt5",
    intervention = list(type = "cap", cap = 5),
    description = "NO2 < 5 ppbv cada semana gestacional"
  ),
  no2_krg_pct20 = list(
    pollutant = "no2_krg",
    intervention_id = "no2_krg_pct20",
    output_stub = "no2_pct20",
    intervention = list(type = "pct_reduce", pct = 0.20),
    description = "NO2 reducción 20% cada semana gestacional"
  ),
  o3_krg_pct20 = list(
    pollutant = "o3_krg",
    intervention_id = "o3_krg_pct20",
    output_stub = "o3_pct20",
    intervention = list(type = "pct_reduce", pct = 0.20),
    description = "O3 reducción 20% cada semana gestacional"
  )
)

`%||%` <- function(x, y) if (is.null(x)) y else x

## Utilidades de formato ancho (alineadas con 9.0) ----

ensure_week_columns <- function(df, prefix, weeks_keep) {
  out <- df
  for (w in weeks_keep) {
    col <- paste0(prefix, "_", w)
    if (!col %in% names(out)) out[[col]] <- NA_real_
  }
  col_order <- c("id", paste0(prefix, "_", weeks_keep))
  out[, col_order[col_order %in% names(out)], drop = FALSE]
}

build_wide_pollutant <- function(df, pollutant, weeks_keep = GFORM_DEFAULTS$weeks_exposure) {
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

build_wide_weekly_var <- function(df, varname, weeks_keep = GFORM_DEFAULTS$weeks_exposure) {
  df |>
    dplyr::select(id, week_gest_num, dplyr::all_of(varname)) |>
    dplyr::filter(week_gest_num %in% weeks_keep) |>
    tidyr::pivot_wider(
      names_from = week_gest_num,
      values_from = dplyr::all_of(varname),
      names_glue = "{.value}_{week_gest_num}"
    )
}

build_wide_raw_exposure <- function(df, pollutant, weeks_keep = GFORM_DEFAULTS$weeks_exposure) {
  out <- df |>
    dplyr::select(id, week_gest_num, dplyr::all_of(pollutant)) |>
    dplyr::filter(week_gest_num %in% weeks_keep) |>
    tidyr::pivot_wider(
      names_from = week_gest_num,
      values_from = dplyr::all_of(pollutant),
      names_prefix = paste0(pollutant, "_")
    )
  ensure_week_columns(out, pollutant, weeks_keep)
}

weighted_lag_from_week_vector <- function(expo, w, week_min = 1L) {
  if (w <= week_min || w > length(expo)) return(NA_real_)
  past <- week_min:(w - 1L)
  sum(expo[past] / (w - past), na.rm = TRUE)
}

apply_exposure_intervention_vec <- function(
    expo,
    intervention,
    single_week = NULL,
    weeks_keep = NULL) {

  if (intervention$type == "none" && is.null(single_week)) return(expo)

  weeks_idx <- if (is.null(weeks_keep)) seq_along(expo) else weeks_keep
  target <- if (is.null(single_week)) {
    rep(TRUE, length(expo))
  } else {
    weeks_idx == single_week
  }

  out <- expo
  if (intervention$type == "cap") {
    out[target] <- pmin(expo[target], intervention$cap, na.rm = FALSE)
  } else if (intervention$type == "pct_reduce") {
    out[target] <- expo[target] * (1 - intervention$pct)
  } else if (intervention$type != "none") {
    stop("Intervención desconocida: ", intervention$type)
  }
  out
}

build_exposicion_wide_from_raw <- function(
    raw_wide,
    pollutant,
    intervention,
    weeks_keep = GFORM_DEFAULTS$weeks_exposure,
    lag_weeks = GFORM_DEFAULTS$lag_weeks,
    single_week = NULL) {

  week_cols <- paste0(pollutant, "_", weeks_keep)
  mat <- as.matrix(raw_wide[, week_cols, drop = FALSE])
  rownames(mat) <- raw_wide$id

  mat_cf <- if (intervention$type == "none" && is.null(single_week)) {
    mat
  } else {
    t(apply(mat, 1L, function(row) {
      apply_exposure_intervention_vec(
        row, intervention,
        single_week = single_week,
        weeks_keep = weeks_keep
      )
    }))
  }

  lag_mat <- matrix(NA_real_, nrow = nrow(mat_cf), ncol = ncol(mat_cf),
                    dimnames = dimnames(mat_cf))
  for (j in seq_len(ncol(mat_cf))) {
    w <- weeks_keep[j]
    if (w %in% lag_weeks) {
      lag_mat[, j] <- vapply(seq_len(nrow(mat_cf)), function(i) {
        weighted_lag_from_week_vector(mat_cf[i, ], w, week_min = min(weeks_keep))
      }, numeric(1))
    }
  }

  expo_df <- data.table::as.data.table(mat_cf)
  data.table::setnames(expo_df, paste0("exposicion_", weeks_keep))
  expo_df[, id := raw_wide$id]

  lag_df <- data.table::as.data.table(lag_mat)
  data.table::setnames(lag_df, paste0("exposicion_lagged_", weeks_keep))
  lag_df[, id := raw_wide$id]

  merge(expo_df, lag_df, by = "id", sort = FALSE)
}

expand_person_weeks <- function(
    births_df,
    risk_weeks,
    max_fu,
    dependent_var = GFORM_DEFAULTS$dependent_var) {

  dt <- data.table::as.data.table(births_df)
  dt[, bw_floor := pmin(floor(pmin(weeks, max_fu)), max_fu)]
  dt[, last_rw := pmin(max(risk_weeks), bw_floor)]

  out <- dt[last_rw >= min(risk_weeks), {
    ws <- risk_weeks[risk_weeks <= last_rw]
    data.table::data.table(
      id = id,
      time = ws,
      event = as.integer(get(dependent_var) == 1L & ws == bw_floor)
    )
  }, by = .(id, bw_floor, last_rw)]

  out[, c("bw_floor", "last_rw") := NULL]
  as.data.frame(out)
}

## Etapa 1: objetos de intervención ----

generate_intervention <- function(
    raw_wide_pollutant,
    pollutant,
    intervention,
    intervention_id,
    description = NULL,
    output_path,
    weeks_keep = GFORM_DEFAULTS$weeks_exposure,
    lag_weeks = GFORM_DEFAULTS$lag_weeks,
    overwrite = FALSE) {

  if (file.exists(output_path) && !overwrite) {
    message("Intervención existente (omitida): ", output_path)
    return(invisible(readRDS(output_path)))
  }

  wide_exposicion <- build_exposicion_wide_from_raw(
    raw_wide = raw_wide_pollutant,
    pollutant = pollutant,
    intervention = intervention,
    weeks_keep = weeks_keep,
    lag_weeks = lag_weeks
  )

  intervention_obj <- list(
    intervention_id = intervention_id,
    pollutant = pollutant,
    intervention = intervention,
    description = description %||% intervention_id,
    weeks_exposure = weeks_keep,
    lag_weeks = lag_weeks,
    wide_exposicion = wide_exposicion,
    created = Sys.time()
  )

  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(intervention_obj, output_path)
  invisible(intervention_obj)
}

build_all_interventions <- function(
    data_long,
    registry = GFORM_INTERVENTION_REGISTRY,
    output_dir,
    weeks_keep = GFORM_DEFAULTS$weeks_exposure,
    lag_weeks = GFORM_DEFAULTS$lag_weeks,
    overwrite = FALSE,
    parallel = TRUE,
    n_workers = NULL) {

  pollutants_needed <- unique(vapply(registry, function(x) x$pollutant, character(1)))
  raw_by_pollutant <- lapply(
    pollutants_needed,
    function(p) build_wide_raw_exposure(data_long, p, weeks_keep = weeks_keep)
  )
  names(raw_by_pollutant) <- pollutants_needed

  build_one <- function(key) {
    spec <- registry[[key]]
    raw_wide <- raw_by_pollutant[[spec$pollutant]]
    output_path <- file.path(output_dir, paste0(spec$intervention_id, ".rds"))
    generate_intervention(
      raw_wide_pollutant = raw_wide,
      pollutant = spec$pollutant,
      intervention = spec$intervention,
      intervention_id = spec$intervention_id,
      description = spec$description,
      output_path = output_path,
      weeks_keep = weeks_keep,
      lag_weeks = lag_weeks,
      overwrite = overwrite
    )
    output_path
  }

  registry_keys <- names(registry)
  if (parallel && length(registry_keys) > 1L && requireNamespace("furrr", quietly = TRUE)) {
    furrr::future_map(registry_keys, build_one, .options = furrr::furrr_options(seed = TRUE))
  } else {
    lapply(registry_keys, build_one)
  }
}

## Cox: utilidades de hazard discreto ----

get_bh_increment <- function(bh, time_point) {
  if (is.null(bh) || nrow(bh) == 0L || !is.finite(time_point)) return(0)
  h_at <- bh$hazard[bh$time <= time_point]
  h_prev <- bh$hazard[bh$time <= (time_point - 1)]
  val_at <- if (length(h_at)) max(h_at, na.rm = TRUE) else 0
  val_prev <- if (length(h_prev)) max(h_prev, na.rm = TRUE) else 0
  if (!is.finite(val_at)) val_at <- 0
  if (!is.finite(val_prev)) val_prev <- 0
  max(val_at - val_prev, 0)
}

lp_to_p_event <- function(lp, bh, time_point) {
  dh0 <- get_bh_increment(bh, time_point)
  if (dh0 <= 0) return(rep(0, length(lp)))
  p_event <- 1 - exp(-dh0 * exp(lp))
  pmin(pmax(p_event, 0), 1)
}

build_cox_model_frame <- function(
    data_base,
    wide_exposicion,
    wide_tad_obs,
    control_vars,
    dependent_var,
    risk_entry_week) {

  base_controls <- data.table::as.data.table(
    data_base[, c("id", "weeks", dependent_var, control_vars), drop = FALSE]
  )
  base_controls[, tstart := risk_entry_week]

  if (!data.table::is.data.table(wide_exposicion)) {
    wide_exposicion <- data.table::as.data.table(wide_exposicion)
  }
  if (!data.table::is.data.table(wide_tad_obs)) {
    wide_tad_obs <- data.table::as.data.table(wide_tad_obs)
  }

  data.table::setkey(base_controls, id)
  data.table::setkey(wide_exposicion, id)
  data.table::setkey(wide_tad_obs, id)

  out <- wide_exposicion[base_controls, on = "id"]
  out <- wide_tad_obs[out, on = "id"]
  out
}

predict_cox_lp <- function(cox_fit, newdata, coef_override = NULL) {
  if (is.null(coef_override)) {
    return(as.numeric(stats::predict(cox_fit, newdata = newdata, type = "lp")))
  }
  mm <- stats::model.matrix(stats::delete.response(stats::terms(cox_fit)), newdata)
  beta <- coef_override
  v <- rep(0, ncol(mm))
  names(v) <- colnames(mm)
  ii <- intersect(names(v), names(beta))
  v[ii] <- unname(beta[ii])
  fit_b <- stats::coef(cox_fit)
  bad <- !is.finite(v)
  if (any(bad)) {
    jj <- intersect(names(v), names(fit_b))
    v[bad] <- fit_b[jj[match(names(v)[bad], jj)]]
    v[!is.finite(v)] <- 0
  }
  drop(mm %*% v)
}

## Etapa 2: ajuste Cox (solo curso natural; referencia 9.0) ----

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

  data_model <- build_cox_model_frame(
    data_base = data_base,
    wide_exposicion = wide_exposicion_natural,
    wide_tad_obs = wide_tad_obs,
    control_vars = control_vars,
    dependent_var = dependent_var,
    risk_entry_week = risk_entry_week
  )
  data_model_df <- as.data.frame(data_model)

  fit_one_week <- function(rw) {
    if (rw < risk_entry_week) return(NULL)

    exp_var <- paste0("exposicion_", rw)
    lag_var <- paste0("exposicion_lagged_", rw)
    tad_var <- paste0("tad_", rw)
    pred_terms <- c(
      exp_var,
      lag_var[lag_var %in% paste0("exposicion_lagged_", lag_weeks)],
      tad_var
    )

    rhs <- paste(c(pred_terms, control_vars), collapse = " + ")
    fml <- stats::as.formula(paste0(
      "Surv(tstart, weeks, ", dependent_var, ") ~ ", rhs
    ))

    vars_needed <- c("id", "weeks", "tstart", dependent_var, pred_terms, control_vars)
    model_df <- data_model_df[, vars_needed, drop = FALSE]
    model_df <- stats::na.omit(model_df)

    if (nrow(model_df) < 50L) {
      warning("Muy pocas filas para la semana ", rw, "; modelo Cox omitido.")
      return(NULL)
    }

    cox_fit <- tryCatch(
      survival::coxph(fml, data = model_df, x = TRUE, y = TRUE),
      error = function(e) NULL
    )
    if (is.null(cox_fit)) {
      warning("Cox no convergió o error en semana ", rw)
      return(NULL)
    }

    list(
      model = cox_fit,
      formula = fml,
      id_order = model_df$id,
      pred_terms = pred_terms,
      control_vars = control_vars,
      baseline_hazard = survival::basehaz(cox_fit, centered = FALSE),
      risk_week = rw
    )
  }

  weeks_to_fit <- risk_weeks_vec[risk_weeks_vec >= risk_entry_week]
  if (parallel && length(weeks_to_fit) > 1L && requireNamespace("furrr", quietly = TRUE)) {
    fitted <- furrr::future_map(weeks_to_fit, fit_one_week, .options = furrr::furrr_options(seed = TRUE))
  } else {
    fitted <- lapply(weeks_to_fit, fit_one_week)
  }

  model_store <- vector("list", length(risk_weeks_vec))
  names(model_store) <- as.character(risk_weeks_vec)
  for (i in seq_along(weeks_to_fit)) {
    model_store[[as.character(weeks_to_fit[i])]] <- fitted[[i]]
  }
  model_store
}

## Predicción de hazards semanales (Cox → probabilidad discreta) ----

predict_weekly_hazards <- function(
    model_store,
    person_weeks,
    data_base,
    wide_exposicion,
    wide_tad_obs,
    risk_weeks_vec,
    control_vars = GFORM_DEFAULTS$control_vars,
    dependent_var = GFORM_DEFAULTS$dependent_var,
    risk_entry_week = GFORM_DEFAULTS$risk_entry_week,
    coef_override = NULL) {

  cox_frame <- build_cox_model_frame(
    data_base = data_base,
    wide_exposicion = wide_exposicion,
    wide_tad_obs = wide_tad_obs,
    control_vars = control_vars,
    dependent_var = dependent_var,
    risk_entry_week = risk_entry_week
  )

  out_list <- vector("list", length(risk_weeks_vec))
  names(out_list) <- as.character(risk_weeks_vec)

  for (rw in risk_weeks_vec) {
    pw_rw <- person_weeks[person_weeks$time == rw, , drop = FALSE]
    if (nrow(pw_rw) == 0L) next

    if (rw < risk_entry_week) {
      out_list[[as.character(rw)]] <- data.table::data.table(
        id = pw_rw$id,
        time = rw,
        p_noevent = 1
      )
      next
    }

    ms <- model_store[[as.character(rw)]]
    if (is.null(ms)) next

    risk_dt <- data.table::as.data.table(pw_rw[, "id", drop = FALSE])
    data.table::setkey(risk_dt, id)
    data.table::setkey(cox_frame, id)
    newdata <- cox_frame[risk_dt, on = "id"]
    newdata <- newdata[id %in% ms$id_order]
    newdata <- newdata[match(ms$id_order, newdata$id)]

    newdata_df <- as.data.frame(newdata)
    beta_rw <- if (is.null(coef_override)) NULL else coef_override[[as.character(rw)]]

    lp <- predict_cox_lp(ms$model, newdata_df, coef_override = beta_rw)
    p_event <- lp_to_p_event(lp, ms$baseline_hazard, rw)

    out_list[[as.character(rw)]] <- data.table::data.table(
      id = newdata$id,
      time = rw,
      p_noevent = 1 - p_event
    )
  }

  result <- data.table::rbindlist(out_list, use.names = TRUE, fill = TRUE)
  if (is.null(result) || nrow(result) == 0L) {
    stop("predict_weekly_hazards: no se generaron predicciones válidas.")
  }
  result
}

## Supervivencia y efectos ----

compute_survival <- function(prob_dt) {
  dt <- data.table::copy(data.table::as.data.table(prob_dt))
  data.table::setorder(dt, id, time)
  dt[, surv := cumprod(p_noevent), by = id]
  dt[, risk := 1 - surv]
  dt
}

compute_weekly_effects <- function(prob_natural, prob_intervention) {
  surv_nat <- compute_survival(prob_natural)
  surv_int <- compute_survival(prob_intervention)

  mean_nat <- surv_nat[, .(risk_natural = mean(risk, na.rm = TRUE)), by = time]
  mean_int <- surv_int[, .(risk_intervention = mean(risk, na.rm = TRUE)), by = time]

  out <- merge(mean_nat, mean_int, by = "time", all = TRUE, sort = TRUE)
  out[, `:=`(
    risk_ratio = data.table::fifelse(
      risk_natural > 0, risk_intervention / risk_natural, NA_real_
    ),
    risk_difference = risk_intervention - risk_natural
  )]
  data.table::setnames(out, "time", "week")
  tibble::as_tibble(out)
}

compute_figure3_curves <- function(
    weekly_effects,
    follow_up_weeks = GFORM_DEFAULTS$follow_up_weeks) {

  dt <- data.table::as.data.table(weekly_effects)
  dt <- dt[week %in% follow_up_weeks]
  dt[, `:=`(
    follow_up_week = week,
    cumulative_risk_observed = risk_natural,
    cumulative_risk_intervention = risk_intervention,
    risk_difference = risk_difference
  )]
  tibble::as_tibble(dt[, .(
    follow_up_week,
    cumulative_risk_observed,
    cumulative_risk_intervention,
    risk_ratio,
    risk_difference
  )])
}

compute_population_effects <- function(
    prob_natural,
    prob_intervention,
    total_births,
    target_week = GFORM_DEFAULTS$population_week,
    baseline_scenario = "observed",
    intervention_scenario = "intervention") {

  build_scenario_probs <- function(prob_dt, scenario_name) {
    dt <- compute_survival(prob_dt)
    dt <- dt[time <= target_week]
    dt <- dt[, .SD[which.max(time)], by = id]
    dt[, scenario := scenario_name]
    dt[, .(id, time, p_noevent, scenario, risk)]
  }

  nat_tail <- build_scenario_probs(prob_natural, baseline_scenario)
  int_tail <- build_scenario_probs(prob_intervention, intervention_scenario)
  combined <- data.table::rbindlist(list(nat_tail, int_tail), use.names = TRUE)

  risk_df <- combined[, .(prevalence = mean(risk, na.rm = TRUE)), by = scenario]
  risk_df[, cases := prevalence * total_births]

  baseline_risk <- risk_df[scenario == baseline_scenario, prevalence]

  risk_df[, `:=`(
    risk_ratio = prevalence / baseline_risk,
    risk_difference = prevalence - baseline_risk,
    attributable_risk = baseline_risk - prevalence
  )]

  tibble::as_tibble(risk_df[, .(
    scenario, prevalence, cases, risk_ratio, risk_difference, attributable_risk
  )])
}

compute_figure4_heatmap <- function(
    model_store,
    person_weeks,
    data_base,
    raw_wide_pollutant,
    pollutant,
    wide_tad_obs,
    wide_exposicion_natural,
    risk_weeks_vec,
    control_vars,
    intervention_weeks = GFORM_DEFAULTS$weeks_exposure,
    follow_up_weeks = GFORM_DEFAULTS$follow_up_weeks,
    pct = GFORM_DEFAULTS$figure4_pct,
    dependent_var = GFORM_DEFAULTS$dependent_var,
    risk_entry_week = GFORM_DEFAULTS$risk_entry_week,
    parallel = TRUE) {

  single_intervention <- list(type = "pct_reduce", pct = pct)

  prob_natural <- predict_weekly_hazards(
    model_store = model_store,
    person_weeks = person_weeks,
    data_base = data_base,
    wide_exposicion = wide_exposicion_natural,
    wide_tad_obs = wide_tad_obs,
    risk_weeks_vec = risk_weeks_vec,
    control_vars = control_vars,
    dependent_var = dependent_var,
    risk_entry_week = risk_entry_week
  )
  surv_nat <- compute_survival(prob_natural)
  nat_mean <- surv_nat[, .(risk_natural = mean(risk, na.rm = TRUE)), by = time]

  compute_one_col <- function(iw) {
    wide_cf <- build_exposicion_wide_from_raw(
      raw_wide = raw_wide_pollutant,
      pollutant = pollutant,
      intervention = single_intervention,
      single_week = iw
    )
    prob_int <- predict_weekly_hazards(
      model_store = model_store,
      person_weeks = person_weeks,
      data_base = data_base,
      wide_exposicion = wide_cf,
      wide_tad_obs = wide_tad_obs,
      risk_weeks_vec = risk_weeks_vec,
      control_vars = control_vars,
      dependent_var = dependent_var,
      risk_entry_week = risk_entry_week
    )
    surv_int <- compute_survival(prob_int)
    int_mean <- surv_int[, .(risk_intervention = mean(risk, na.rm = TRUE)), by = time]
    merged <- merge(nat_mean, int_mean, by = "time", all = TRUE, sort = TRUE)
    merged[, `:=`(
      intervention_week = iw,
      follow_up_week = time,
      risk_difference = risk_intervention - risk_natural
    )]
    merged[follow_up_week %in% follow_up_weeks,
             .(intervention_week, follow_up_week, risk_difference, risk_natural, risk_intervention)]
  }

  if (parallel && length(intervention_weeks) > 1L && requireNamespace("furrr", quietly = TRUE)) {
    pieces <- furrr::future_map(
      intervention_weeks,
      compute_one_col,
      .options = furrr::furrr_options(seed = TRUE)
    )
  } else {
    pieces <- lapply(intervention_weeks, compute_one_col)
  }

  long_df <- data.table::rbindlist(pieces, use.names = TRUE, fill = TRUE)
  wide_mat <- data.table::dcast(
    long_df,
    follow_up_week ~ intervention_week,
    value.var = "risk_difference"
  )
  list(
    long = tibble::as_tibble(long_df),
    wide = tibble::as_tibble(wide_mat)
  )
}

## Bootstrap paramétrico (coeficientes Cox) ----

simulate_coefs_list <- function(model_store, risk_weeks_vec) {
  coefs_list <- vector("list", length(risk_weeks_vec))
  names(coefs_list) <- as.character(risk_weeks_vec)
  for (rw in risk_weeks_vec) {
    ms <- model_store[[as.character(rw)]]
    if (is.null(ms)) next
    b <- stats::coef(ms$model)
    V <- stats::vcov(ms$model)
    if (any(!is.finite(V))) {
      sim_b <- b
    } else {
      V_stable <- tryCatch(
        as.matrix(Matrix::nearPD(V, corr = FALSE)$mat),
        error = function(e) V + diag(1e-6, nrow(V))
      )
      sim_b <- tryCatch(
        MASS::mvrnorm(1L, mu = b, Sigma = V_stable),
        error = function(e) b
      )
    }
    sim_b <- as.numeric(sim_b)
    names(sim_b) <- names(b)
    bad <- !is.finite(sim_b)
    if (any(bad)) sim_b[bad] <- b[bad]
    coefs_list[[as.character(rw)]] <- sim_b
  }
  coefs_list
}

bootstrap_gformula_effects <- function(
    model_store,
    person_weeks,
    data_base,
    wide_exposicion_natural,
    wide_exposicion_intervention,
    wide_tad_obs,
    risk_weeks_vec,
    control_vars,
    total_births,
    dependent_var = GFORM_DEFAULTS$dependent_var,
    risk_entry_week = GFORM_DEFAULTS$risk_entry_week,
    boot_iter = GFORM_DEFAULTS$boot_iter,
    boot_seed = GFORM_DEFAULTS$boot_seed,
    target_week = GFORM_DEFAULTS$population_week,
    baseline_scenario = GFORM_DEFAULTS$baseline_scenario,
    intervention_scenario = "intervention",
    parallel = TRUE) {

  boot_fn <- function(iter) {
    coefs_list <- simulate_coefs_list(model_store, risk_weeks_vec)
    prob_nat <- predict_weekly_hazards(
      model_store, person_weeks, data_base,
      wide_exposicion_natural, wide_tad_obs, risk_weeks_vec,
      control_vars = control_vars,
      dependent_var = dependent_var,
      risk_entry_week = risk_entry_week,
      coef_override = coefs_list
    )
    prob_int <- predict_weekly_hazards(
      model_store, person_weeks, data_base,
      wide_exposicion_intervention, wide_tad_obs, risk_weeks_vec,
      control_vars = control_vars,
      dependent_var = dependent_var,
      risk_entry_week = risk_entry_week,
      coef_override = coefs_list
    )
    list(
      iter = iter,
      weekly = compute_weekly_effects(prob_nat, prob_int),
      population = compute_population_effects(
        prob_nat, prob_int, total_births,
        target_week = target_week,
        baseline_scenario = baseline_scenario,
        intervention_scenario = intervention_scenario
      )
    )
  }

  if (parallel && requireNamespace("furrr", quietly = TRUE)) {
    boot_list <- furrr::future_map(
      seq_len(boot_iter),
      boot_fn,
      .options = furrr::furrr_options(seed = TRUE)
    )
  } else {
    set.seed(boot_seed)
    boot_list <- lapply(seq_len(boot_iter), boot_fn)
  }

  weekly_boot <- data.table::rbindlist(
    lapply(boot_list, function(x) {
      dt <- data.table::as.data.table(x$weekly)
      dt[, iter := x$iter]
      dt
    }),
    use.names = TRUE, fill = TRUE
  )

  pop_boot <- data.table::rbindlist(
    lapply(boot_list, function(x) {
      dt <- data.table::as.data.table(x$population)
      dt[, iter := x$iter]
      dt
    }),
    use.names = TRUE, fill = TRUE
  )

  weekly_ci <- weekly_boot[, .(
    risk_natural_lcl = stats::quantile(risk_natural, 0.025, na.rm = TRUE),
    risk_natural_ucl = stats::quantile(risk_natural, 0.975, na.rm = TRUE),
    risk_intervention_lcl = stats::quantile(risk_intervention, 0.025, na.rm = TRUE),
    risk_intervention_ucl = stats::quantile(risk_intervention, 0.975, na.rm = TRUE),
    risk_ratio_lcl = stats::quantile(risk_ratio, 0.025, na.rm = TRUE),
    risk_ratio_ucl = stats::quantile(risk_ratio, 0.975, na.rm = TRUE),
    risk_difference_lcl = stats::quantile(risk_difference, 0.025, na.rm = TRUE),
    risk_difference_ucl = stats::quantile(risk_difference, 0.975, na.rm = TRUE)
  ), by = week]

  pop_ci <- pop_boot[, .(
    prevalence_lcl = stats::quantile(prevalence, 0.025, na.rm = TRUE),
    prevalence_ucl = stats::quantile(prevalence, 0.975, na.rm = TRUE),
    risk_ratio_lcl = stats::quantile(risk_ratio, 0.025, na.rm = TRUE),
    risk_ratio_ucl = stats::quantile(risk_ratio, 0.975, na.rm = TRUE),
    risk_difference_lcl = stats::quantile(risk_difference, 0.025, na.rm = TRUE),
    risk_difference_ucl = stats::quantile(risk_difference, 0.975, na.rm = TRUE),
    cases_lcl = stats::quantile(cases, 0.025, na.rm = TRUE),
    cases_ucl = stats::quantile(cases, 0.975, na.rm = TRUE),
    attributable_risk_lcl = stats::quantile(attributable_risk, 0.025, na.rm = TRUE),
    attributable_risk_ucl = stats::quantile(attributable_risk, 0.975, na.rm = TRUE)
  ), by = scenario]

  list(
    weekly_boot = weekly_boot,
    population_boot = pop_boot,
    weekly_ci = tibble::as_tibble(weekly_ci),
    population_ci = tibble::as_tibble(pop_ci)
  )
}

## Orquestador por intervención ----

run_gform_intervention <- function(
    intervention_spec,
    intervention_path,
    data_base,
    data_long,
    wide_tad_obs,
    model_store,
    person_weeks,
    risk_weeks_vec,
    control_vars,
    total_births,
    dependent_var = GFORM_DEFAULTS$dependent_var,
    risk_entry_week = GFORM_DEFAULTS$risk_entry_week,
    follow_up_weeks = GFORM_DEFAULTS$follow_up_weeks,
    boot_iter = GFORM_DEFAULTS$boot_iter,
    boot_seed = GFORM_DEFAULTS$boot_seed,
    target_week = GFORM_DEFAULTS$population_week,
    run_bootstrap = TRUE,
    run_figure4 = FALSE,
    parallel = TRUE) {

  intervention_obj <- readRDS(intervention_path)
  pollutant <- intervention_spec$pollutant

  raw_pm <- build_wide_raw_exposure(
    data_long, pollutant, weeks_keep = GFORM_DEFAULTS$weeks_exposure
  )
  wide_exposicion_natural <- build_exposicion_wide_from_raw(
    raw_wide = raw_pm,
    pollutant = pollutant,
    intervention = list(type = "none")
  )
  wide_exposicion_intervention <- intervention_obj$wide_exposicion
  wide_exposicion_intervention <- wide_exposicion_intervention[id %in% data_base$id]

  prob_natural <- predict_weekly_hazards(
    model_store = model_store,
    person_weeks = person_weeks,
    data_base = data_base,
    wide_exposicion = wide_exposicion_natural,
    wide_tad_obs = wide_tad_obs,
    risk_weeks_vec = risk_weeks_vec,
    control_vars = control_vars,
    dependent_var = dependent_var,
    risk_entry_week = risk_entry_week
  )

  prob_intervention <- predict_weekly_hazards(
    model_store = model_store,
    person_weeks = person_weeks,
    data_base = data_base,
    wide_exposicion = wide_exposicion_intervention,
    wide_tad_obs = wide_tad_obs,
    risk_weeks_vec = risk_weeks_vec,
    control_vars = control_vars,
    dependent_var = dependent_var,
    risk_entry_week = risk_entry_week
  )

  weekly_effects <- compute_weekly_effects(prob_natural, prob_intervention)
  figure3 <- compute_figure3_curves(weekly_effects, follow_up_weeks = follow_up_weeks)
  population_effects <- compute_population_effects(
    prob_natural = prob_natural,
    prob_intervention = prob_intervention,
    total_births = total_births,
    target_week = target_week
  )

  boot_out <- NULL
  weekly_ci <- weekly_effects
  population_ci <- population_effects
  if (run_bootstrap && boot_iter > 0L) {
    boot_out <- bootstrap_gformula_effects(
      model_store = model_store,
      person_weeks = person_weeks,
      data_base = data_base,
      wide_exposicion_natural = wide_exposicion_natural,
      wide_exposicion_intervention = wide_exposicion_intervention,
      wide_tad_obs = wide_tad_obs,
      risk_weeks_vec = risk_weeks_vec,
      control_vars = control_vars,
      total_births = total_births,
      dependent_var = dependent_var,
      risk_entry_week = risk_entry_week,
      boot_iter = boot_iter,
      boot_seed = boot_seed,
      target_week = target_week,
      parallel = parallel
    )
    weekly_ci <- dplyr::left_join(weekly_effects, boot_out$weekly_ci, by = "week")
    population_ci <- dplyr::left_join(population_effects, boot_out$population_ci, by = "scenario")
  }

  figure4 <- NULL
  if (run_figure4) {
    figure4 <- compute_figure4_heatmap(
      model_store = model_store,
      person_weeks = person_weeks,
      data_base = data_base,
      raw_wide_pollutant = raw_pm,
      pollutant = pollutant,
      wide_tad_obs = wide_tad_obs,
      wide_exposicion_natural = wide_exposicion_natural,
      risk_weeks_vec = risk_weeks_vec,
      control_vars = control_vars,
      parallel = parallel
    )
  }

  list(
    intervention_spec = intervention_spec,
    weekly_effects = weekly_ci,
    population_effects = population_ci,
    figure3 = figure3,
    figure4 = figure4,
    bootstrap = boot_out
  )
}

## Guardado de resultados ----

save_results <- function(
    weekly_effects,
    population_effects,
    weekly_path,
    population_path,
    figure3 = NULL,
    figure3_path = NULL,
    figure4 = NULL,
    figure4_path = NULL,
    weekly_boot = NULL,
    population_boot = NULL,
    metadata = list()) {

  weekly_obj <- c(
    list(point = weekly_effects, bootstrap = weekly_boot),
    metadata
  )
  population_obj <- c(
    list(point = population_effects, bootstrap = population_boot),
    metadata
  )

  dir.create(dirname(weekly_path), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(population_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(weekly_obj, weekly_path)
  saveRDS(population_obj, population_path)

  if (!is.null(figure3) && !is.null(figure3_path)) {
    dir.create(dirname(figure3_path), recursive = TRUE, showWarnings = FALSE)
    saveRDS(c(list(point = figure3), metadata), figure3_path)
  }
  if (!is.null(figure4) && !is.null(figure4_path)) {
    dir.create(dirname(figure4_path), recursive = TRUE, showWarnings = FALSE)
    saveRDS(c(list(point = figure4), metadata), figure4_path)
  }

  invisible(list(weekly = weekly_path, population = population_path))
}

save_gform_excel <- function(
    results,
    excel_path,
    intervention_id) {

  sheets <- list(
    weekly_effects = results$weekly_effects,
    population_effects = results$population_effects,
    figure3 = results$figure3
  )
  if (!is.null(results$figure4)) {
    sheets$figure4_long <- results$figure4$long
    sheets$figure4_wide <- results$figure4$wide
  }

  dir.create(dirname(excel_path), recursive = TRUE, showWarnings = FALSE)
  writexl::write_xlsx(sheets, path = excel_path)
  invisible(excel_path)
}

gform_setup_parallel <- function(n_workers = NULL) {
  if (is.null(n_workers)) {
    n_detect <- parallel::detectCores(logical = TRUE)
    if (!is.finite(n_detect) || n_detect < 1L) n_detect <- 4L
    n_workers <- max(1L, n_detect - 4L)
  }
  future::plan(future::multisession, workers = n_workers)
  options(future.globals.maxSize = 16 * 1024^3)
  invisible(n_workers)
}
