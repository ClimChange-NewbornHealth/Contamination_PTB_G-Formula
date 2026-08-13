# 12.0 G-Formula por períodos — funciones (trimestre conjunto + período completo) ----
#
# Evalúa intervenciones sobre exposiciones agregadas por trimestre (T1+T2+T3
# en un mismo modelo Cox ajustado, como 11.x) o sobre el período completo
# (variable *_full + covariables). Reutiliza bootstrap y g-computación de 10.1.
#
# Alcance:
#   - application_scope = "trimester": misma regla en pm25_krg_t1, t2 y t3
#   - application_scope = "full": regla sobre pm25_krg_full
# Solo modelos ajustados (covariables fijas al nivel del nacimiento).

source("00_Code/10.1 G-Form_functions.R")

## Parámetros período ----

GFORM_PERIOD_TRIMESTER_WEEKS <- list(
  t1 = 1:13,
  t2 = 14:26,
  t3 = 27:44
)

GFORM_PERIOD_CONTROL_VARS <- c(
  "sex", "age_group_mom", "educ_group_mom", "job_group_mom",
  "age_group_dad", "educ_group_dad", "job_group_dad",
  "month_week1", "year_week1", "covid", "vulnerability",
  "tad_full", "ndvi_full"
)

GFORM_PERIOD_DEFAULTS <- utils::modifyList(GFORM_DEFAULTS, list(
  control_vars = GFORM_PERIOD_CONTROL_VARS,
  risk_weeks = 28:36,
  output_cumulative_risk_suffix = "curvas_riesgo_acumulado_periodo",
  run_singleweek_heatmap = FALSE,
  minimal_output = TRUE,
  boot_iter = 0L
))

## Registro de intervenciones (34 = 17 trimestre + 17 full) ----

gform_period_pollutant_short <- function(pollutant) {
  sub("_krg$", "", pollutant)
}

gform_period_exposure_var_names <- function(pollutant, application_scope) {
  short <- gform_period_pollutant_short(pollutant)
  if (identical(application_scope, "full")) {
    paste0(short, "_krg_full")
  } else if (identical(application_scope, "trimester")) {
    c(
      paste0(short, "_krg_t1"),
      paste0(short, "_krg_t2"),
      paste0(short, "_krg_t3")
    )
  } else {
    stop("application_scope desconocido: ", application_scope)
  }
}

gform_period_pollutant_label <- function(pollutant) {
  switch(
    pollutant,
    pm25_krg = "PM2.5",
    no2_krg = "NO2",
    o3_krg = "O3",
    pollutant
  )
}

gform_period_intervention_rule_label <- function(pollutant, intervention) {
  pol <- gform_period_pollutant_label(pollutant)
  if (identical(intervention$type, "cap")) {
    unit <- switch(
      pollutant,
      no2_krg = " ppbv",
      o3_krg = " ppbv",
      " µg/m³"
    )
    paste0(pol, " cap < ", intervention$cap, unit)
  } else if (identical(intervention$type, "pct_reduce")) {
    paste0(pol, " reducción ", round(100 * intervention$pct), "%")
  } else {
    pol
  }
}

gform_period_intervention_description <- function(
    pollutant,
    intervention,
    application_scope) {

  rule <- gform_period_intervention_rule_label(pollutant, intervention)
  if (identical(application_scope, "trimester")) {
    paste0(
      rule,
      " — misma regla en medias T1, T2 y T3; ",
      "Cox conjunto (T1+T2+T3 mutuamente ajustados + covariables)"
    )
  } else if (identical(application_scope, "full")) {
    paste0(
      rule,
      " — media del período gestacional completo; ",
      "Cox (exposición full + covariables)"
    )
  } else {
    stop("application_scope desconocido: ", application_scope)
  }
}

gform_period_cox_formula_text <- function(pollutant, application_scope) {
  dep <- GFORM_PERIOD_DEFAULTS$dependent_var
  cov <- paste(GFORM_PERIOD_CONTROL_VARS, collapse = " + ")
  short <- gform_period_pollutant_short(pollutant)
  if (identical(application_scope, "trimester")) {
    expo <- paste(
      paste0(short, "_krg_t1"),
      paste0(short, "_krg_t2"),
      paste0(short, "_krg_t3"),
      sep = " + "
    )
    paste0(
      "Surv(tstart, weeks, ", dep, ") ~ ", expo, " + ", cov,
      "  [T1|T2|T3 mutuamente ajustados]"
    )
  } else {
    expo <- paste0(short, "_krg_full")
    paste0("Surv(tstart, weeks, ", dep, ") ~ ", expo, " + ", cov)
  }
}

print_gform_period_cox_specs <- function() {
  message("\n=== Especificaciones Cox (sección 12) ===")
  message("Por contaminante se estiman 2 curvas naturales (cache en Models/):")
  for (pol in GFORM_DEFAULTS$krg_contaminants) {
    message("  ", gform_period_pollutant_label(pol), ":")
    message("    • Trimestre (_tri): ", gform_period_cox_formula_text(pol, "trimester"))
    message("    • Período completo (_full): ", gform_period_cox_formula_text(pol, "full"))
  }
  message(
    "El modelo trimestral conjunto equivale a T1 ajustado por T2+T3, ",
    "T2 por T1+T3 y T3 por T1+T2 (un solo ajuste, como 11.x t1_t2_t3)."
  )
  message("34 escenarios = 17 reglas × 2 alcances (orden: % tri → % full → caps tri → caps full).")
  message(
    "Salida: efectos poblacionales (prevalencia, riesgo, RD, RR, PAF); ",
    "sin curvas semanales ni heatmap."
  )
  message(
    "Bootstrap (IC): opcional — GFORM_RUN_BOOTSTRAP=true GFORM_BOOT_ITER=250"
  )
  invisible(NULL)
}

build_gform_period_intervention_registry <- function() {
  reg <- list()
  keys <- gform_intervention_keys()
  for (key in keys) {
    base <- GFORM_INTERVENTION_REGISTRY[[key]]
    tri_id <- paste0(base$intervention_id, "_tri")
    reg[[tri_id]] <- modifyList(base, list(
      intervention_id = tri_id,
      output_stub = paste0(base$output_stub, "_tri"),
      application_scope = "trimester",
      exposure_model = "aggregate_trimester_joint",
      cox_formula = gform_period_cox_formula_text(base$pollutant, "trimester"),
      description = gform_period_intervention_description(
        base$pollutant, base$intervention, "trimester"
      )
    ))
    full_id <- paste0(base$intervention_id, "_full")
    reg[[full_id]] <- modifyList(base, list(
      intervention_id = full_id,
      output_stub = paste0(base$output_stub, "_full"),
      application_scope = "full",
      exposure_model = "aggregate_full",
      cox_formula = gform_period_cox_formula_text(base$pollutant, "full"),
      description = gform_period_intervention_description(
        base$pollutant, base$intervention, "full"
      )
    ))
  }
  reg
}

GFORM_PERIOD_INTERVENTION_REGISTRY <- build_gform_period_intervention_registry()

gform_period_pct_base_keys <- function() {
  c(
    "pm25_krg_pct10", "pm25_krg_pct20", "pm25_krg_pct30",
    "no2_krg_pct10", "no2_krg_pct20", "no2_krg_pct30",
    "o3_krg_pct10", "o3_krg_pct20", "o3_krg_pct30"
  )
}

gform_period_cap_base_keys <- function() {
  c(
    "pm25_krg_lt20", "pm25_krg_lt15", "pm25_krg_lt10", "pm25_krg_lt5",
    "no2_krg_lt20", "no2_krg_lt15", "no2_krg_lt10", "no2_krg_lt5"
  )
}

build_gform_period_intervention_order <- function() {
  pct <- gform_period_pct_base_keys()
  caps <- gform_period_cap_base_keys()
  with_scope <- function(key, scope) {
    paste0(GFORM_INTERVENTION_REGISTRY[[key]]$intervention_id, "_", scope)
  }
  c(
    vapply(pct, with_scope, character(1), scope = "tri", USE.NAMES = FALSE),
    vapply(pct, with_scope, character(1), scope = "full", USE.NAMES = FALSE),
    vapply(caps, with_scope, character(1), scope = "tri", USE.NAMES = FALSE),
    vapply(caps, with_scope, character(1), scope = "full", USE.NAMES = FALSE)
  )
}

GFORM_PERIOD_INTERVENTION_ORDER <- build_gform_period_intervention_order()

gform_period_intervention_keys <- function() {
  missing <- setdiff(
    GFORM_PERIOD_INTERVENTION_ORDER,
    names(GFORM_PERIOD_INTERVENTION_REGISTRY)
  )
  if (length(missing)) {
    stop(
      "GFORM_PERIOD_INTERVENTION_ORDER contiene IDs ausentes: ",
      paste(missing, collapse = ", ")
    )
  }
  GFORM_PERIOD_INTERVENTION_ORDER
}

print_gform_period_intervention_menu <- function(selected = NULL) {
  keys <- gform_period_intervention_keys()
  message("Intervenciones período (cambiar intervention_number en 12.2):")
  for (i in seq_along(keys)) {
    spec <- GFORM_PERIOD_INTERVENTION_REGISTRY[[keys[i]]]
    marker <- if (!is.null(selected) && i == selected) "  <-- seleccionada" else ""
    message(sprintf(
      "  %2d. %s [%s] — %s%s",
      i, spec$intervention_id, spec$application_scope,
      spec$description, marker
    ))
  }
  invisible(keys)
}

resolve_gform_period_intervention <- function(n) {
  keys <- gform_period_intervention_keys()
  n <- as.integer(n[[1L]])
  if (length(n) != 1L || is.na(n) || n < 1L || n > length(keys)) {
    stop(
      "intervention_number debe ser un entero entre 1 y ", length(keys),
      ". Ejecute print_gform_period_intervention_menu()."
    )
  }
  keys[[n]]
}

## Intervención sobre columnas de exposición agregada ----

gform_period_apply_intervention_cols <- function(df, var_names, intervention) {
  out <- df
  for (v in var_names) {
    if (!v %in% names(out)) next
    if (identical(intervention$type, "cap")) {
      out[[v]] <- pmin(out[[v]], intervention$cap, na.rm = FALSE)
    } else if (identical(intervention$type, "pct_reduce")) {
      out[[v]] <- out[[v]] * (1 - intervention$pct)
    } else if (!identical(intervention$type, "none")) {
      stop("Intervención desconocida: ", intervention$type)
    }
  }
  out
}

gform_period_extract_exposure <- function(births_df, pollutant, application_scope) {
  vars <- gform_period_exposure_var_names(pollutant, application_scope)
  missing <- setdiff(vars, names(births_df))
  if (length(missing)) {
    stop("Faltan variables de exposición: ", paste(missing, collapse = ", "))
  }
  births_df[, c("id", vars), drop = FALSE]
}

## Etapa 1: objetos RDS de intervención (nivel nacimiento) ----

generate_period_intervention <- function(
    births_df,
    spec,
    output_path,
    overwrite = FALSE) {

  if (file.exists(output_path) && !overwrite) {
    message("Intervención período existente (omitida): ", output_path)
    return(invisible(readRDS(output_path)))
  }

  vars <- gform_period_exposure_var_names(spec$pollutant, spec$application_scope)
  natural <- gform_period_extract_exposure(births_df, spec$pollutant, spec$application_scope)
  counterfactual <- gform_period_apply_intervention_cols(
    natural, vars, spec$intervention
  )

  intervention_obj <- list(
    intervention_id = spec$intervention_id,
    pollutant = spec$pollutant,
    application_scope = spec$application_scope,
    exposure_model = spec$exposure_model,
    intervention = spec$intervention,
    description = spec$description,
    exposure_vars = vars,
    exposure_natural = natural,
    exposure_intervention = counterfactual[, c("id", vars), drop = FALSE],
    created = Sys.time()
  )

  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(intervention_obj, output_path)
  invisible(intervention_obj)
}

build_all_period_interventions <- function(
    births_df,
    registry = GFORM_PERIOD_INTERVENTION_REGISTRY,
    output_dir,
    overwrite = FALSE,
    parallel = TRUE,
    n_workers = NULL) {

  build_one <- function(key) {
    spec <- registry[[key]]
    label <- paste0("RDS período — ", spec$intervention_id)
    block <- gform_time_block(label, {
      output_path <- file.path(output_dir, paste0(spec$intervention_id, ".rds"))
      generate_period_intervention(
        births_df = births_df,
        spec = spec,
        output_path = output_path,
        overwrite = overwrite
      )
      output_path
    })
    list(path = block$result, timing = block$timing)
  }

  registry_keys <- names(registry)
  if (parallel && length(registry_keys) > 1L && requireNamespace("furrr", quietly = TRUE)) {
    gform_setup_parallel(task = "build", n_workers = n_workers)
    built <- furrr::future_map(registry_keys, build_one, .options = gform_furrr_options())
  } else {
    built <- lapply(registry_keys, build_one)
  }

  timing <- gform_timing_log_init()
  for (item in built) {
    timing <- gform_timing_log_add(timing, item$timing)
  }

  list(
    paths = vapply(built, function(x) x$path, character(1)),
    timing = timing
  )
}

## Cox agregado por período (referencia 11.x, g-computación 10.1) ----

build_period_cox_frame <- function(
    data_base,
    exposure_df,
    control_vars = GFORM_PERIOD_CONTROL_VARS,
    dependent_var = GFORM_PERIOD_DEFAULTS$dependent_var,
    risk_entry_week = GFORM_PERIOD_DEFAULTS$risk_entry_week) {

  base_controls <- data.table::as.data.table(
    data_base[, c("id", "weeks", dependent_var, control_vars), drop = FALSE]
  )
  base_controls[, tstart := risk_entry_week]

  exp_dt <- data.table::as.data.table(exposure_df)
  data.table::setkey(base_controls, id)
  data.table::setkey(exp_dt, id)
  base_controls[exp_dt, on = "id"]
}

fit_period_cox_one_week <- function(
    rw,
    data_model_df,
    exposure_vars,
    control_vars,
    dependent_var,
    risk_entry_week) {

  if (rw < risk_entry_week) return(NULL)

  pred_terms <- exposure_vars[exposure_vars %in% names(data_model_df)]
  if (!length(pred_terms)) {
    warning("Sin predictores de exposición en semana ", rw)
    return(NULL)
  }

  rhs <- paste(c(pred_terms, control_vars), collapse = " + ")
  fml <- stats::as.formula(paste0(
    "Surv(tstart, weeks, ", dependent_var, ") ~ ", rhs
  ))

  vars_needed <- c("id", "weeks", "tstart", dependent_var, pred_terms, control_vars)
  model_df <- data_model_df[, vars_needed, drop = FALSE]
  model_df <- stats::na.omit(model_df)

  if (nrow(model_df) < 50L) {
    warning("Muy pocas filas para la semana ", rw, "; Cox período omitido.")
    return(NULL)
  }

  cox_fit <- tryCatch(
    survival::coxph(fml, data = model_df, x = FALSE, y = FALSE),
    error = function(e) NULL
  )
  if (is.null(cox_fit)) {
    warning("Cox período no convergió en semana ", rw)
    return(NULL)
  }

  slim_coxph_model_store_entry(list(
    model = cox_fit,
    formula = fml,
    id_order = model_df$id,
    pred_terms = pred_terms,
    control_vars = control_vars,
    baseline_hazard = survival::basehaz(cox_fit, centered = FALSE),
    risk_week = rw,
    exposure_model = "period_aggregate"
  ))
}

fit_period_natural_course_models <- function(
    data_base,
    exposure_natural,
    exposure_vars,
    risk_weeks_vec,
    control_vars = GFORM_PERIOD_CONTROL_VARS,
    dependent_var = GFORM_PERIOD_DEFAULTS$dependent_var,
    risk_entry_week = GFORM_PERIOD_DEFAULTS$risk_entry_week,
    parallel = TRUE) {

  cox_frame_natural <- build_period_cox_frame(
    data_base = data_base,
    exposure_df = exposure_natural,
    control_vars = control_vars,
    dependent_var = dependent_var,
    risk_entry_week = risk_entry_week
  )
  data_model_df <- as.data.frame(cox_frame_natural)

  weeks_to_fit <- risk_weeks_vec[risk_weeks_vec >= risk_entry_week]
  if (parallel && length(weeks_to_fit) > 1L && requireNamespace("furrr", quietly = TRUE)) {
    gform_setup_parallel(task = "cox")
    fitted <- furrr::future_map(
      weeks_to_fit,
      fit_period_cox_one_week,
      data_model_df = data_model_df,
      exposure_vars = exposure_vars,
      control_vars = control_vars,
      dependent_var = dependent_var,
      risk_entry_week = risk_entry_week,
      .options = gform_furrr_options()
    )
  } else {
    fitted <- lapply(weeks_to_fit, function(rw) {
      fit_period_cox_one_week(
        rw, data_model_df, exposure_vars, control_vars,
        dependent_var, risk_entry_week
      )
    })
  }

  model_store <- vector("list", length(risk_weeks_vec))
  names(model_store) <- as.character(risk_weeks_vec)
  for (i in seq_along(weeks_to_fit)) {
    model_store[[as.character(weeks_to_fit[i])]] <- fitted[[i]]
  }
  model_store <- slim_model_store(model_store)

  list(
    model_store = model_store,
    cox_frame_natural = cox_frame_natural,
    exposure_vars = exposure_vars
  )
}

gform_period_model_cache_path <- function(pollutant, application_scope, dir_models) {
  file.path(
    dir_models,
    paste0("natural_course_", pollutant, "_", application_scope, ".rds")
  )
}

gform_period_natural_course_cache_key <- function(
    data_base,
    pollutant,
    application_scope,
    exposure_vars,
    risk_weeks_vec,
    control_vars,
    dependent_var,
    risk_entry_week,
    sample_frac = NULL) {

  list(
    n_births = nrow(data_base),
    pollutant = pollutant,
    application_scope = application_scope,
    exposure_vars = exposure_vars,
    risk_weeks = risk_weeks_vec,
    control_vars = control_vars,
    dependent_var = dependent_var,
    risk_entry_week = risk_entry_week,
    sample_frac = sample_frac,
    pipeline = "gform_period_v1"
  )
}

load_or_fit_period_natural_course_models <- function(
    pollutant,
    application_scope,
    exposure_vars,
    dir_models,
    data_base,
    exposure_natural,
    risk_weeks_vec,
    control_vars = GFORM_PERIOD_CONTROL_VARS,
    dependent_var = GFORM_PERIOD_DEFAULTS$dependent_var,
    risk_entry_week = GFORM_PERIOD_DEFAULTS$risk_entry_week,
    sample_frac = NULL,
    parallel = FALSE,
    force_refit = FALSE) {

  dir.create(dir_models, recursive = TRUE, showWarnings = FALSE)
  cache_path <- gform_period_model_cache_path(pollutant, application_scope, dir_models)
  cache_key <- gform_period_natural_course_cache_key(
    data_base = data_base,
    pollutant = pollutant,
    application_scope = application_scope,
    exposure_vars = exposure_vars,
    risk_weeks_vec = risk_weeks_vec,
    control_vars = control_vars,
    dependent_var = dependent_var,
    risk_entry_week = risk_entry_week,
    sample_frac = sample_frac
  )

  if (!force_refit && file.exists(cache_path)) {
    cached <- tryCatch(readRDS(cache_path), error = function(e) NULL)
    if (!is.null(cached) && identical(cached$cache_key, cache_key)) {
      message("Modelos Cox período cacheados: ", cache_path)
      cached$model_store <- slim_model_store(cached$model_store)
      return(list(
        model_store = cached$model_store,
        cox_frame_natural = cached$cox_frame_natural,
        exposure_vars = cached$exposure_vars,
        from_cache = TRUE
      ))
    }
    if (!is.null(cached)) {
      message("Cache Cox período desactualizado; re-ajustando...")
    }
  }

  fit <- fit_period_natural_course_models(
    data_base = data_base,
    exposure_natural = exposure_natural,
    exposure_vars = exposure_vars,
    risk_weeks_vec = risk_weeks_vec,
    control_vars = control_vars,
    dependent_var = dependent_var,
    risk_entry_week = risk_entry_week,
    parallel = parallel
  )

  cache_obj <- c(
    fit,
    list(cache_key = cache_key, fitted_at = Sys.time(), slimmed = TRUE)
  )
  tmp_path <- paste0(cache_path, ".tmp")
  saveRDS(cache_obj, tmp_path)
  file.rename(tmp_path, cache_path)
  message("Modelos Cox período guardados: ", cache_path)
  message("  Fórmula: ", gform_period_cox_formula_text(pollutant, application_scope))
  c(fit, list(from_cache = FALSE))
}

## Orquestador por intervención (sin heatmap semanal) ----

run_period_gform_intervention <- function(
    intervention_spec,
    intervention_path,
    data_base,
    model_store,
    person_weeks,
    risk_weeks_vec,
    control_vars,
    total_births,
    cox_frame_natural = NULL,
    exposure_natural = NULL,
    dependent_var = GFORM_PERIOD_DEFAULTS$dependent_var,
    risk_entry_week = GFORM_PERIOD_DEFAULTS$risk_entry_week,
    follow_up_weeks = GFORM_PERIOD_DEFAULTS$follow_up_weeks,
    boot_iter = GFORM_PERIOD_DEFAULTS$boot_iter,
    boot_seed = GFORM_PERIOD_DEFAULTS$boot_seed,
    target_week = GFORM_PERIOD_DEFAULTS$population_week,
    sample_frac = NULL,
    run_bootstrap = FALSE,
    parallel_bootstrap = FALSE,
    dir_bootstrap = NULL,
    bootstrap_resume = TRUE,
    minimal_output = GFORM_PERIOD_DEFAULTS$minimal_output,
    dir_weekly = NULL,
    dir_population = NULL,
    dir_other = NULL) {

  timing <- gform_timing_log_init()

  block <- gform_time_block("Construir frame Cox intervención (período)", {
    intervention_obj <- readRDS(intervention_path)
    exposure_vars <- intervention_obj$exposure_vars

    if (is.null(cox_frame_natural)) {
      if (is.null(exposure_natural)) {
        exposure_natural <- intervention_obj$exposure_natural
      }
      cox_frame_natural <- build_period_cox_frame(
        data_base = data_base,
        exposure_df = exposure_natural,
        control_vars = control_vars,
        dependent_var = dependent_var,
        risk_entry_week = risk_entry_week
      )
    }

    cox_frame_intervention <- build_period_cox_frame(
      data_base = data_base,
      exposure_df = intervention_obj$exposure_intervention,
      control_vars = control_vars,
      dependent_var = dependent_var,
      risk_entry_week = risk_entry_week
    )
    rm(intervention_obj)
    gc()

    list(
      exposure_vars = exposure_vars,
      cox_frame_natural = cox_frame_natural,
      cox_frame_intervention = cox_frame_intervention
    )
  })
  timing <- gform_timing_log_add(timing, block$timing)
  cox_frame_natural <- block$result$cox_frame_natural
  cox_frame_intervention <- block$result$cox_frame_intervention

  fingerprint <- gform_run_fingerprint(
    output_stub = intervention_spec$output_stub,
    boot_iter = boot_iter,
    boot_seed = boot_seed,
    total_births = total_births,
    sample_frac = sample_frac
  )
  point_ck_path <- if (!is.null(dir_bootstrap)) {
    gform_point_checkpoint_path(intervention_spec$output_stub, dir_bootstrap)
  } else {
    NULL
  }
  skip_point_estimate <- FALSE
  if (isTRUE(bootstrap_resume) && !is.null(point_ck_path)) {
    pt_ck <- gform_read_point_checkpoint(point_ck_path, fingerprint)
    if (!is.null(pt_ck)) {
      skip_point_estimate <- TRUE
      message("Reanudación período: punto estimado en disco; omitiendo predicciones.")
    }
  }

  if (isTRUE(skip_point_estimate)) {
    pt_ck <- gform_read_point_checkpoint(point_ck_path, fingerprint)
    population_effects <- pt_ck$population_effects
    weekly_effects <- if (isTRUE(minimal_output)) NULL else pt_ck$weekly_effects
    cumulative_risk_curves <- if (isTRUE(minimal_output)) NULL else pt_ck$cumulative_risk_curves
    nat_mean <- pt_ck$nat_mean
  } else {
    block <- gform_time_block("Predicción hazards — curso natural (período)", {
      predict_weekly_hazards(
        model_store = model_store,
        person_weeks = person_weeks,
        risk_weeks_vec = risk_weeks_vec,
        cox_frame = cox_frame_natural,
        risk_entry_week = risk_entry_week
      )
    })
    timing <- gform_timing_log_add(timing, block$timing)
    prob_natural <- block$result

    block <- gform_time_block("Predicción hazards — intervención (período)", {
      predict_weekly_hazards(
        model_store = model_store,
        person_weeks = person_weeks,
        risk_weeks_vec = risk_weeks_vec,
        cox_frame = cox_frame_intervention,
        risk_entry_week = risk_entry_week
      )
    })
    timing <- gform_timing_log_add(timing, block$timing)
    prob_intervention <- block$result

    block <- gform_time_block("Efectos puntuales (período)", {
      if (isTRUE(minimal_output)) {
        population_effects <- compute_population_effects(
          prob_natural = prob_natural,
          prob_intervention = prob_intervention,
          total_births = total_births,
          target_week = target_week
        )
        list(
          nat_mean = NULL,
          weekly_effects = NULL,
          cumulative_risk_curves = NULL,
          population_effects = population_effects
        )
      } else {
        surv_nat <- compute_survival(prob_natural)
        nat_mean <- surv_nat[, .(risk_natural = mean(risk, na.rm = TRUE)), by = time]
        weekly_effects <- compute_weekly_effects(prob_natural, prob_intervention)
        cumulative_risk_curves <- compute_cumulative_risk_curves_global(
          weekly_effects, follow_up_weeks = follow_up_weeks
        )
        population_effects <- compute_population_effects(
          prob_natural = prob_natural,
          prob_intervention = prob_intervention,
          total_births = total_births,
          target_week = target_week
        )
        list(
          nat_mean = nat_mean,
          weekly_effects = weekly_effects,
          cumulative_risk_curves = cumulative_risk_curves,
          population_effects = population_effects
        )
      }
    })
    timing <- gform_timing_log_add(timing, block$timing)
    nat_mean <- block$result$nat_mean
    weekly_effects <- block$result$weekly_effects
    cumulative_risk_curves <- block$result$cumulative_risk_curves
    population_effects <- block$result$population_effects
    rm(prob_natural, prob_intervention)
    gc()

    if (run_bootstrap && boot_iter > 0L && !is.null(dir_bootstrap)) {
      gform_save_point_checkpoint(
        path = point_ck_path,
        fingerprint = fingerprint,
        weekly_effects = weekly_effects,
        population_effects = population_effects,
        cumulative_risk_curves = cumulative_risk_curves,
        nat_mean = nat_mean
      )
    }
  }

  boot_out <- NULL
  weekly_ci <- weekly_effects
  population_ci <- population_effects

  if (run_bootstrap && boot_iter > 0L) {
    if (is.null(dir_bootstrap)) {
      stop("run_period_gform_intervention: se requiere dir_bootstrap.")
    }
    block <- gform_time_block(paste0("Bootstrap paramétrico (", boot_iter, " iter)"), {
      boot_out <- bootstrap_gformula_effects(
        model_store = model_store,
        person_weeks = person_weeks,
        cox_frame_natural = cox_frame_natural,
        cox_frame_intervention = cox_frame_intervention,
        risk_weeks_vec = risk_weeks_vec,
        total_births = total_births,
        output_stub = intervention_spec$output_stub,
        dir_bootstrap = dir_bootstrap,
        dependent_var = dependent_var,
        risk_entry_week = risk_entry_week,
        boot_iter = boot_iter,
        boot_seed = boot_seed,
        target_week = target_week,
        sample_frac = sample_frac,
        resume = bootstrap_resume,
        parallel = parallel_bootstrap
      )
      weekly_ci <- if (is.null(weekly_effects)) {
        NULL
      } else {
        dplyr::left_join(weekly_effects, boot_out$weekly_ci, by = "week")
      }
      population_ci <- dplyr::left_join(
        population_effects,
        boot_out$population_ci |>
          dplyr::select(
            "scenario",
            dplyr::ends_with("_lcl"),
            dplyr::ends_with("_ucl")
          ),
        by = "scenario"
      )
      boot_out$weekly_boot <- NULL
      boot_out$population_boot <- NULL
      gc(verbose = FALSE)
      list(
        boot_out = boot_out,
        weekly_ci = weekly_ci,
        population_ci = population_ci
      )
    })
    timing <- gform_timing_log_add(timing, block$timing)
    boot_out <- block$result$boot_out
    weekly_ci <- block$result$weekly_ci
    population_ci <- block$result$population_ci
  }

  rm(cox_frame_intervention)
  gc(verbose = FALSE)

  list(
    intervention_spec = intervention_spec,
    weekly_effects = weekly_ci,
    population_effects = population_ci,
    cumulative_risk_curves = cumulative_risk_curves,
    singleweek_intervention_heatmap = NULL,
    bootstrap = boot_out,
    timing = timing
  )
}

gform_period_infer_application_scope <- function(output_stub) {
  if (grepl("_tri$", output_stub)) {
    return("trimester")
  }
  if (grepl("_full$", output_stub)) {
    return("full")
  }
  NA_character_
}

gform_period_read_excel_application_scope <- function(excel_path) {
  if (!file.exists(excel_path)) {
    return(NA_character_)
  }
  sheets <- tryCatch(readxl::excel_sheets(excel_path), error = function(e) character())
  if (!"metadata" %in% sheets) {
    return(NA_character_)
  }
  meta <- tryCatch(
    readxl::read_excel(excel_path, sheet = "metadata"),
    error = function(e) NULL
  )
  if (is.null(meta) || !all(c("field", "value") %in% names(meta))) {
    return(NA_character_)
  }
  as.character(meta$value[meta$field == "application_scope"][1L])
}

gform_period_stub_matches_scope <- function(stub, expected_scope, excel_path) {
  if (!identical(expected_scope, gform_period_infer_application_scope(stub))) {
    scope <- gform_period_read_excel_application_scope(excel_path)
    if (is.na(scope) || !nzchar(scope)) {
      return(FALSE)
    }
    return(identical(scope, expected_scope))
  }
  TRUE
}

gform_period_legacy_tri_stub <- function(output_stub) {
  if (grepl("_tri$", output_stub)) {
    sub("_tri$", "", output_stub)
  } else {
    NULL
  }
}

gform_period_output_stub_variants <- function(output_stub) {
  legacy <- gform_period_legacy_tri_stub(output_stub)
  if (is.null(legacy)) {
    output_stub
  } else {
    unique(c(output_stub, legacy))
  }
}

gform_period_intervention_rds_path <- function(intervention_id, dir_interventions) {
  primary <- file.path(dir_interventions, paste0(intervention_id, ".rds"))
  if (file.exists(primary)) {
    return(primary)
  }
  if (grepl("_tri$", intervention_id)) {
    legacy_id <- sub("_tri$", "", intervention_id)
    legacy_path <- file.path(dir_interventions, paste0(legacy_id, ".rds"))
    if (file.exists(legacy_path)) {
      message("RDS legacy (sin _tri): ", legacy_path)
      return(legacy_path)
    }
  }
  primary
}

gform_period_bootstrap_ck_ok <- function(
    ck,
    output_stub,
    boot_iter,
    boot_seed,
    total_births,
    sample_frac) {

  if (is.null(ck) || is.null(ck$last_completed) || ck$last_completed < boot_iter) {
    return(FALSE)
  }
  for (stub in gform_period_output_stub_variants(output_stub)) {
    fp <- gform_run_fingerprint(stub, boot_iter, boot_seed, total_births, sample_frac)
    if (gform_bootstrap_ck_matches(ck, fp)) {
      return(TRUE)
    }
  }
  FALSE
}

gform_period_intervention_is_complete <- function(
    output_stub,
    dir_population,
    dir_summary,
    dir_bootstrap = NULL,
    boot_iter = GFORM_PERIOD_DEFAULTS$boot_iter,
    boot_seed = GFORM_PERIOD_DEFAULTS$boot_seed,
    total_births = NULL,
    sample_frac = NULL,
    run_bootstrap = FALSE,
    expected_application_scope = NULL) {

  if (is.null(expected_application_scope)) {
    expected_application_scope <- gform_period_infer_application_scope(output_stub)
  }

  for (stub in gform_period_output_stub_variants(output_stub)) {
    population_path <- file.path(dir_population, paste0(stub, "_population_effects.rds"))
    excel_path <- file.path(dir_summary, paste0(stub, "_point_estimates.xlsx"))

    if (!all(file.exists(c(population_path, excel_path)))) {
      next
    }
    if (!gform_period_stub_matches_scope(stub, expected_application_scope, excel_path)) {
      next
    }

    if (isTRUE(run_bootstrap) && boot_iter > 0L && !is.null(dir_bootstrap)) {
      boot_ok <- FALSE
      for (boot_stub in gform_period_output_stub_variants(output_stub)) {
        boot_ck_path <- gform_bootstrap_paths(boot_stub, dir_bootstrap)$checkpoint
        if (!file.exists(boot_ck_path)) next
        ck <- tryCatch(readRDS(boot_ck_path), error = function(e) NULL)
        if (gform_period_bootstrap_ck_ok(
          ck, output_stub, boot_iter, boot_seed, total_births, sample_frac
        )) {
          boot_ok <- TRUE
          break
        }
      }
      if (!boot_ok) next
    }

    return(TRUE)
  }

  FALSE
}

gform_period_output_dirs <- function(data_out_g = "02_Output/G-Form-Period") {
  list(
    population = file.path(data_out_g, "PopulationEffects"),
    summary = file.path(data_out_g, "Summary_results"),
    bootstrap = file.path(data_out_g, "Bootstrap")
  )
}

gform_period_missing_intervention_numbers <- function(
    data_out_g = "02_Output/G-Form-Period",
    boot_iter = GFORM_PERIOD_DEFAULTS$boot_iter,
    boot_seed = GFORM_PERIOD_DEFAULTS$boot_seed,
    run_bootstrap = FALSE) {

  dirs <- gform_period_output_dirs(data_out_g)
  keys <- gform_period_intervention_keys()
  missing <- integer()

  for (i in seq_along(keys)) {
    spec <- GFORM_PERIOD_INTERVENTION_REGISTRY[[keys[i]]]
    if (gform_period_intervention_is_complete(
      output_stub = spec$output_stub,
      expected_application_scope = spec$application_scope,
      dir_population = dirs$population,
      dir_summary = dirs$summary,
      dir_bootstrap = dirs$bootstrap,
      boot_iter = boot_iter,
      boot_seed = boot_seed,
      run_bootstrap = run_bootstrap
    )) {
      next
    }
    missing <- c(missing, i)
  }

  missing
}

print_gform_period_missing_interventions <- function(
    data_out_g = "02_Output/G-Form-Period",
    run_bootstrap = FALSE) {

  boot_iter <- if (isTRUE(run_bootstrap)) 250L else 0L
  missing <- gform_period_missing_intervention_numbers(
    data_out_g = data_out_g,
    boot_iter = boot_iter,
    run_bootstrap = run_bootstrap
  )
  keys <- gform_period_intervention_keys()

  message("\n=== Intervenciones período pendientes ===")
  if (!length(missing)) {
    message("Ninguna (34/34 completas con criterio actual).")
    return(invisible(integer()))
  }

  message("Faltan ", length(missing), " / ", length(keys), ":")
  for (i in missing) {
    spec <- GFORM_PERIOD_INTERVENTION_REGISTRY[[keys[i]]]
    message(sprintf(
      "  %2d. %s -> %s [%s]",
      i, spec$intervention_id, spec$output_stub, spec$application_scope
    ))
  }
  message("\nIDs para GFORM_PERIOD_INTERVENTIONS=",
          paste(missing, collapse = ","))
  invisible(missing)
}

save_period_results <- function(
    population_effects,
    population_path,
    metadata = list()) {

  obj <- c(list(point = population_effects), metadata)
  dir.create(dirname(population_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(obj, population_path)
  invisible(population_path)
}

save_period_excel <- function(
    population_effects,
    excel_path,
    metadata = NULL) {

  sheets <- list(population_effects = population_effects)
  if (!is.null(metadata) && length(metadata)) {
    meta_df <- data.frame(
      field = names(metadata),
      value = vapply(metadata, function(x) {
        if (is.null(x)) ""
        else if (length(x) > 1L) paste(x, collapse = ", ")
        else as.character(x)
      }, character(1)),
      stringsAsFactors = FALSE
    )
    sheets$metadata <- meta_df
  }
  dir.create(dirname(excel_path), recursive = TRUE, showWarnings = FALSE)
  writexl::write_xlsx(sheets, path = excel_path)
  invisible(excel_path)
}

get_period_exposure_natural <- function(births_df, pollutant, application_scope, cache) {
  cache_key <- paste(pollutant, application_scope, sep = "__")
  if (!is.null(cache[[cache_key]])) {
    return(cache[[cache_key]])
  }
  cache[[cache_key]] <- gform_period_extract_exposure(
    births_df, pollutant, application_scope
  )
  cache[[cache_key]]
}

prepare_period_births_cohort <- function(data_inp) {
  data_inp <- sub("/+$", "", data_inp)
  lagged_path <- file.path(data_inp, "births_2010_2020_exposure_weeks_lagged.RData")
  metrics_path <- file.path(data_inp, "births_exposure_period_metrics_full_d30_d4_tri.RData")

  if (!file.exists(lagged_path)) {
    stop(
      "Falta ", basename(lagged_path),
      " (mismo archivo que usa 10.2 G-Form_models.R)."
    )
  }
  if (!file.exists(metrics_path)) {
    stop(
      "Falta ", basename(metrics_path),
      " (mismo archivo que usa 10.2 G-Form_models.R)."
    )
  }

  period_expo <- rio::import(metrics_path) |>
    dplyr::select(
      "id",
      dplyr::matches("_krg_(full|t1|t2|t3)$"),
      dplyr::any_of(c("tad_full", "ndvi_full"))
    ) |>
    dplyr::distinct(.data$id, .keep_all = TRUE)

  load(lagged_path)

  controls_from_long <- setdiff(
    GFORM_PERIOD_CONTROL_VARS,
    c("tad_full", "ndvi_full")
  )

  births_base <- data_long |>
    dplyr::select(dplyr::any_of(c(
      "id", "weeks", GFORM_PERIOD_DEFAULTS$dependent_var, controls_from_long
    ))) |>
    dplyr::distinct(.data$id, .keep_all = TRUE) |>
    dplyr::filter(!is.na(.data$weeks), .data$weeks >= 28L)

  out <- births_base |>
    dplyr::left_join(period_expo, by = "id") |>
    dplyr::mutate(
      tstart = 27L,
      month_week1 = factor(.data$month_week1),
      year_week1 = factor(.data$year_week1),
      covid = factor(.data$covid)
    )

  rm(data_long)
  gc(verbose = FALSE)

  missing_controls <- setdiff(GFORM_PERIOD_CONTROL_VARS, names(out))
  if (length(missing_controls)) {
    stop("Faltan covariables período: ", paste(missing_controls, collapse = ", "))
  }

  out
}
