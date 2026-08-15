# 16.0 G-Formula por períodos multicontaminante — funciones ----
#
# Extensión de 12.0: un Cox multicontaminante por alcance temporal
#   - trimester: pm25_t1+t2+t3 + no2_t1+t2+t3 + o3_t1+t2+t3 + covariables
#   - full:      pm25_full + no2_full + o3_full + covariables
#
# Seis intervenciones pct20 (3 trimestre + 3 overall): en cada una se reduce
# el contaminante focal y los co-contaminantes permanecen en curso natural.
# El mismo modelo Cox cacheado sirve para las 3 intervenciones del alcance.

source("00_Code/12.0 G-Form_period_functions.R")

GFORM_PERIOD_MULTI_DATA_OUT <- "02_Output/G-Form-Period-Multi/"

GFORM_PERIOD_MULTI_POLLUTANTS <- c("pm25_krg", "no2_krg", "o3_krg")

GFORM_PERIOD_MULTI_DEFAULTS <- utils::modifyList(GFORM_PERIOD_DEFAULTS, list(
  output_cumulative_risk_suffix = "curvas_riesgo_acumulado_periodo_multi"
))

gform_period_multi_pollutant_short <- function(pollutant) {
  sub("_krg$", "", pollutant)
}

gform_period_multi_all_exposure_vars <- function(application_scope) {
  pollutants_short <- c("pm25", "no2", "o3")
  if (identical(application_scope, "trimester")) {
    unlist(lapply(pollutants_short, function(p) {
      paste0(p, "_krg_t", 1:3)
    }), use.names = FALSE)
  } else if (identical(application_scope, "full")) {
    paste0(pollutants_short, "_krg_full")
  } else {
    stop("application_scope desconocido: ", application_scope)
  }
}

gform_period_multi_target_exposure_vars <- function(pollutant, application_scope) {
  short <- gform_period_multi_pollutant_short(pollutant)
  if (identical(application_scope, "trimester")) {
    paste0(short, "_krg_t", 1:3)
  } else if (identical(application_scope, "full")) {
    paste0(short, "_krg_full")
  } else {
    stop("application_scope desconocido: ", application_scope)
  }
}

gform_period_multi_cox_formula_text <- function(application_scope) {
  dep <- GFORM_PERIOD_MULTI_DEFAULTS$dependent_var
  cov <- paste(GFORM_PERIOD_CONTROL_VARS, collapse = " + ")
  expo <- paste(gform_period_multi_all_exposure_vars(application_scope), collapse = " + ")
  suffix <- if (identical(application_scope, "trimester")) {
    "  [multicontaminante: T1+T2+T3 pm25+no2+o3 mutuamente ajustados]"
  } else {
    "  [multicontaminante: full pm25+no2+o3 mutuamente ajustados]"
  }
  paste0("Surv(tstart, weeks, ", dep, ") ~ ", expo, " + ", cov, suffix)
}

gform_period_multi_intervention_description <- function(pollutant, application_scope) {
  pol <- gform_period_pollutant_label(pollutant)
  rule <- paste0(pol, " reducción 20%")
  if (identical(application_scope, "trimester")) {
    paste0(
      rule,
      " — misma regla en medias T1, T2 y T3 de ", pol,
      "; co-contaminantes en curso natural; ",
      "Cox trimestral multicontaminante (un modelo compartido)"
    )
  } else if (identical(application_scope, "full")) {
    paste0(
      rule,
      " — media del período completo de ", pol,
      "; co-contaminantes en curso natural; ",
      "Cox overall multicontaminante (un modelo compartido)"
    )
  } else {
    stop("application_scope desconocido: ", application_scope)
  }
}

build_gform_period_multi_intervention_registry <- function() {
  reg <- list()
  intervention <- list(type = "pct_reduce", pct = 0.20)

  for (pol in GFORM_PERIOD_MULTI_POLLUTANTS) {
    short <- gform_period_multi_pollutant_short(pol)
    for (scope_info in list(
      list(scope = "trimester", suffix = "tri"),
      list(scope = "full", suffix = "full")
    )) {
      scope <- scope_info$scope
      suffix <- scope_info$suffix
      id <- paste0(short, "_krg_pct20_", suffix, "_multi")
      reg[[id]] <- list(
        intervention_id = id,
        output_stub = paste0(short, "_pct20_", suffix, "_multi"),
        pollutant = pol,
        application_scope = scope,
        exposure_model = if (identical(scope, "trimester")) {
          "aggregate_trimester_multi"
        } else {
          "aggregate_full_multi"
        },
        intervention = intervention,
        target_exposure_vars = gform_period_multi_target_exposure_vars(pol, scope),
        exposure_vars = gform_period_multi_all_exposure_vars(scope),
        cox_formula = gform_period_multi_cox_formula_text(scope),
        description = gform_period_multi_intervention_description(pol, scope)
      )
    }
  }

  reg
}

GFORM_PERIOD_MULTI_INTERVENTION_REGISTRY <- build_gform_period_multi_intervention_registry()

GFORM_PERIOD_MULTI_INTERVENTION_ORDER <- names(GFORM_PERIOD_MULTI_INTERVENTION_REGISTRY)

gform_period_multi_intervention_keys <- function() {
  missing <- setdiff(
    GFORM_PERIOD_MULTI_INTERVENTION_ORDER,
    names(GFORM_PERIOD_MULTI_INTERVENTION_REGISTRY)
  )
  if (length(missing)) {
    stop(
      "GFORM_PERIOD_MULTI_INTERVENTION_ORDER contiene IDs ausentes: ",
      paste(missing, collapse = ", ")
    )
  }
  GFORM_PERIOD_MULTI_INTERVENTION_ORDER
}

print_gform_period_multi_cox_specs <- function() {
  message("\n=== Especificaciones Cox multicontaminante (sección 16) ===")
  message("Dos modelos cacheados (Models/):")
  message("  • Trimestre: ", gform_period_multi_cox_formula_text("trimester"))
  message("  • Overall:   ", gform_period_multi_cox_formula_text("full"))
  message("Seis intervenciones pct20: 3 trimestre + 3 overall (comparten Cox por alcance).")
  message("Salida mínima: efectos poblacionales (prevalencia, riesgo, RD, RR, PAF).")
  message("Bootstrap (IC): opcional — GFORM_RUN_BOOTSTRAP=true GFORM_BOOT_ITER=250")
  invisible(NULL)
}

print_gform_period_multi_intervention_menu <- function(selected = NULL) {
  keys <- gform_period_multi_intervention_keys()
  message("Intervenciones período multicontaminante (cambiar intervention_number en 16.2):")
  for (i in seq_along(keys)) {
    spec <- GFORM_PERIOD_MULTI_INTERVENTION_REGISTRY[[keys[i]]]
    marker <- if (!is.null(selected) && i == selected) "  <-- seleccionada" else ""
    message(sprintf(
      "  %d. %s [%s] focal=%s — %s%s",
      i, spec$intervention_id, spec$application_scope,
      gform_period_pollutant_label(spec$pollutant),
      spec$description, marker
    ))
  }
  invisible(keys)
}

resolve_gform_period_multi_intervention <- function(n) {
  keys <- gform_period_multi_intervention_keys()
  n <- as.integer(n[[1L]])
  if (length(n) != 1L || is.na(n) || n < 1L || n > length(keys)) {
    stop(
      "intervention_number debe ser un entero entre 1 y ", length(keys),
      ". Ejecute print_gform_period_multi_intervention_menu()."
    )
  }
  keys[[n]]
}

gform_period_multi_extract_exposure <- function(births_df, application_scope) {
  vars <- gform_period_multi_all_exposure_vars(application_scope)
  missing <- setdiff(vars, names(births_df))
  if (length(missing)) {
    stop("Faltan variables de exposición: ", paste(missing, collapse = ", "))
  }
  births_df[, c("id", vars), drop = FALSE]
}

generate_period_multi_intervention <- function(
    births_df,
    spec,
    output_path,
    overwrite = FALSE) {

  if (file.exists(output_path) && !overwrite) {
    message("Intervención período multi existente (omitida): ", output_path)
    return(invisible(readRDS(output_path)))
  }

  natural <- gform_period_multi_extract_exposure(births_df, spec$application_scope)
  counterfactual <- gform_period_apply_intervention_cols(
    natural,
    spec$target_exposure_vars,
    spec$intervention
  )

  intervention_obj <- list(
    intervention_id = spec$intervention_id,
    pollutant = spec$pollutant,
    application_scope = spec$application_scope,
    exposure_model = spec$exposure_model,
    intervention = spec$intervention,
    description = spec$description,
    exposure_vars = spec$exposure_vars,
    target_exposure_vars = spec$target_exposure_vars,
    exposure_natural = natural,
    exposure_intervention = counterfactual[, c("id", spec$exposure_vars), drop = FALSE],
    created = Sys.time()
  )

  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(intervention_obj, output_path)
  invisible(intervention_obj)
}

build_all_period_multi_interventions <- function(
    births_df,
    registry = GFORM_PERIOD_MULTI_INTERVENTION_REGISTRY,
    output_dir,
    overwrite = FALSE,
    parallel = TRUE,
    n_workers = NULL) {

  build_one <- function(key) {
    spec <- registry[[key]]
    label <- paste0("RDS período multi — ", spec$intervention_id)
    block <- gform_time_block(label, {
      output_path <- file.path(output_dir, paste0(spec$intervention_id, ".rds"))
      generate_period_multi_intervention(
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

gform_period_multi_model_cache_path <- function(application_scope, dir_models) {
  file.path(
    dir_models,
    paste0("natural_course_period_multi_", application_scope, ".rds")
  )
}

gform_period_multi_natural_course_cache_key <- function(
    data_base,
    application_scope,
    exposure_vars,
    risk_weeks_vec,
    control_vars,
    dependent_var,
    risk_entry_week,
    sample_frac = NULL) {

  list(
    n_births = nrow(data_base),
    application_scope = application_scope,
    exposure_vars = exposure_vars,
    risk_weeks = risk_weeks_vec,
    control_vars = control_vars,
    dependent_var = dependent_var,
    risk_entry_week = risk_entry_week,
    sample_frac = sample_frac,
    pipeline = "gform_period_multi_v1"
  )
}

load_or_fit_period_multi_natural_course_models <- function(
    application_scope,
    exposure_vars,
    dir_models,
    data_base,
    exposure_natural,
    risk_weeks_vec,
    control_vars = GFORM_PERIOD_CONTROL_VARS,
    dependent_var = GFORM_PERIOD_MULTI_DEFAULTS$dependent_var,
    risk_entry_week = GFORM_PERIOD_MULTI_DEFAULTS$risk_entry_week,
    sample_frac = NULL,
    parallel = FALSE,
    force_refit = FALSE) {

  dir.create(dir_models, recursive = TRUE, showWarnings = FALSE)
  cache_path <- gform_period_multi_model_cache_path(application_scope, dir_models)
  cache_key <- gform_period_multi_natural_course_cache_key(
    data_base = data_base,
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
      message("Modelos Cox período multi cacheados: ", cache_path)
      cached$model_store <- slim_model_store(cached$model_store)
      return(list(
        model_store = cached$model_store,
        cox_frame_natural = cached$cox_frame_natural,
        exposure_vars = cached$exposure_vars,
        from_cache = TRUE
      ))
    }
    if (!is.null(cached)) {
      message("Cache Cox período multi desactualizado; re-ajustando...")
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
  message("Modelos Cox período multi guardados: ", cache_path)
  message("  Fórmula: ", gform_period_multi_cox_formula_text(application_scope))
  c(fit, list(from_cache = FALSE))
}

get_period_multi_exposure_natural <- function(births_df, application_scope, cache) {
  cache_key <- paste("multi", application_scope, sep = "__")
  if (!is.null(cache[[cache_key]])) {
    return(cache[[cache_key]])
  }
  cache[[cache_key]] <- gform_period_multi_extract_exposure(births_df, application_scope)
  cache[[cache_key]]
}

gform_period_multi_output_dirs <- function(
    data_out_g = GFORM_PERIOD_MULTI_DATA_OUT) {
  list(
    population = file.path(data_out_g, "PopulationEffects"),
    summary = file.path(data_out_g, "Summary_results"),
    bootstrap = file.path(data_out_g, "Bootstrap")
  )
}

gform_period_multi_intervention_is_complete <- function(
    output_stub,
    dir_population,
    dir_summary,
    dir_bootstrap = NULL,
    boot_iter = GFORM_PERIOD_MULTI_DEFAULTS$boot_iter,
    boot_seed = GFORM_PERIOD_MULTI_DEFAULTS$boot_seed,
    total_births = NULL,
    sample_frac = NULL,
    run_bootstrap = FALSE,
    expected_application_scope = NULL) {

  gform_period_intervention_is_complete(
    output_stub = output_stub,
    dir_population = dir_population,
    dir_summary = dir_summary,
    dir_bootstrap = dir_bootstrap,
    boot_iter = boot_iter,
    boot_seed = boot_seed,
    total_births = total_births,
    sample_frac = sample_frac,
    run_bootstrap = run_bootstrap,
    expected_application_scope = expected_application_scope
  )
}

gform_period_multi_missing_intervention_numbers <- function(
    data_out_g = GFORM_PERIOD_MULTI_DATA_OUT,
    boot_iter = GFORM_PERIOD_MULTI_DEFAULTS$boot_iter,
    boot_seed = GFORM_PERIOD_MULTI_DEFAULTS$boot_seed,
    run_bootstrap = FALSE) {

  dirs <- gform_period_multi_output_dirs(data_out_g)
  keys <- gform_period_multi_intervention_keys()
  missing <- integer()

  for (i in seq_along(keys)) {
    spec <- GFORM_PERIOD_MULTI_INTERVENTION_REGISTRY[[keys[i]]]
    if (gform_period_multi_intervention_is_complete(
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

print_gform_period_multi_missing_interventions <- function(
    data_out_g = GFORM_PERIOD_MULTI_DATA_OUT,
    run_bootstrap = FALSE) {

  boot_iter <- if (isTRUE(run_bootstrap)) 250L else 0L
  missing <- gform_period_multi_missing_intervention_numbers(
    data_out_g = data_out_g,
    boot_iter = boot_iter,
    run_bootstrap = run_bootstrap
  )
  keys <- gform_period_multi_intervention_keys()

  message("\n=== Intervenciones período multi pendientes ===")
  if (!length(missing)) {
    message("Ninguna (6/6 completas con criterio actual).")
    return(invisible(integer()))
  }

  message("Faltan ", length(missing), " / ", length(keys), ":")
  for (i in missing) {
    spec <- GFORM_PERIOD_MULTI_INTERVENTION_REGISTRY[[keys[i]]]
    message(sprintf(
      "  %d. %s -> %s [%s]",
      i, spec$intervention_id, spec$output_stub, spec$application_scope
    ))
  }
  message("\nIDs para GFORM_PERIOD_MULTI_INTERVENTIONS=",
          paste(missing, collapse = ","))
  invisible(missing)
}

gform_period_multi_intervention_rds_path <- function(intervention_id, dir_interventions) {
  file.path(dir_interventions, paste0(intervention_id, ".rds"))
}
