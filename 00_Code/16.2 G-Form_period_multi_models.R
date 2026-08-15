# 16.2 G-Formula por períodos multicontaminante — modelos ajustados ----
#
# Requisito previo: 16.1 G-Form_period_multi_build_interventions.R
#
# Modelos Cox cacheados (2):
#   - trimester: T1+T2+T3 de pm25+no2+o3 mutuamente ajustados + covariables
#   - full:      full de pm25+no2+o3 mutuamente ajustados + covariables
#
# Seis intervenciones pct20 comparten el Cox del alcance correspondiente.
#
# Uso:
#   GFORM_EXEC_MODE=server Rscript "00_Code/16.2 G-Form_period_multi_models.R"
#   GFORM_PERIOD_MULTI_INTERVENTIONS=1,2,3,4,5,6 GFORM_EXEC_MODE=server Rscript "00_Code/16.2 G-Form_period_multi_models.R"

rm(list = ls())
gc(verbose = FALSE)

source("00_Code/0.2 Packages_gform.R")
source("00_Code/0.1 Settings.R")
source("00_Code/0.3 Functions.R")
source("00_Code/16.0 G-Form_period_multi_functions.R")

data_inp <- "01_Data/Output"
data_out_g <- GFORM_PERIOD_MULTI_DATA_OUT
output_subdir <- trimws(Sys.getenv("GFORM_PERIOD_MULTI_OUTPUT_SUBDIR", unset = ""))
if (nzchar(output_subdir)) {
  data_out_g <- file.path(GFORM_PERIOD_MULTI_DATA_OUT, output_subdir)
  message("Outputs período multi en subdirectorio: ", data_out_g)
}

dir_interventions <- file.path(GFORM_PERIOD_MULTI_DATA_OUT, "Interventions")
dir_population <- file.path(data_out_g, "PopulationEffects")
dir_models <- file.path(data_out_g, "Models")
dir_bootstrap <- file.path(data_out_g, "Bootstrap")
dir_summary <- file.path(data_out_g, "Summary_results")
dir_timing <- file.path(data_out_g, "Timing")

execution_mode <- tolower(Sys.getenv("GFORM_EXEC_MODE", "auto"))
if (identical(execution_mode, "auto")) {
  execution_mode <- if (gform_is_linux_server()) "server" else "local"
}

parallel_config <- gform_parallel_config()
options(gform.parallel = parallel_config)

run_parallel_cox <- execution_mode == "server"
run_parallel_bootstrap <- execution_mode == "server" &&
  gform_env_bool("GFORM_BOOTSTRAP_PARALLEL", default = FALSE)

sample_frac_env <- gform_env_num("GFORM_SAMPLE_FRAC", NA_real_)
sample_frac <- if (is.finite(sample_frac_env) && sample_frac_env > 0 && sample_frac_env < 1) {
  sample_frac_env
} else {
  NULL
}
sample_seed <- as.integer(gform_env_num("GFORM_SAMPLE_SEED", 2026L))

intervention_env <- gform_env_int_vec("GFORM_PERIOD_MULTI_INTERVENTIONS", NULL)
if (is.null(intervention_env)) {
  intervention_env <- gform_env_int_vec("GFORM_INTERVENTIONS", NULL)
}
intervention_numbers <- if (!is.null(intervention_env)) {
  intervention_env
} else {
  seq_len(length(GFORM_PERIOD_MULTI_INTERVENTION_ORDER))
}

max_batch_hours <- as.numeric(gform_env_num(
  "GFORM_MAX_BATCH_HOURS",
  if (execution_mode == "server") 168 else 12
))

run_bootstrap <- gform_env_bool("GFORM_RUN_BOOTSTRAP", default = FALSE)
boot_iter_env <- gform_env_num("GFORM_BOOT_ITER", NA_real_)
boot_iter <- if (isTRUE(run_bootstrap)) {
  if (is.finite(boot_iter_env)) as.integer(boot_iter_env) else 250L
} else {
  0L
}
boot_seed <- as.integer(gform_env_num("GFORM_BOOT_SEED", GFORM_PERIOD_MULTI_DEFAULTS$boot_seed))
skip_completed <- gform_env_bool("GFORM_SKIP_COMPLETED", default = TRUE)

max_follow_up <- GFORM_PERIOD_MULTI_DEFAULTS$max_follow_up
risk_weeks <- GFORM_PERIOD_MULTI_DEFAULTS$risk_weeks
risk_entry_week <- GFORM_PERIOD_MULTI_DEFAULTS$risk_entry_week
follow_up_weeks <- GFORM_PERIOD_MULTI_DEFAULTS$follow_up_weeks
population_week <- GFORM_PERIOD_MULTI_DEFAULTS$population_week
control_vars <- GFORM_PERIOD_CONTROL_VARS
dependent_var <- GFORM_PERIOD_MULTI_DEFAULTS$dependent_var
baseline_scenario <- GFORM_PERIOD_MULTI_DEFAULTS$baseline_scenario

run_one_period_multi_intervention <- function(intervention_number, shared_data) {
  timing_log <- gform_timing_log_init()

  registry_key <- resolve_gform_period_multi_intervention(intervention_number)
  spec <- GFORM_PERIOD_MULTI_INTERVENTION_REGISTRY[[registry_key]]
  application_scope <- spec$application_scope
  exposure_vars <- spec$exposure_vars

  message("\n", strrep("=", 72))
  message(
    "Intervención ", intervention_number, ": ", spec$intervention_id,
    " [", application_scope, "] focal=", gform_period_pollutant_label(spec$pollutant)
  )
  message(spec$description)
  message("Cox compartido: ", spec$cox_formula)
  message(strrep("=", 72))

  data_base <- shared_data$data_base
  person_weeks <- shared_data$person_weeks
  total_births <- shared_data$total_births
  n_total_original <- shared_data$n_original
  exposure_natural <- get_period_multi_exposure_natural(
    shared_data$births_df,
    application_scope,
    shared_data$exposure_cache
  )

  block <- gform_time_block(
    paste0(
      "Cox curso natural multi — ", application_scope,
      " (cache: ", dir_models, ")"
    ),
    load_or_fit_period_multi_natural_course_models(
      application_scope = application_scope,
      exposure_vars = exposure_vars,
      dir_models = dir_models,
      data_base = data_base,
      exposure_natural = exposure_natural,
      risk_weeks_vec = risk_weeks,
      control_vars = control_vars,
      dependent_var = dependent_var,
      risk_entry_week = risk_entry_week,
      sample_frac = sample_frac,
      parallel = run_parallel_cox
    )
  )
  timing_log <- gform_timing_log_add(timing_log, block$timing)
  fit_nat <- block$result

  model_store <- fit_nat$model_store
  cox_frame_natural <- fit_nat$cox_frame_natural
  rm(fit_nat)
  gc(verbose = FALSE)

  intervention_path <- gform_period_multi_intervention_rds_path(
    spec$intervention_id,
    dir_interventions
  )
  if (!file.exists(intervention_path)) {
    stop("RDS no encontrado (ejecutar 16.1 primero): ", intervention_path)
  }

  population_path <- file.path(dir_population, paste0(spec$output_stub, "_population_effects.rds"))

  res <- run_period_gform_intervention(
    intervention_spec = spec,
    intervention_path = intervention_path,
    data_base = data_base,
    model_store = model_store,
    person_weeks = person_weeks,
    risk_weeks_vec = risk_weeks,
    control_vars = control_vars,
    total_births = total_births,
    cox_frame_natural = cox_frame_natural,
    exposure_natural = exposure_natural,
    dependent_var = dependent_var,
    risk_entry_week = risk_entry_week,
    follow_up_weeks = follow_up_weeks,
    boot_iter = boot_iter,
    boot_seed = boot_seed,
    target_week = population_week,
    run_bootstrap = run_bootstrap,
    parallel_bootstrap = run_parallel_bootstrap,
    dir_bootstrap = dir_bootstrap,
    bootstrap_resume = TRUE,
    minimal_output = TRUE,
    dir_population = dir_population,
    sample_frac = sample_frac
  )
  if (!is.null(res$timing)) {
    timing_log <- gform_timing_log_merge(timing_log, res$timing)
  }

  excel_path <- file.path(dir_summary, paste0(spec$output_stub, "_point_estimates.xlsx"))
  boot_paths <- gform_bootstrap_paths(spec$output_stub, dir_bootstrap)

  metadata <- list(
    intervention_id = spec$intervention_id,
    pollutant = spec$pollutant,
    application_scope = spec$application_scope,
    exposure_model = spec$exposure_model,
    cox_formula = spec$cox_formula,
    exposure_vars = exposure_vars,
    target_exposure_vars = spec$target_exposure_vars,
    description = spec$description,
    adjustment = "Adjusted",
    multicontaminant = TRUE,
    execution_mode = execution_mode,
    sample_frac = sample_frac,
    n_births = total_births,
    n_original = n_total_original,
    risk_weeks = risk_weeks,
    risk_entry_week = risk_entry_week,
    follow_up_weeks = follow_up_weeks,
    population_week = population_week,
    boot_iter = boot_iter,
    model_type = "coxph_period_multi",
    run_parallel_cox = run_parallel_cox,
    run_parallel_bootstrap = run_parallel_bootstrap,
    parallel_config = parallel_config,
    bootstrap_population_csv = boot_paths$population,
    run_time = Sys.time()
  )

  block <- gform_time_block("Guardar resultados (poblacional + Excel)", {
    save_period_results(
      population_effects = res$population_effects,
      population_path = population_path,
      metadata = metadata
    )
    save_period_excel(
      population_effects = res$population_effects,
      excel_path = excel_path,
      metadata = metadata
    )
  })
  timing_log <- gform_timing_log_add(timing_log, block$timing)

  nat_row <- res$population_effects |> dplyr::filter(.data$scenario == baseline_scenario)
  int_row <- res$population_effects |> dplyr::filter(.data$scenario == "intervention")

  rm(res, model_store, cox_frame_natural)
  gc(verbose = FALSE)

  timing_log$intervention_number <- intervention_number
  timing_log$intervention_id <- spec$intervention_id
  timing_log$application_scope <- application_scope
  timing_log$finished_at <- Sys.time()
  timing_log$wall_sec <- as.numeric(
    difftime(timing_log$finished_at, timing_log$started_at, units = "secs")
  )
  timing_log$total_sec <- gform_timing_total_sec(timing_log)

  gform_print_timing_summary(
    timing_log,
    paste0("Intervención ", intervention_number, ": ", spec$intervention_id)
  )

  timing_path <- file.path(
    dir_timing,
    paste0(spec$output_stub, "_timing_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds")
  )
  saveRDS(timing_log, timing_path)

  if (file.exists(gform_point_checkpoint_path(spec$output_stub, dir_bootstrap))) {
    file.remove(gform_point_checkpoint_path(spec$output_stub, dir_bootstrap))
  }

  list(
    intervention_number = intervention_number,
    intervention_id = spec$intervention_id,
    application_scope = application_scope,
    description = spec$description,
    timing_log = timing_log,
    timing_path = timing_path,
    prevalence_natural = nat_row$prevalence,
    prevalence_intervention = int_row$prevalence,
    risk_difference = int_row$risk_difference,
    output_files = c(population_path, excel_path, boot_paths$population)
  )
}

## ===== Lote de intervenciones =====
batch_start <- Sys.time()
batch_deadline <- batch_start + max_batch_hours * 3600
log_path <- file.path(
  dir_timing,
  paste0(execution_mode, "_run_", format(batch_start, "%Y%m%d_%H%M%S"), ".log")
)
dir.create(dir_timing, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_population, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_models, recursive = TRUE, showWarnings = FALSE)
if (run_bootstrap && boot_iter > 0L) {
  dir.create(dir_bootstrap, recursive = TRUE, showWarnings = FALSE)
}
dir.create(dir_summary, recursive = TRUE, showWarnings = FALSE)
writeLines(log_path, file.path(dir_timing, "batch_run.logpath"))

on.exit(gform_finalize_run(), add = TRUE)

sink(log_path, split = TRUE)
on.exit(sink(), add = TRUE)

message("\n=== G-Formula Período Multi — corrida (modo: ", execution_mode, ") ===")
print_gform_period_multi_cox_specs()
message("Inicio: ", batch_start)
message("Límite: ", max_batch_hours, " h (hasta ", batch_deadline, ")")
if (!is.null(sample_frac)) {
  message("Submuestra: ", round(100 * sample_frac, 2), "% (seed ", sample_seed, ")")
}
message(
  "Intervenciones: ",
  paste(intervention_numbers, collapse = ", "),
  " | bootstrap: ",
  if (run_bootstrap) paste0(boot_iter, " iter (",
    if (run_parallel_bootstrap) "paralelo" else "secuencial", ")") else "no"
)
message(
  "CPUs: ", parallel_config$n_cores,
  " | RAM: ", round(parallel_config$ram_gb, 1), " GiB",
  " | fork: ", parallel_config$use_fork
)
if (execution_mode == "server") {
  gform_setup_parallel(task = "default", config = parallel_config)
  message("Paralelo activo en ajuste Cox por semana de riesgo.")
} else {
  message("Modo local: intervenciones secuenciales, sin paralelo interno.")
}
print_gform_period_multi_intervention_menu()

batch_timing <- gform_timing_log_init()

block <- gform_time_block("Cargar cohorte período multi", {
  births_df <- prepare_period_births_cohort(data_inp)
  if (!is.null(sample_frac)) {
    set.seed(sample_seed)
    ids_keep <- sample(
      unique(births_df$id),
      size = floor(sample_frac * length(unique(births_df$id)))
    )
    births_df <- births_df[births_df$id %in% ids_keep, , drop = FALSE]
  }
  births_df
})
batch_timing <- gform_timing_log_add(batch_timing, block$timing)
births_df <- block$result
n_total_original <- nrow(births_df)

block <- gform_time_block("Preparar person_weeks", {
  data_base <- births_df |>
    dplyr::select(dplyr::any_of(c("id", "weeks", dependent_var, control_vars)))
  person_weeks <- expand_person_weeks(
    births_df = data_base,
    risk_weeks = risk_weeks,
    max_fu = max_follow_up,
    dependent_var = dependent_var
  )
  total_births <- nrow(data_base)
  list(
    data_base = data_base,
    person_weeks = person_weeks,
    total_births = total_births
  )
})
batch_timing <- gform_timing_log_add(batch_timing, block$timing)

block2 <- gform_time_block("Inicializar cache de exposiciones naturales multi", {
  list()
})
batch_timing <- gform_timing_log_add(batch_timing, block2$timing)

shared_data <- list(
  data_base = block$result$data_base,
  person_weeks = block$result$person_weeks,
  total_births = block$result$total_births,
  n_original = n_total_original,
  births_df = births_df,
  exposure_cache = block2$result
)

message("Nacimientos (análisis): ", shared_data$total_births)
if (!is.null(sample_frac)) {
  message("Nacimientos (cohorte original): ", n_total_original)
}

rm(births_df)
gc(verbose = FALSE)

gform_print_timing_summary(batch_timing, "Preparación del lote período multi")

batch_results <- list()

for (int_num in intervention_numbers) {
  if (Sys.time() >= batch_deadline) {
    message("\nLímite de ", max_batch_hours, " h alcanzado; deteniendo lote.")
    break
  }

  spec <- GFORM_PERIOD_MULTI_INTERVENTION_REGISTRY[[
    resolve_gform_period_multi_intervention(int_num)
  ]]
  if (skip_completed && gform_period_multi_intervention_is_complete(
    output_stub = spec$output_stub,
    expected_application_scope = spec$application_scope,
    dir_population = dir_population,
    dir_summary = dir_summary,
    dir_bootstrap = dir_bootstrap,
    boot_iter = boot_iter,
    boot_seed = boot_seed,
    total_births = shared_data$total_births,
    sample_frac = sample_frac,
    run_bootstrap = run_bootstrap
  )) {
    message("\nIntervención ", int_num, " (", spec$intervention_id, ") ya completa — omitida.")
    next
  }

  result <- tryCatch(
    run_one_period_multi_intervention(int_num, shared_data),
    error = function(e) {
      call_txt <- if (!is.null(e$call)) {
        paste(deparse(e$call, width.cutoff = 120L), collapse = "\n")
      } else {
        ""
      }
      message("\n*** ERROR intervención ", int_num, ": ", conditionMessage(e))
      if (nzchar(call_txt)) message("En: ", call_txt)
      tb <- utils::capture.output(base::traceback(max.lines = 25L))
      if (length(tb)) message(paste(tb, collapse = "\n"))
      gc(verbose = FALSE)
      NULL
    }
  )

  if (is.null(result)) {
    message("Intervención ", int_num, " falló; continuando con la siguiente.")
    next
  }

  batch_results[[as.character(int_num)]] <- result
  gc(verbose = TRUE)
}

batch_end <- Sys.time()
batch_log <- list(
  started_at = batch_start,
  finished_at = batch_end,
  total_sec = as.numeric(difftime(batch_end, batch_start, units = "secs")),
  execution_mode = execution_mode,
  parallel_config = parallel_config,
  max_batch_hours = max_batch_hours,
  intervention_numbers = intervention_numbers,
  completed = batch_results,
  log_path = log_path
)
saveRDS(batch_log, file.path(
  dir_timing,
  paste0(execution_mode, "_batch_", format(batch_start, "%Y%m%d_%H%M%S"), ".rds")
))

message("\n=== Resumen del lote período multi ===")
if (length(batch_results) == 0L) {
  message("Ninguna intervención completada en esta sesión.")
} else {
  for (nm in names(batch_results)) {
    r <- batch_results[[nm]]
    message(sprintf(
      "  %2s. %-28s [%s] | RD=% .6f | %.1f min",
      r$intervention_number, r$intervention_id, r$application_scope,
      r$risk_difference, r$timing_log$wall_sec / 60
    ))
  }
}

message("Log: ", log_path)
message(
  "Tiempo total del lote: ",
  gform_format_duration(as.numeric(difftime(batch_end, batch_start, units = "secs"))),
  " (inicio ", gform_format_timestamp(batch_start),
  " → fin ", gform_format_timestamp(batch_end), ")"
)

beepr::beep(8)
