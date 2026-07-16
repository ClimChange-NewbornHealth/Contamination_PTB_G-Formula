# 10.2 G-Formula — corrida de intervenciones (local secuencial / servidor paralelo) ----
#
# Requisito previo: 10.0 G-Form_build_interventions.R (una vez)
#
# Uso (desde la raíz del proyecto):
#   bash "00_Code/run_gform_server.sh"     # servidor Linux (paralelo, primer plano)
#   Rscript "00_Code/10.2 G-Form_models.R" # respeta GFORM_EXEC_MODE

rm(list = ls())
gc(verbose = FALSE)

source("00_Code/0.2 Packages_gform.R")
source("00_Code/0.1 Settings.R")
source("00_Code/0.3 Functions.R")
source("00_Code/10.1 G-Form_functions.R")

data_inp <- "01_Data/Output/"
data_out_g <- "02_Output/G-Form/"
output_subdir <- trimws(Sys.getenv("GFORM_OUTPUT_SUBDIR", unset = ""))
if (nzchar(output_subdir)) {
  data_out_g <- file.path("02_Output/G-Form", output_subdir)
  message("Outputs en subdirectorio: ", data_out_g)
}
dir_interventions <- file.path("02_Output/G-Form", "Interventions")
dir_weekly <- file.path(data_out_g, "WeeklyEffects")
dir_population <- file.path(data_out_g, "PopulationEffects")
dir_other <- file.path(data_out_g, "Other")
dir_models <- file.path(data_out_g, "Models")
dir_bootstrap <- file.path(data_out_g, "Bootstrap")
dir_heatmap <- file.path(data_out_g, "Heatmap")
dir_summary <- file.path(data_out_g, "Summary_results")
dir_timing <- file.path(data_out_g, "Timing")

## ===== Modo de ejecución =====
# GFORM_EXEC_MODE: "server" | "local" | "auto" (default)
execution_mode <- tolower(Sys.getenv("GFORM_EXEC_MODE", "auto"))
if (identical(execution_mode, "auto")) {
  execution_mode <- if (gform_is_linux_server()) "server" else "local"
}

parallel_config <- gform_parallel_config()
options(gform.parallel = parallel_config)

run_parallel_cox <- execution_mode == "server"
run_parallel_bootstrap <- execution_mode == "server" &&
  gform_env_bool("GFORM_BOOTSTRAP_PARALLEL", default = FALSE)
run_parallel_singleweek_heatmap <- execution_mode == "server" &&
  gform_env_bool("GFORM_HEATMAP_PARALLEL", default = FALSE)

## ===== Configuración (variables de entorno opcionales) =====
sample_frac_env <- gform_env_num("GFORM_SAMPLE_FRAC", NA_real_)
sample_frac <- if (is.finite(sample_frac_env) && sample_frac_env > 0 && sample_frac_env < 1) {
  sample_frac_env
} else {
  NULL
}
sample_seed <- as.integer(gform_env_num("GFORM_SAMPLE_SEED", 2026L))

intervention_env <- gform_env_int_vec("GFORM_INTERVENTIONS", NULL)
intervention_numbers <- if (!is.null(intervention_env)) {
  intervention_env
} else {
  seq_len(length(GFORM_INTERVENTION_ORDER))
}

max_batch_hours <- as.numeric(gform_env_num(
  "GFORM_MAX_BATCH_HOURS",
  if (execution_mode == "server") 168 else 12
))

run_bootstrap <- gform_env_bool("GFORM_RUN_BOOTSTRAP", default = TRUE)
boot_iter_env <- gform_env_num("GFORM_BOOT_ITER", NA_real_)
boot_iter <- if (is.finite(boot_iter_env)) {
  as.integer(boot_iter_env)
} else {
  GFORM_DEFAULTS$boot_iter
}
boot_seed <- as.integer(gform_env_num("GFORM_BOOT_SEED", GFORM_DEFAULTS$boot_seed))
run_singleweek_heatmap <- gform_env_bool("GFORM_RUN_HEATMAP", default = TRUE)
skip_completed <- gform_env_bool("GFORM_SKIP_COMPLETED", default = TRUE)
heatmap_only <- gform_env_bool("GFORM_HEATMAP_ONLY", default = FALSE)
if (heatmap_only) {
  message(
    "GFORM_HEATMAP_ONLY activo: solo recalcula heatmaps cuando ya existen RDS puntuales."
  )
}

max_follow_up <- GFORM_DEFAULTS$max_follow_up
weeks_exposure <- GFORM_DEFAULTS$weeks_exposure
lag_weeks <- GFORM_DEFAULTS$lag_weeks
risk_weeks <- GFORM_DEFAULTS$risk_weeks
risk_entry_week <- GFORM_DEFAULTS$risk_entry_week
follow_up_weeks <- GFORM_DEFAULTS$follow_up_weeks
population_week <- GFORM_DEFAULTS$population_week
control_vars <- GFORM_DEFAULTS$control_vars
dependent_var <- GFORM_DEFAULTS$dependent_var
baseline_scenario <- GFORM_DEFAULTS$baseline_scenario

load_pollutant_long <- function(pollutant, ids = NULL) {
  load(paste0(data_inp, "births_2010_2020_exposure_weeks_lagged.RData"))
  out <- data_long |>
    dplyr::select(id, week_gest_num, tad, dplyr::all_of(pollutant))
  if (!is.null(ids)) {
    out <- out[out$id %in% ids, , drop = FALSE]
  }
  rm(data_long)
  gc(verbose = FALSE)
  out
}

run_one_gform_intervention <- function(intervention_number, shared_data) {
  timing_log <- gform_timing_log_init()

  registry_key <- resolve_gform_intervention(intervention_number)
  spec <- GFORM_INTERVENTION_REGISTRY[[registry_key]]
  pollutant <- spec$pollutant

  message("\n", strrep("=", 72))
  message("Intervención ", intervention_number, ": ", spec$intervention_id)
  message(spec$description)
  message(strrep("=", 72))

  data_base <- shared_data$data_base
  wide_tad_obs <- shared_data$wide_tad_obs
  person_weeks <- shared_data$person_weeks
  total_births <- shared_data$total_births
  n_total_original <- shared_data$n_original

  block <- gform_time_block(
    paste0("Cargar exposición semanal — ", pollutant),
    load_pollutant_long(pollutant, ids = data_base$id)
  )
  timing_log <- gform_timing_log_add(timing_log, block$timing)
  data_long_weekly <- block$result

  block <- gform_time_block("Preparar matrices de exposición (natural)", {
    raw_pollutant <- build_wide_raw_exposure(
      data_long_weekly, pollutant, weeks_keep = weeks_exposure
    )
    wide_exposicion_natural <- build_exposicion_wide_from_raw(
      raw_wide = raw_pollutant,
      pollutant = pollutant,
      intervention = list(type = "none"),
      weeks_keep = weeks_exposure,
      lag_weeks = lag_weeks
    )
    list(
      raw_pollutant = raw_pollutant,
      wide_exposicion_natural = wide_exposicion_natural
    )
  })
  timing_log <- gform_timing_log_add(timing_log, block$timing)
  raw_pollutant <- block$result$raw_pollutant
  wide_exposicion_natural <- block$result$wide_exposicion_natural

  block <- gform_time_block(
    paste0("Cox bajo curso natural — ", pollutant, " (cache: ", dir_models, ")"),
    load_or_fit_natural_course_models(
      pollutant = pollutant,
      dir_models = dir_models,
      data_base = data_base,
      wide_exposicion_natural = wide_exposicion_natural,
      wide_tad_obs = wide_tad_obs,
      risk_weeks_vec = risk_weeks,
      control_vars = control_vars,
      lag_weeks = lag_weeks,
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
  rm(fit_nat, wide_exposicion_natural)
  gc(verbose = FALSE)

  intervention_path <- file.path(dir_interventions, paste0(spec$intervention_id, ".rds"))
  if (!file.exists(intervention_path)) {
    stop("RDS no encontrado (ejecutar 10.0 primero): ", intervention_path)
  }

  output_paths <- gform_output_paths(spec$output_stub, dir_other)
  weekly_path <- file.path(dir_weekly, paste0(spec$output_stub, "_weekly_effects.rds"))
  population_path <- file.path(dir_population, paste0(spec$output_stub, "_population_effects.rds"))
  use_heatmap_only <- isTRUE(heatmap_only) &&
    file.exists(weekly_path) &&
    file.exists(population_path) &&
    file.exists(output_paths$cumulative_risk_curves) &&
    !gform_heatmap_is_complete(spec$output_stub, dir_heatmap, dir_other)
  if (use_heatmap_only) {
    message("Modo heatmap_only para ", spec$output_stub, " (RDS previos + heatmap pendiente).")
  }

  res <- run_gform_intervention(
    intervention_spec = spec,
    intervention_path = intervention_path,
    data_base = data_base,
    wide_tad_obs = wide_tad_obs,
    model_store = model_store,
    person_weeks = person_weeks,
    risk_weeks_vec = risk_weeks,
    control_vars = control_vars,
    total_births = total_births,
    cox_frame_natural = cox_frame_natural,
    raw_wide_pollutant = raw_pollutant,
    dependent_var = dependent_var,
    risk_entry_week = risk_entry_week,
    follow_up_weeks = follow_up_weeks,
    boot_iter = boot_iter,
    boot_seed = boot_seed,
    target_week = population_week,
    run_bootstrap = if (use_heatmap_only) FALSE else run_bootstrap,
    run_singleweek_heatmap = run_singleweek_heatmap,
    parallel_singleweek_heatmap = run_parallel_singleweek_heatmap,
    parallel_bootstrap = run_parallel_bootstrap,
    dir_bootstrap = dir_bootstrap,
    bootstrap_resume = TRUE,
    dir_heatmap = dir_heatmap,
    heatmap_resume = TRUE,
    heatmap_only = use_heatmap_only,
    dir_weekly = dir_weekly,
    dir_population = dir_population,
    dir_other = dir_other,
    sample_frac = sample_frac
  )
  if (!is.null(res$timing)) {
    timing_log <- gform_timing_log_merge(timing_log, res$timing)
  }

  output_paths <- gform_output_paths(spec$output_stub, dir_other)
  excel_path <- file.path(dir_summary, paste0(spec$output_stub, "_point_estimates.xlsx"))
  boot_paths <- gform_bootstrap_paths(spec$output_stub, dir_bootstrap)

  metadata <- list(
    intervention_id = spec$intervention_id,
    pollutant = spec$pollutant,
    description = spec$description,
    execution_mode = execution_mode,
    sample_frac = sample_frac,
    n_births = total_births,
    n_original = n_total_original,
    risk_weeks = risk_weeks,
    risk_entry_week = risk_entry_week,
    follow_up_weeks = follow_up_weeks,
    population_week = population_week,
    boot_iter = if (run_bootstrap) boot_iter else 0L,
    model_type = "coxph",
    run_parallel_cox = run_parallel_cox,
    run_parallel_bootstrap = run_parallel_bootstrap,
    run_singleweek_heatmap = run_singleweek_heatmap,
    run_parallel_singleweek_heatmap = run_parallel_singleweek_heatmap,
    parallel_config = parallel_config,
    bootstrap_weekly_csv = boot_paths$weekly,
    bootstrap_population_csv = boot_paths$population,
    run_time = Sys.time()
  )

  block <- gform_time_block("Guardar resultados (RDS + Excel)", {
    save_results(
      weekly_effects = res$weekly_effects,
      population_effects = res$population_effects,
      weekly_path = weekly_path,
      population_path = population_path,
      cumulative_risk_curves = res$cumulative_risk_curves,
      cumulative_risk_curves_path = output_paths$cumulative_risk_curves,
      singleweek_intervention_heatmap = res$singleweek_intervention_heatmap,
      singleweek_intervention_heatmap_path = output_paths$singleweek_heatmap,
      weekly_boot = NULL,
      population_boot = NULL,
      metadata = metadata
    )
    save_gform_excel(results = res, excel_path = excel_path, intervention_id = spec$intervention_id)
  })
  timing_log <- gform_timing_log_add(timing_log, block$timing)

  nat_row <- res$population_effects |> dplyr::filter(.data$scenario == baseline_scenario)
  int_row <- res$population_effects |> dplyr::filter(.data$scenario == "intervention")

  rm(res, model_store, cox_frame_natural, raw_pollutant, data_long_weekly)
  gc(verbose = FALSE)

  timing_log$intervention_number <- intervention_number
  timing_log$intervention_id <- spec$intervention_id
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
  heatmap_paths <- gform_heatmap_paths(spec$output_stub, dir_heatmap)
  if (file.exists(heatmap_paths$checkpoint)) {
    file.remove(heatmap_paths$checkpoint)
  }

  list(
    intervention_number = intervention_number,
    intervention_id = spec$intervention_id,
    description = spec$description,
    timing_log = timing_log,
    timing_path = timing_path,
    prevalence_natural = nat_row$prevalence,
    prevalence_intervention = int_row$prevalence,
    risk_difference = int_row$risk_difference,
    output_files = c(
      weekly_path, population_path, excel_path,
      output_paths$cumulative_risk_curves, output_paths$singleweek_heatmap,
      boot_paths$weekly, boot_paths$population
    )
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
writeLines(log_path, file.path(dir_timing, "batch_run.logpath"))

on.exit(gform_finalize_run(), add = TRUE)

sink(log_path, split = TRUE)
on.exit(sink(), add = TRUE)

message("\n=== G-Formula corrida (modo: ", execution_mode, ") ===")
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
    if (run_parallel_bootstrap) "paralelo" else "secuencial", ")") else "no",
  if (heatmap_only) " | heatmap_only: auto por intervención" else ""
)
message(
  "CPUs: ", parallel_config$n_cores,
  " | RAM: ", round(parallel_config$ram_gb, 1), " GiB",
  " | fork: ", parallel_config$use_fork
)
message(
  "Workers Cox/bootstrap/heatmap: ",
  parallel_config$n_workers_cox, " / ",
  parallel_config$n_workers_bootstrap, " / ",
  parallel_config$n_workers_heatmap,
  " (bootstrap lote: ", parallel_config$bootstrap_batch_size,
  " | heatmap lote: ", parallel_config$heatmap_batch_size, ")"
)
message(
  "Bootstrap RAM: reserva padre ", parallel_config$bootstrap_parent_reserve_gb,
  " GiB | ~", parallel_config$bootstrap_ram_per_worker_gb, " GiB/worker"
)
message(
  "Heatmap RAM: reserva padre ", parallel_config$heatmap_parent_reserve_gb,
  " GiB | ~", parallel_config$heatmap_ram_per_worker_gb, " GiB/worker",
  " | paralelo: ", run_parallel_singleweek_heatmap
)
if (execution_mode == "server") {
  gform_setup_parallel(task = "default", config = parallel_config)
  message("Paralelo activo en Cox, bootstrap por lotes y mapa calor.")
} else {
  message("Modo local: intervenciones secuenciales, sin paralelo interno.")
}
print_gform_intervention_menu()

pollutants_needed <- unique(vapply(intervention_numbers, function(n) {
  GFORM_INTERVENTION_REGISTRY[[resolve_gform_intervention(n)]]$pollutant
}, character(1)))

batch_timing <- gform_timing_log_init()

block <- gform_time_block("Cargar cohorte base", {
  births <- rio::import(paste0(data_inp, "births_exposure_period_metrics_full_d30_d4_tri.RData")) |>
    dplyr::select(id, ndvi_full) |>
    dplyr::distinct(id, .keep_all = TRUE)
  load(paste0(data_inp, "births_2010_2020_exposure_weeks_lagged.RData"))

  data_long <- data_long |>
    dplyr::select(dplyr::any_of(unique(c(
      "id", "weeks", "week_gest_num", "tad", dependent_var, control_vars
    )))) |>
    dplyr::left_join(births, by = "id") |>
    dplyr::mutate(
      month_week1 = factor(month_week1),
      year_week1 = factor(year_week1),
      covid = factor(covid)
    )
  rm(births)
  gc(verbose = FALSE)
  data_long
})
batch_timing <- gform_timing_log_add(batch_timing, block$timing)
data_long <- block$result
n_total_original <- dplyr::n_distinct(data_long$id)

if (!is.null(sample_frac)) {
  block <- gform_time_block(
    paste0("Submuestra ", round(100 * sample_frac, 2), "%"),
    gform_subsample_data_long(data_long, sample_frac, sample_seed)
  )
  batch_timing <- gform_timing_log_add(batch_timing, block$timing)
  data_long <- block$result
}

message("Nacimientos (análisis): ", dplyr::n_distinct(data_long$id))
if (!is.null(sample_frac)) {
  message("Nacimientos (cohorte original): ", n_total_original)
}
message("Contaminantes en lote: ", paste(pollutants_needed, collapse = ", "))

block <- gform_time_block("Preparar person_weeks y TAD semanal", {
  data_base <- data_long |>
    dplyr::select(dplyr::any_of(c("id", "weeks", dependent_var, control_vars))) |>
    dplyr::distinct(id, .keep_all = TRUE)
  wide_tad_obs <- build_wide_weekly_var(data_long, "tad", weeks_keep = weeks_exposure)
  wide_tad_obs <- data.table::as.data.table(wide_tad_obs)
  person_weeks <- expand_person_weeks(
    births_df = data_base,
    risk_weeks = risk_weeks,
    max_fu = max_follow_up,
    dependent_var = dependent_var
  )
  total_births <- nrow(data_base)
  rm(data_long)
  gc(verbose = FALSE)
  list(
    data_base = data_base,
    wide_tad_obs = wide_tad_obs,
    person_weeks = person_weeks,
    total_births = total_births
  )
})
batch_timing <- gform_timing_log_add(batch_timing, block$timing)

shared_data <- list(
  data_base = block$result$data_base,
  wide_tad_obs = block$result$wide_tad_obs,
  person_weeks = block$result$person_weeks,
  total_births = block$result$total_births,
  n_original = n_total_original
)

gform_print_timing_summary(batch_timing, "Preparación del lote")

batch_results <- list()

for (int_num in intervention_numbers) {
  if (Sys.time() >= batch_deadline) {
    message("\nLímite de ", max_batch_hours, " h alcanzado; deteniendo lote.")
    break
  }

  spec <- GFORM_INTERVENTION_REGISTRY[[resolve_gform_intervention(int_num)]]
  if (skip_completed && gform_intervention_is_complete(
    output_stub = spec$output_stub,
    dir_weekly = dir_weekly,
    dir_population = dir_population,
    dir_summary = dir_summary,
    dir_bootstrap = dir_bootstrap,
    boot_iter = boot_iter,
    boot_seed = boot_seed,
    total_births = shared_data$total_births,
    sample_frac = sample_frac,
    run_heatmap = run_singleweek_heatmap,
    dir_heatmap = dir_heatmap,
    dir_other = dir_other
  )) {
    message("\nIntervención ", int_num, " (", spec$intervention_id, ") ya completa — omitida.")
    next
  }

  result <- tryCatch(
    run_one_gform_intervention(int_num, shared_data),
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

message("\n=== Resumen del lote ===")
if (length(batch_results) == 0L) {
  message("Ninguna intervención completada en esta sesión.")
} else {
  for (nm in names(batch_results)) {
    r <- batch_results[[nm]]
    message(sprintf(
      "  %2s. %-18s | RD=% .6f | %.1f min (reloj)",
      r$intervention_number, r$intervention_id,
      r$risk_difference, r$timing_log$wall_sec / 60
    ))
  }
}
batch_timing$intervention_totals <- vapply(
  batch_results,
  function(r) r$timing_log$wall_sec,
  numeric(1)
)
batch_timing$finished_at <- batch_end
batch_timing$wall_sec <- as.numeric(
  difftime(batch_end, batch_start, units = "secs")
)
batch_timing$total_sec <- batch_timing$wall_sec
gform_print_timing_summary(
  list(
    phases = lapply(names(batch_results), function(nm) {
      r <- batch_results[[nm]]
      list(
        label = sprintf("%2s. %s", r$intervention_number, r$intervention_id),
        sec = r$timing_log$wall_sec,
        start_str = gform_format_timestamp(r$timing_log$started_at),
        end_str = gform_format_timestamp(r$timing_log$finished_at)
      )
    }),
    started_at = batch_start,
    finished_at = batch_end,
    wall_sec = batch_timing$wall_sec
  ),
  "Tiempo total por intervención completada"
)
message("Log: ", log_path)
message(
  "Tiempo total del lote: ",
  gform_format_duration(batch_timing$wall_sec),
  " (inicio ", gform_format_timestamp(batch_start),
  " → fin ", gform_format_timestamp(batch_end), ")"
)

beepr::beep(8)
