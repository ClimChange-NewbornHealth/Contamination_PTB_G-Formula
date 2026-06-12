# 10.2 G-Formula — análisis principal (Cox + g-computación) ----
#
# Pipeline armonizado:
#   1. Carga datos y submuestra opcional
#   2. Ajuste Cox bajo curso natural (referencia 9.0, tstart = 28)
#   3. Iteración sobre intervenciones (RDS en Interventions/)
#   4. Salidas: efectos semanales, poblacionales, tablas Fig. 3/4 (Other), Summary_results
#
# Requisito previo: ejecutar 10.0 G-Form_build_interventions.R (una vez)
#
# Para muestra completa: sample_frac <- NULL, interventions_to_run <- NULL
# Para prueba rápida:   sample_frac <- 0.01, interventions_to_run <- NULL (solo la 1.ª)
#
# Uso (desde la raíz del proyecto):
#   Rscript "00_Code/10.2 G-Form_models.R"

rm(list = ls())

## Settings ----
source("00_Code/0.1 Settings.R")
source("00_Code/0.2 Packages.R")
source("00_Code/0.3 Functions.R")
source("00_Code/10.1 G-Form_functions.R")

data_inp <- "01_Data/Output/"
data_out_g <- "02_Output/G-Form/"
dir_interventions <- file.path(data_out_g, "Interventions")
dir_weekly <- file.path(data_out_g, "WeeklyEffects")
dir_population <- file.path(data_out_g, "PopulationEffects")
dir_other <- file.path(data_out_g, "Other")
dir_summary <- file.path(data_out_g, "Summary_results")
dir_timing <- file.path(data_out_g, "Timing")

## ===== Configuración de ejecución =====
# NULL = cohorte completa; 0.01 = 1% de la muestra original
sample_frac <- 0.01
sample_seed <- 2026L

# NULL = todas las intervenciones del registro; vector = subset por intervention_id
# En modo prueba se usa solo la primera intervención del registro
interventions_to_run <- NULL
test_first_intervention_only <- TRUE

run_bootstrap <- TRUE
boot_iter <- GFORM_DEFAULTS$boot_iter
boot_seed <- GFORM_DEFAULTS$boot_seed

# Figura 4 (heatmap intervención semana a semana): costoso; activar en prueba
run_figure4 <- TRUE

n_workers <- 10L  # detectCores() - 4 en máquina local (ajustar si es necesario)

## Parámetros del modelo ----
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

## Inicio ----
timing_log <- list()
run_start <- Sys.time()

message("=== G-Formula análisis principal (Cox) ===")
message("Workers: ", n_workers)
message("sample_frac: ", if (is.null(sample_frac)) "cohorte completa" else sample_frac)
message("risk_weeks: ", min(risk_weeks), "-", max(risk_weeks),
        " | tstart: ", risk_entry_week)

gform_setup_parallel(n_workers = n_workers)

## Carga de datos ----
message("\nCargando datos...")

births <- rio::import(paste0(data_inp, "births_exposure_period_metrics_full_d30_d4_tri.RData")) |>
  dplyr::select(id, ndvi_full) |>
  dplyr::distinct(id, .keep_all = TRUE)

load(paste0(data_inp, "births_2010_2020_exposure_weeks_lagged.RData"))

data_long <- data_long |>
  dplyr::left_join(births, by = "id") |>
  dplyr::mutate(
    month_week1 = factor(month_week1),
    year_week1 = factor(year_week1),
    covid = factor(covid)
  )

n_total_original <- dplyr::n_distinct(data_long$id)
message("Nacimientos cohorte original: ", n_total_original)

if (!is.null(sample_frac)) {
  set.seed(sample_seed)
  n_sample <- max(50L, floor(n_total_original * sample_frac))
  ids_s <- sample(unique(data_long$id), size = min(n_sample, n_total_original))
  data_long <- dplyr::filter(data_long, id %in% ids_s)
  message("Submuestra: n = ", length(ids_s),
          " (", round(100 * length(ids_s) / n_total_original, 2), "%)")
}

base_vars <- c("id", "weeks", dependent_var, control_vars)
data_base <- data_long |>
  dplyr::select(dplyr::any_of(base_vars)) |>
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

## Selección de intervenciones ----
registry_keys <- names(GFORM_INTERVENTION_REGISTRY)
if (!is.null(interventions_to_run)) {
  registry_keys <- interventions_to_run
  stopifnot(all(registry_keys %in% names(GFORM_INTERVENTION_REGISTRY)))
} else if (test_first_intervention_only && !is.null(sample_frac)) {
  registry_keys <- registry_keys[1L]
  message("Modo prueba: solo intervención '", registry_keys, "'")
}

## Iteración por intervención (agrupada por contaminante) ----
all_results <- list()

interventions_by_pollutant <- split(registry_keys, vapply(registry_keys, function(k) {
  GFORM_INTERVENTION_REGISTRY[[k]]$pollutant
}, character(1)))

for (pollutant in names(interventions_by_pollutant)) {
  keys_poll <- interventions_by_pollutant[[pollutant]]

  message("\n=== Contaminante: ", pollutant, " (", length(keys_poll), " intervenciones) ===")

  raw_pollutant_fit <- build_wide_raw_exposure(data_long, pollutant, weeks_keep = weeks_exposure)
  wide_exposicion_natural <- build_exposicion_wide_from_raw(
    raw_wide = raw_pollutant_fit,
    pollutant = pollutant,
    intervention = list(type = "none"),
    weeks_keep = weeks_exposure,
    lag_weeks = lag_weeks
  )

  message("Ajustando modelos Cox bajo curso natural...")
  tictoc::tic(paste0("fit_cox_", pollutant))
  model_store <- fit_natural_course_models(
    data_base = data_base,
    wide_exposicion_natural = wide_exposicion_natural,
    wide_tad_obs = wide_tad_obs,
    risk_weeks_vec = risk_weeks,
    control_vars = control_vars,
    lag_weeks = lag_weeks,
    dependent_var = dependent_var,
    risk_entry_week = risk_entry_week,
    parallel = TRUE
  )
  fit_toc <- tictoc::toc(quiet = TRUE)
  timing_log[[paste0("fit_cox_", pollutant)]] <- fit_toc$toc - fit_toc$tic

  n_models <- sum(!vapply(model_store, is.null, logical(1)))
  message("Modelos Cox ajustados: ", n_models)

  for (key in keys_poll) {
    spec <- GFORM_INTERVENTION_REGISTRY[[key]]
    intervention_path <- file.path(dir_interventions, paste0(spec$intervention_id, ".rds"))

    if (!file.exists(intervention_path)) {
      warning("RDS no encontrado (ejecutar 10.0 primero): ", intervention_path)
      next
    }

    message("\n--- Intervención: ", spec$description, " ---")

    tictoc::tic(paste0("run_", spec$output_stub))
    res <- run_gform_intervention(
      intervention_spec = spec,
      intervention_path = intervention_path,
      data_base = data_base,
      data_long = data_long,
      wide_tad_obs = wide_tad_obs,
      model_store = model_store,
      person_weeks = person_weeks,
      risk_weeks_vec = risk_weeks,
      control_vars = control_vars,
      total_births = total_births,
      dependent_var = dependent_var,
      risk_entry_week = risk_entry_week,
      follow_up_weeks = follow_up_weeks,
      boot_iter = boot_iter,
      boot_seed = boot_seed,
      target_week = population_week,
      run_bootstrap = run_bootstrap,
      run_figure4 = run_figure4,
      parallel = TRUE
    )
    run_toc <- tictoc::toc(quiet = TRUE)
    timing_log[[paste0("run_", spec$output_stub)]] <- run_toc$toc - run_toc$tic

    metadata <- list(
      intervention_id = spec$intervention_id,
      pollutant = spec$pollutant,
      description = spec$description,
      sample_frac = sample_frac,
      n_births = total_births,
      n_original = n_total_original,
      risk_weeks = risk_weeks,
      risk_entry_week = risk_entry_week,
      follow_up_weeks = follow_up_weeks,
      population_week = population_week,
      boot_iter = if (run_bootstrap) boot_iter else 0L,
      model_type = "coxph",
      run_time = Sys.time()
    )

    weekly_path <- file.path(dir_weekly, paste0(spec$output_stub, "_weekly_effects.rds"))
    population_path <- file.path(dir_population, paste0(spec$output_stub, "_population_effects.rds"))
    figure3_path <- file.path(dir_other, paste0(spec$output_stub, "_figure3.rds"))
    figure4_path <- file.path(dir_other, paste0(spec$output_stub, "_figure4.rds"))
    excel_path <- file.path(dir_summary, paste0(spec$output_stub, "_point_estimates.xlsx"))

    save_results(
      weekly_effects = res$weekly_effects,
      population_effects = res$population_effects,
      weekly_path = weekly_path,
      population_path = population_path,
      figure3 = res$figure3,
      figure3_path = figure3_path,
      figure4 = res$figure4,
      figure4_path = figure4_path,
      weekly_boot = if (run_bootstrap) res$bootstrap$weekly_boot else NULL,
      population_boot = if (run_bootstrap) res$bootstrap$population_boot else NULL,
      metadata = metadata
    )

    save_gform_excel(
      results = res,
      excel_path = excel_path,
      intervention_id = spec$intervention_id
    )

    all_results[[key]] <- res

    nat_row <- res$population_effects |> dplyr::filter(.data$scenario == baseline_scenario)
    int_row <- res$population_effects |> dplyr::filter(.data$scenario == "intervention")
    message("Prevalencia natural: ", round(nat_row$prevalence, 6),
            " | intervención: ", round(int_row$prevalence, 6),
            " | RD: ", round(int_row$risk_difference, 6))
  }
}

## Reporte de tiempos ----
run_end <- Sys.time()
timing_log$total_sec <- as.numeric(difftime(run_end, run_start, units = "secs"))
timing_log$n_births <- total_births
timing_log$n_original <- n_total_original
timing_log$sample_frac <- sample_frac
timing_log$n_workers <- n_workers
timing_log$interventions_run <- registry_keys

dir.create(dir_timing, recursive = TRUE, showWarnings = FALSE)
timing_path <- file.path(
  dir_timing,
  paste0("timing_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds")
)
saveRDS(timing_log, timing_path)

future::plan(future::sequential)

message("\n=== Tiempos de ejecución (segundos) ===")
for (nm in names(timing_log)) {
  if (nm %in% c("n_births", "n_original", "sample_frac", "n_workers", "interventions_run")) next
  message(sprintf("  %-30s %8.1f s", nm, timing_log[[nm]]))
}

if (!is.null(sample_frac)) {
  scale_factor <- n_total_original / total_births
  est_total <- timing_log$total_sec * scale_factor
  message("\nEstimado cohorte completa (extrapolación lineal en N): ",
          round(est_total / 60, 1), " min (", round(est_total / 3600, 2), " h)")
  message("  (factor de escala N: ", round(scale_factor, 1), "x)")
}

message("\nOutputs guardados en: ", data_out_g)
message("Log de tiempos: ", timing_path)

beepr::beep(8)
