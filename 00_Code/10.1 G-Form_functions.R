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
  boot_iter = 500L,
  boot_seed = 2026L,
  figure4_pct = 0.20,
  output_cumulative_risk_suffix = "curvas_riesgo_acumulado_global",
  output_singleweek_heatmap_suffix = "mapa_calor_rd_semana_intervencion"
)

## Paralelismo (servidor Linux / local) ----

gform_detect_ram_gb <- function() {
  parse_kb <- function(line) {
    if (!length(line) || !nzchar(line)) return(NA_real_)
    m <- regexpr("[0-9]+", line, perl = TRUE)
    if (m[1L] < 0) return(NA_real_)
    kb <- suppressWarnings(as.numeric(substring(line, m[1L], m[1L] + attr(m, "match.length") - 1L)))
    if (!is.finite(kb)) return(NA_real_)
    kb / 1024 / 1024
  }

  if (.Platform$OS.type == "unix" && file.exists("/proc/meminfo")) {
    lines <- readLines("/proc/meminfo", warn = FALSE)
    memtotal <- grep("^MemTotal:", lines, value = TRUE)
    gb <- parse_kb(memtotal[1L])
    if (is.finite(gb)) return(gb)
  }

  gb <- tryCatch({
    out <- system2("getconf", c("_PHYS_PAGES"), stdout = TRUE, stderr = FALSE)
    psize <- system2("getconf", c("PAGE_SIZE"), stdout = TRUE, stderr = FALSE)
    if (length(out) && length(psize)) {
      as.numeric(out[1]) * as.numeric(psize[1]) / 1024^3
    } else {
      NA_real_
    }
  }, error = function(e) NA_real_)
  if (is.finite(gb)) return(gb)

  gb <- tryCatch({
    out <- system2("free", c("-g", "--si"), stdout = TRUE, stderr = FALSE)
    mem_line <- grep("^Mem:", out, value = TRUE)
    if (!length(mem_line)) return(NA_real_)
    parts <- strsplit(trimws(mem_line[1L]), "\\s+")[[1L]]
    if (length(parts) >= 2L) as.numeric(parts[2L]) else NA_real_
  }, error = function(e) NA_real_)
  if (is.finite(gb)) return(gb)

  NA_real_
}

gform_env_num <- function(name, default = NA_real_) {
  val <- Sys.getenv(name, unset = "")
  if (!nzchar(val)) return(default)
  out <- suppressWarnings(as.numeric(val))
  if (is.finite(out)) out else default
}

gform_env_bool <- function(name, default = FALSE) {
  val <- tolower(trimws(Sys.getenv(name, unset = "")))
  if (!nzchar(val)) return(isTRUE(default))
  val %in% c("1", "true", "yes", "on")
}

gform_env_int_vec <- function(name, default = NULL) {
  val <- trimws(Sys.getenv(name, unset = ""))
  if (!nzchar(val)) return(default)
  parts <- strsplit(val, "[,;\\s]+")[[1L]]
  parts <- parts[nzchar(parts)]
  out <- suppressWarnings(as.integer(parts))
  out <- out[!is.na(out)]
  if (!length(out)) default else out
}

gform_run_fingerprint <- function(
    output_stub,
    boot_iter,
    boot_seed,
    total_births,
    sample_frac = NULL) {
  list(
    output_stub = output_stub,
    boot_iter = as.integer(boot_iter),
    boot_seed = as.integer(boot_seed),
    n_births = as.integer(total_births),
    sample_frac = sample_frac
  )
}

gform_fingerprint_matches <- function(saved, current) {
  if (is.null(saved) || is.null(current)) return(FALSE)
  keys <- c("output_stub", "boot_iter", "boot_seed", "n_births", "sample_frac")
  for (k in keys) {
    a <- saved[[k]]
    b <- current[[k]]
    if (identical(a, NULL) && identical(b, NULL)) next
    if (is.numeric(a) && is.numeric(b)) {
      if (!isTRUE(all.equal(a, b))) return(FALSE)
    } else if (!identical(a, b)) {
      return(FALSE)
    }
  }
  TRUE
}

gform_point_checkpoint_path <- function(output_stub, dir_bootstrap) {
  file.path(dir_bootstrap, output_stub, "point_checkpoint.rds")
}

gform_save_point_checkpoint <- function(
    path,
    fingerprint,
    weekly_effects,
    population_effects,
    cumulative_risk_curves,
    nat_mean) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(c(
    fingerprint,
    list(
      weekly_effects = weekly_effects,
      population_effects = population_effects,
      cumulative_risk_curves = cumulative_risk_curves,
      nat_mean = nat_mean,
      saved_at = Sys.time()
    )
  ), path)
  invisible(path)
}

gform_read_point_checkpoint <- function(path, fingerprint) {
  if (!file.exists(path)) return(NULL)
  pt <- tryCatch(readRDS(path), error = function(e) NULL)
  if (is.null(pt)) return(NULL)
  fp <- pt[c("output_stub", "boot_iter", "boot_seed", "n_births", "sample_frac")]
  if (!gform_fingerprint_matches(fp, fingerprint)) {
    message("Checkpoint de punto estimado incompatible; se recalculará.")
    return(NULL)
  }
  pt
}

gform_bootstrap_ck_matches <- function(ck, fingerprint) {
  if (is.null(ck)) return(FALSE)
  if (!"output_stub" %in% names(ck)) {
    return(
      identical(as.integer(ck$boot_iter), as.integer(fingerprint$boot_iter)) &&
        identical(as.integer(ck$boot_seed), as.integer(fingerprint$boot_seed))
    )
  }
  fp <- ck[c("output_stub", "boot_iter", "boot_seed", "n_births", "sample_frac")]
  gform_fingerprint_matches(fp, fingerprint)
}

gform_save_bootstrap_checkpoint <- function(path, fingerprint, last_completed) {
  saveRDS(c(
    fingerprint,
    list(
      last_completed = as.integer(last_completed),
      updated_at = Sys.time()
    )
  ), path)
  invisible(path)
}

gform_clear_intervention_checkpoints <- function(output_stub, dir_bootstrap) {
  paths <- gform_bootstrap_paths(output_stub, dir_bootstrap)
  for (p in c(
    paths$checkpoint,
    paths$weekly,
    paths$population,
    gform_point_checkpoint_path(output_stub, dir_bootstrap)
  )) {
    if (file.exists(p)) file.remove(p)
  }
  invisible(NULL)
}

gform_subsample_data_long <- function(data_long, sample_frac, sample_seed = 2026L) {
  if (is.null(sample_frac) || !is.finite(sample_frac) ||
      sample_frac <= 0 || sample_frac >= 1) {
    return(data_long)
  }
  ids <- unique(data_long$id)
  n_keep <- max(1L, floor(length(ids) * sample_frac))
  set.seed(sample_seed)
  ids_keep <- sample(ids, n_keep)
  out <- data_long[data_long$id %in% ids_keep, , drop = FALSE]
  message(
    "Submuestra ", round(100 * sample_frac, 2), "%: ",
    n_keep, " / ", length(ids), " nacimientos"
  )
  out
}

gform_parallel_config <- function(
    n_cores = parallel::detectCores(logical = TRUE),
    ram_gb = gform_detect_ram_gb(),
    reserve_cores = 4L,
    reserve_ram_gb = NULL) {

  if (!is.finite(n_cores) || n_cores < 1L) n_cores <- 4L
  if (!is.finite(ram_gb)) {
    ram_gb <- gform_env_num("GFORM_RAM_GB", 16L)
    if (!is.finite(ram_gb)) ram_gb <- 16L
    message("RAM no detectada; usando fallback ", ram_gb, " GiB (export GFORM_RAM_GB para fijar).")
  } else {
    ram_override <- gform_env_num("GFORM_RAM_GB", NA_real_)
    if (is.finite(ram_override)) ram_gb <- ram_override
  }
  if (is.null(reserve_ram_gb)) {
    reserve_ram_gb <- as.integer(gform_env_num("GFORM_RESERVE_RAM_GB", 24L))
  }

  n_usable <- max(1L, n_cores - reserve_cores)
  ram_usable <- max(8L, ram_gb - reserve_ram_gb)

  n_cox <- min(9L, n_usable)
  n_build <- min(11L, n_usable)

  # Heatmap (fork): cada worker corre predict_weekly_hazards sobre cohorte completa (~12+ GiB).
  heatmap_ram_per_worker <- as.integer(
    gform_env_num("GFORM_HEATMAP_RAM_PER_WORKER_GB", 14L)
  )
  heatmap_max_workers <- as.integer(
    gform_env_num("GFORM_HEATMAP_MAX_WORKERS", 4L)
  )
  heatmap_workers_override <- gform_env_num("GFORM_HEATMAP_WORKERS", NA_real_)
  heatmap_parent_reserve_gb <- max(
    as.integer(reserve_ram_gb),
    40L,
    floor(ram_gb * 0.40)
  )
  n_heatmap_by_ram <- max(
    1L,
    floor((ram_gb - heatmap_parent_reserve_gb) / heatmap_ram_per_worker)
  )
  n_heatmap <- min(
    heatmap_max_workers,
    n_heatmap_by_ram,
    max(1L, floor(n_usable * 0.25))
  )
  if (is.finite(heatmap_workers_override)) {
    n_heatmap <- max(1L, min(as.integer(heatmap_workers_override), n_usable))
  }
  heatmap_batch_size <- as.integer(
    gform_env_num("GFORM_HEATMAP_BATCH_SIZE", n_heatmap)
  )
  heatmap_batch_size <- max(1L, min(heatmap_batch_size, n_heatmap))

  # Bootstrap (fork): reservar RAM para el padre (model_store + frames; COW bajo carga).
  bootstrap_ram_per_worker <- as.integer(
    gform_env_num("GFORM_BOOTSTRAP_RAM_PER_WORKER_GB", 8L)
  )
  bootstrap_max_workers <- as.integer(
    gform_env_num("GFORM_BOOTSTRAP_MAX_WORKERS", 8L)
  )
  bootstrap_workers_override <- gform_env_num("GFORM_BOOTSTRAP_WORKERS", NA_real_)
  parent_reserve_gb <- max(
    as.integer(reserve_ram_gb),
    36L,
    floor(ram_gb * 0.35)
  )
  n_bootstrap_by_ram <- max(
    1L,
    floor((ram_gb - parent_reserve_gb) / bootstrap_ram_per_worker)
  )
  n_bootstrap <- min(
    bootstrap_max_workers,
    n_bootstrap_by_ram,
    max(2L, floor(n_usable * 0.75))
  )
  if (is.finite(bootstrap_workers_override)) {
    n_bootstrap <- max(1L, min(as.integer(bootstrap_workers_override), n_usable))
  }

  globals_env <- gform_env_num("GFORM_GLOBALS_MAX_GB", NA_real_)
  globals_max_gb <- if (is.finite(globals_env)) {
    as.integer(globals_env)
  } else if (gform_is_linux_server()) {
    min(128L, max(64L, ceiling(ram_gb * 0.85)))
  } else {
    min(80L, max(48L, floor(ram_usable * 0.65)))
  }
  globals_max_gb_bootstrap <- if (is.finite(globals_env)) {
    as.integer(globals_env)
  } else if (gform_is_linux_server()) {
    # Con multicore/fork el chequeo de tamaño es conservador; objetos ~40 GiB son normales.
    min(128L, max(64L, ceiling(ram_gb * 0.90)))
  } else {
    globals_max_gb
  }

  list(
    n_cores = n_cores,
    ram_gb = ram_gb,
    n_workers_cox = n_cox,
    n_workers_bootstrap = n_bootstrap,
    n_workers_heatmap = n_heatmap,
    n_workers_build = n_build,
    n_workers_default = min(n_bootstrap, n_usable),
    globals_max_gb = globals_max_gb,
    globals_max_gb_bootstrap = globals_max_gb_bootstrap,
    bootstrap_batch_size = n_bootstrap,
    bootstrap_parent_reserve_gb = parent_reserve_gb,
    bootstrap_ram_per_worker_gb = bootstrap_ram_per_worker,
    heatmap_batch_size = heatmap_batch_size,
    heatmap_parent_reserve_gb = heatmap_parent_reserve_gb,
    heatmap_ram_per_worker_gb = heatmap_ram_per_worker,
    use_fork = gform_is_linux_server(),
    dt_threads = max(1L, n_cores - 2L)
  )
}

gform_is_linux_server <- function() {
  .Platform$OS.type == "unix" &&
    identical(Sys.info()[["sysname"]], "Linux") &&
    parallel::detectCores(logical = TRUE) >= 8L
}

gform_setup_parallel <- function(
    n_workers = NULL,
    globals_max_gb = NULL,
    task = c("default", "cox", "bootstrap", "heatmap", "build"),
    config = getOption("gform.parallel", NULL)) {

  task <- match.arg(task)
  if (is.null(config)) config <- gform_parallel_config()
  options(gform.parallel = config)

  if (is.null(n_workers)) {
    n_workers <- switch(
      task,
      cox = config$n_workers_cox,
      bootstrap = config$n_workers_bootstrap,
      heatmap = config$n_workers_heatmap,
      build = config$n_workers_build,
      default = config$n_workers_default
    )
  }
  if (is.null(globals_max_gb)) {
    globals_max_gb <- if (task == "bootstrap" && !is.null(config$globals_max_gb_bootstrap)) {
      config$globals_max_gb_bootstrap
    } else {
      config$globals_max_gb
    }
  }

  # Linux fork: el límite de 16 GiB de future bloquea model_store ~40 GiB sin copiarlos.
  if (task == "bootstrap" && isTRUE(config$use_fork)) {
    globals_max_gb <- max(globals_max_gb, 64L)
    options(future.globals.maxSize = globals_max_gb * 1024^3)
  } else {
    options(future.globals.maxSize = globals_max_gb * 1024^3)
  }

  if (isTRUE(config$use_fork) && requireNamespace("future", quietly = TRUE)) {
    future::plan(future::multicore, workers = n_workers)
    plan_label <- paste0("multicore/fork (", n_workers, " workers)")
  } else if (requireNamespace("future", quietly = TRUE)) {
    future::plan(future::multisession, workers = n_workers)
    plan_label <- paste0("multisession (", n_workers, " workers)")
  }

  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::setDTthreads(config$dt_threads)
  }

  message(
    "Paralelo [", task, "]: ", plan_label,
    " | globals max: ", globals_max_gb, " GiB",
    " | data.table threads: ", config$dt_threads
  )
  invisible(list(config = config, n_workers = n_workers, task = task))
}

## Registro de tiempos por fase ----

gform_format_timestamp <- function(x = Sys.time()) {
  format(x, "%Y-%m-%d %H:%M:%S")
}

gform_format_duration <- function(sec) {
  if (!is.finite(sec) || sec < 0) return("NA")
  if (sec < 60) return(sprintf("%.1f s", sec))
  if (sec < 3600) return(sprintf("%.1f min (%.0f s)", sec / 60, sec))
  sprintf("%.2f h (%.1f min)", sec / 3600, sec / 60)
}

gform_timing_log_init <- function() {
  list(phases = list(), started_at = Sys.time())
}

gform_phase_start <- function(label) {
  t0 <- Sys.time()
  message("[", gform_format_timestamp(t0), "] INICIO — ", label)
  list(label = label, start = t0)
}

gform_phase_end <- function(phase) {
  t1 <- Sys.time()
  sec <- as.numeric(difftime(t1, phase$start, units = "secs"))
  message(
    "[", gform_format_timestamp(t1), "] FIN — ", phase$label,
    " | duración: ", gform_format_duration(sec)
  )
  list(
    label = phase$label,
    start = phase$start,
    end = t1,
    start_str = gform_format_timestamp(phase$start),
    end_str = gform_format_timestamp(t1),
    sec = sec
  )
}

gform_timing_log_add <- function(log, record) {
  log$phases <- c(log$phases, list(record))
  log
}

gform_time_block <- function(label, expr) {
  phase <- gform_phase_start(label)
  result <- force(expr)
  record <- gform_phase_end(phase)
  list(result = result, timing = record)
}

gform_timing_total_sec <- function(log) {
  if (!length(log$phases)) return(0)
  sum(vapply(log$phases, function(p) as.numeric(p$sec), numeric(1)))
}

gform_timing_log_merge <- function(parent, child) {
  if (is.null(child) || !length(child$phases)) return(parent)
  for (p in child$phases) {
    parent <- gform_timing_log_add(parent, p)
  }
  parent
}

gform_print_timing_summary <- function(log, title = "Resumen de tiempos") {
  if (!length(log$phases)) {
    message(title, ": (sin fases registradas)")
    return(invisible(0))
  }
  total_sec <- gform_timing_total_sec(log)
  message("\n", strrep("-", 72))
  message(title)
  message(strrep("-", 72))
  for (p in log$phases) {
    message(sprintf(
      "  %-42s %10s  (%s → %s)",
      p$label,
      gform_format_duration(p$sec),
      p$start_str,
      p$end_str
    ))
  }
  message(strrep("-", 72))
  wall_sec <- log$wall_sec
  if (!is.null(wall_sec) && is.finite(wall_sec)) {
    message(sprintf("  %-42s %10s", "TOTAL (suma fases)", gform_format_duration(total_sec)))
    message(sprintf(
      "  %-42s %10s  (%s → %s)",
      "TOTAL (reloj)",
      gform_format_duration(wall_sec),
      gform_format_timestamp(log$started_at),
      gform_format_timestamp(if (!is.null(log$finished_at)) log$finished_at else Sys.time())
    ))
    if (total_sec > wall_sec * 1.05) {
      message("  Nota: en ejecución paralela, la suma de fases puede superar el reloj.")
    }
  } else {
    message(sprintf("  %-42s %10s", "TOTAL", gform_format_duration(total_sec)))
  }
  message(strrep("-", 72))
  invisible(if (!is.null(wall_sec) && is.finite(wall_sec)) wall_sec else total_sec)
}

## Orden de ejecución (una intervención por corrida; ver 10.2) ----

GFORM_INTERVENTION_ORDER <- c(
  "pm25_krg_pct20",
  "no2_krg_pct20",
  "o3_krg_pct20",
  "pm25_krg_lt20",
  "pm25_krg_lt5",
  "no2_krg_lt20",
  "no2_krg_lt5",
  "pm25_krg_lt15",
  "pm25_krg_lt10",
  "no2_krg_lt15",
  "no2_krg_lt10"
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

gform_intervention_keys <- function() {
  missing <- setdiff(GFORM_INTERVENTION_ORDER, names(GFORM_INTERVENTION_REGISTRY))
  if (length(missing)) {
    stop("GFORM_INTERVENTION_ORDER contiene IDs ausentes en el registro: ",
         paste(missing, collapse = ", "))
  }
  GFORM_INTERVENTION_ORDER
}

print_gform_intervention_menu <- function(selected = NULL) {
  keys <- gform_intervention_keys()
  message("Intervenciones disponibles (cambiar intervention_number en 10.2):")
  for (i in seq_along(keys)) {
    spec <- GFORM_INTERVENTION_REGISTRY[[keys[i]]]
    marker <- if (!is.null(selected) && i == selected) "  <-- seleccionada" else ""
    message(sprintf(
      "  %2d. %s — %s%s",
      i, spec$intervention_id, spec$description, marker
    ))
  }
  invisible(keys)
}

resolve_gform_intervention <- function(n) {
  keys <- gform_intervention_keys()
  n <- as.integer(n[[1L]])
  if (length(n) != 1L || is.na(n) || n < 1L || n > length(keys)) {
    stop(
      "intervention_number debe ser un entero entre 1 y ", length(keys),
      ". Ejecute print_gform_intervention_menu() para ver la lista."
    )
  }
  keys[[n]]
}

gform_output_paths <- function(output_stub, dir_other) {
  list(
    cumulative_risk_curves = file.path(
      dir_other,
      paste0(output_stub, "_", GFORM_DEFAULTS$output_cumulative_risk_suffix, ".rds")
    ),
    singleweek_heatmap = file.path(
      dir_other,
      paste0(output_stub, "_", GFORM_DEFAULTS$output_singleweek_heatmap_suffix, ".rds")
    )
  )
}

gform_model_cache_path <- function(pollutant, dir_models) {
  file.path(dir_models, paste0("natural_course_", pollutant, ".rds"))
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

  list(
    n_births = nrow(data_base),
    pollutant = pollutant,
    risk_weeks = risk_weeks_vec,
    control_vars = control_vars,
    lag_weeks = lag_weeks,
    dependent_var = dependent_var,
    risk_entry_week = risk_entry_week,
    sample_frac = sample_frac
  )
}

load_or_fit_natural_course_models <- function(
    pollutant,
    dir_models,
    data_base,
    wide_exposicion_natural,
    wide_tad_obs,
    risk_weeks_vec,
    control_vars = GFORM_DEFAULTS$control_vars,
    lag_weeks = GFORM_DEFAULTS$lag_weeks,
    dependent_var = GFORM_DEFAULTS$dependent_var,
    risk_entry_week = GFORM_DEFAULTS$risk_entry_week,
    sample_frac = NULL,
    parallel = FALSE,
    force_refit = FALSE) {

  dir.create(dir_models, recursive = TRUE, showWarnings = FALSE)
  cache_path <- gform_model_cache_path(pollutant, dir_models)
  cache_key <- gform_natural_course_cache_key(
    data_base = data_base,
    pollutant = pollutant,
    risk_weeks_vec = risk_weeks_vec,
    control_vars = control_vars,
    lag_weeks = lag_weeks,
    dependent_var = dependent_var,
    risk_entry_week = risk_entry_week,
    sample_frac = sample_frac
  )

  if (!force_refit && file.exists(cache_path)) {
    cached <- tryCatch(
      readRDS(cache_path),
      error = function(e) {
        message("Cache Cox ilegible (", conditionMessage(e), "); re-ajustando modelos...")
        NULL
      }
    )
    if (!is.null(cached) && identical(cached$cache_key, cache_key)) {
      message("Modelos Cox cacheados cargados: ", cache_path)
      cached$model_store <- slim_model_store(cached$model_store)
      if (!isTRUE(cached$slimmed)) {
        message("Cache Cox marcado slim (re-guardado omitido para ahorrar RAM/tiempo).")
      }
      return(list(
        model_store = cached$model_store,
        cox_frame_natural = cached$cox_frame_natural,
        from_cache = TRUE
      ))
    }
    if (!is.null(cached)) {
      message("Cache Cox desactualizado; re-ajustando modelos...")
    }
  }

  fit <- fit_natural_course_models(
    data_base = data_base,
    wide_exposicion_natural = wide_exposicion_natural,
    wide_tad_obs = wide_tad_obs,
    risk_weeks_vec = risk_weeks_vec,
    control_vars = control_vars,
    lag_weeks = lag_weeks,
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
  message("Modelos Cox guardados en cache: ", cache_path)
  c(fit, list(from_cache = FALSE))
}

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
  prep_block <- gform_time_block(
    "Preparar matrices raw por contaminante",
    {
      raw_by_pollutant <- lapply(
        pollutants_needed,
        function(p) build_wide_raw_exposure(data_long, p, weeks_keep = weeks_keep)
      )
      names(raw_by_pollutant) <- pollutants_needed
      raw_by_pollutant
    }
  )
  raw_by_pollutant <- prep_block$result

  build_one <- function(key) {
    spec <- registry[[key]]
    label <- paste0("RDS — ", spec$intervention_id)
    block <- gform_time_block(label, {
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
  timing <- gform_timing_log_add(timing, prep_block$timing)
  for (item in built) {
    timing <- gform_timing_log_add(timing, item$timing)
  }

  list(
    paths = vapply(built, function(x) x$path, character(1)),
    timing = timing
  )
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

align_newdata_to_cox_fit <- function(cox_fit, newdata_df) {
  out <- newdata_df
  xlevels <- cox_fit$xlevels
  if (length(xlevels)) {
    for (nm in names(xlevels)) {
      if (!nm %in% names(out)) next
      out[[nm]] <- factor(out[[nm]], levels = xlevels[[nm]])
    }
  }
  out
}

predict_cox_lp <- function(cox_fit, newdata, coef_override = NULL) {
  newdata_df <- align_newdata_to_cox_fit(cox_fit, as.data.frame(newdata))
  b <- stats::coef(cox_fit)
  if (!is.null(coef_override)) {
    ii <- intersect(names(b), names(coef_override))
    b[ii] <- coef_override[ii]
  }
  b[is.na(b)] <- 0
  mm <- stats::model.matrix(
    stats::delete.response(stats::terms(cox_fit)),
    newdata_df,
    xlev = cox_fit$xlevels
  )
  if ("(Intercept)" %in% colnames(mm)) {
    mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
  }
  cols <- intersect(colnames(mm), names(b))
  as.numeric(drop(mm[, cols, drop = FALSE] %*% b[cols]))
}

## Etapa 2: ajuste Cox (solo curso natural; referencia 9.0) ----

fit_cox_one_week <- function(
    rw,
    data_model_df,
    control_vars,
    lag_weeks,
    dependent_var,
    risk_entry_week) {

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

slim_coxph_model_store_entry <- function(entry) {
  if (is.null(entry) || is.null(entry$model)) return(entry)
  m <- entry$model
  m$x <- NULL
  m$y <- NULL
  m$residuals <- NULL
  m$linear.predictors <- NULL
  m$means <- NULL
  m$weights <- NULL
  m$offset <- NULL
  m$model <- NULL
  m$na.action <- NULL
  m$call <- NULL
  entry$model <- m
  entry
}

slim_model_store <- function(model_store) {
  if (is.null(model_store)) return(model_store)
  out <- lapply(model_store, slim_coxph_model_store_entry)
  names(out) <- names(model_store)
  out
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
  data_model_df <- as.data.frame(cox_frame_natural)

  weeks_to_fit <- risk_weeks_vec[risk_weeks_vec >= risk_entry_week]
  if (parallel && length(weeks_to_fit) > 1L && requireNamespace("furrr", quietly = TRUE)) {
    gform_setup_parallel(task = "cox")
    fitted <- furrr::future_map(
      weeks_to_fit,
      fit_cox_one_week,
      data_model_df = data_model_df,
      control_vars = control_vars,
      lag_weeks = lag_weeks,
      dependent_var = dependent_var,
      risk_entry_week = risk_entry_week,
      .options = gform_furrr_options()
    )
  } else {
    fitted <- lapply(weeks_to_fit, function(rw) {
      fit_cox_one_week(
        rw, data_model_df, control_vars, lag_weeks, dependent_var, risk_entry_week
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
    cox_frame_natural = cox_frame_natural
  )
}

## Predicción de hazards semanales (Cox → probabilidad discreta) ----

predict_weekly_hazards <- function(
    model_store,
    person_weeks,
    risk_weeks_vec,
    risk_entry_week = GFORM_DEFAULTS$risk_entry_week,
    coef_override = NULL,
    cox_frame = NULL,
    data_base = NULL,
    wide_exposicion = NULL,
    wide_tad_obs = NULL,
    control_vars = GFORM_DEFAULTS$control_vars,
    dependent_var = GFORM_DEFAULTS$dependent_var) {

  if (is.null(cox_frame)) {
    if (is.null(data_base) || is.null(wide_exposicion) || is.null(wide_tad_obs)) {
      stop("predict_weekly_hazards: se requiere cox_frame o data_base + wide_exposicion + wide_tad_obs.")
    }
    cox_frame <- build_cox_model_frame(
      data_base = data_base,
      wide_exposicion = wide_exposicion,
      wide_tad_obs = wide_tad_obs,
      control_vars = control_vars,
      dependent_var = dependent_var,
      risk_entry_week = risk_entry_week
    )
  } else if (!data.table::is.data.table(cox_frame)) {
    cox_frame <- data.table::as.data.table(cox_frame)
    data.table::setkey(cox_frame, id)
  }

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

compute_cumulative_risk_curves_global <- function(
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
  out <- tibble::as_tibble(dt[, .(
    follow_up_week,
    cumulative_risk_observed,
    cumulative_risk_intervention,
    risk_ratio,
    risk_difference
  )])
  attr(out, "description") <- paste(
    "Riesgo acumulado medio de PTB por semana de seguimiento:",
    "curso natural (exposicion observada) vs intervencion global en todas las semanas gestacionales."
  )
  out
}

compute_figure3_curves <- compute_cumulative_risk_curves_global

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

  baseline_risk <- as.numeric(risk_df[scenario == baseline_scenario, prevalence][1L])
  paf <- if (baseline_risk > 0) {
    (baseline_risk - risk_df$prevalence) / baseline_risk
  } else {
    rep(NA_real_, nrow(risk_df))
  }

  risk_df[, `:=`(
    risk_ratio = prevalence / baseline_risk,
    risk_difference = prevalence - baseline_risk,
    attributable_risk = baseline_risk - prevalence,
    attributable_fraction = data.table::fifelse(
      scenario == baseline_scenario,
      0,
      paf
    )
  )]

  tibble::as_tibble(risk_df[, .(
    scenario,
    prevalence,
    cases,
    risk_ratio,
    risk_difference,
    attributable_risk,
    attributable_fraction
  )])
}

compute_singleweek_intervention_heatmap <- function(
    model_store,
    person_weeks,
    data_base,
    raw_wide_pollutant,
    pollutant,
    wide_tad_obs,
    risk_weeks_vec,
    control_vars,
    nat_mean = NULL,
    intervention_weeks = GFORM_DEFAULTS$weeks_exposure,
    follow_up_weeks = GFORM_DEFAULTS$follow_up_weeks,
    pct = GFORM_DEFAULTS$figure4_pct,
    dependent_var = GFORM_DEFAULTS$dependent_var,
    risk_entry_week = GFORM_DEFAULTS$risk_entry_week,
    output_stub,
    dir_heatmap,
    total_births,
    sample_frac = NULL,
    resume = TRUE,
    parallel = FALSE,
    heatmap_batch_size = NULL) {

  single_intervention <- list(type = "pct_reduce", pct = pct)
  n_cols <- length(intervention_weeks)

  if (is.null(nat_mean)) {
    stop("compute_singleweek_intervention_heatmap: se requiere nat_mean (riesgo acumulado medio bajo curso natural).")
  }
  if (is.null(output_stub) || is.null(dir_heatmap)) {
    stop("compute_singleweek_intervention_heatmap: se requiere output_stub y dir_heatmap.")
  }

  paths <- gform_heatmap_paths(output_stub, dir_heatmap)
  dir.create(paths$dir, recursive = TRUE, showWarnings = FALSE)
  fingerprint <- gform_heatmap_fingerprint(
    output_stub = output_stub,
    pct = pct,
    total_births = total_births,
    sample_frac = sample_frac,
    follow_up_weeks = follow_up_weeks
  )

  compute_one_col <- function(iw) {
    wide_cf <- build_exposicion_wide_from_raw(
      raw_wide = raw_wide_pollutant,
      pollutant = pollutant,
      intervention = single_intervention,
      single_week = iw
    )
    cox_frame_iw <- build_cox_model_frame(
      data_base = data_base,
      wide_exposicion = wide_cf,
      wide_tad_obs = wide_tad_obs,
      control_vars = control_vars,
      dependent_var = dependent_var,
      risk_entry_week = risk_entry_week
    )
    prob_int <- predict_weekly_hazards(
      model_store = model_store,
      person_weeks = person_weeks,
      risk_weeks_vec = risk_weeks_vec,
      cox_frame = cox_frame_iw,
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

  finalize_heatmap <- function(long_df) {
    wide_mat <- data.table::dcast(
      long_df,
      follow_up_week ~ intervention_week,
      value.var = "risk_difference"
    )
    out <- list(
      long = tibble::as_tibble(long_df),
      wide = tibble::as_tibble(wide_mat)
    )
    attr(out, "description") <- paste0(
      "Mapa de calor RD acumulado de PTB: reduccion del ",
      round(100 * pct), "% en una sola semana gestacional (eje X) ",
      "evaluado en semanas de seguimiento ", min(follow_up_weeks), "-", max(follow_up_weeks),
      " (eje Y), referencia curso natural comun."
    )
    out
  }

  append_one_col <- function(piece) {
    append_mode <- file.exists(paths$long)
    data.table::fwrite(
      piece, paths$long,
      append = append_mode,
      col.names = !append_mode
    )
  }

  completed_weeks <- integer(0)
  if (isTRUE(resume) && file.exists(paths$checkpoint)) {
    ck <- tryCatch(readRDS(paths$checkpoint), error = function(e) NULL)
    if (!is.null(ck) && gform_heatmap_ck_matches(ck, fingerprint)) {
      completed_weeks <- as.integer(ck$completed_weeks)
      if (length(completed_weeks) >= n_cols &&
          all(intervention_weeks %in% completed_weeks) &&
          file.exists(paths$long)) {
        message("Heatmap ya completo (", n_cols, " columnas); leyendo desde disco.")
        long_df <- data.table::fread(paths$long)
        return(finalize_heatmap(long_df))
      }
      if (length(completed_weeks)) {
        next_col <- length(completed_weeks) + 1L
        message(
          "Heatmap: reanudando desde columna ", next_col, " / ", n_cols,
          " (", length(completed_weeks), " columnas en checkpoint)"
        )
      }
    } else {
      message("Checkpoint heatmap incompatible; reiniciando mapa de calor.")
      completed_weeks <- integer(0)
    }
  }

  weeks_todo <- setdiff(intervention_weeks, completed_weeks)
  if (!length(completed_weeks)) {
    if (file.exists(paths$long)) file.remove(paths$long)
    if (file.exists(paths$checkpoint)) file.remove(paths$checkpoint)
  }

  if (!length(weeks_todo)) {
    long_df <- data.table::fread(paths$long)
    return(finalize_heatmap(long_df))
  }

  cfg <- getOption("gform.parallel", gform_parallel_config())
  if (is.null(heatmap_batch_size)) {
    heatmap_batch_size <- if (parallel) cfg$heatmap_batch_size else 1L
  }
  heatmap_batch_size <- max(1L, as.integer(heatmap_batch_size))

  col_index <- function(iw) {
    match(iw, intervention_weeks)
  }

  if (parallel && length(weeks_todo) > 0L && requireNamespace("furrr", quietly = TRUE)) {
    gform_setup_parallel(task = "heatmap", config = cfg)
    message(
      "Heatmap paralelo: lotes de ", heatmap_batch_size,
      " columnas (", cfg$n_workers_heatmap, " workers)"
    )
    batch_starts <- seq(1L, length(weeks_todo), by = heatmap_batch_size)
    for (batch_start in batch_starts) {
      batch_end <- min(batch_start + heatmap_batch_size - 1L, length(weeks_todo))
      batch_weeks <- weeks_todo[batch_start:batch_end]
      batch_results <- furrr::future_map(
        batch_weeks,
        compute_one_col,
        .options = gform_furrr_options()
      )
      for (k in seq_along(batch_weeks)) {
        append_one_col(batch_results[[k]])
        message("Heatmap: columna ", col_index(batch_weeks[k]), " / ", n_cols)
      }
      rm(batch_results)
      gc(verbose = FALSE)
      completed_weeks <- sort(unique(c(completed_weeks, batch_weeks)))
      gform_save_heatmap_checkpoint(paths$checkpoint, fingerprint, completed_weeks)
    }
    future::plan(future::sequential)
  } else {
    if (isTRUE(parallel)) {
      message("Heatmap paralelo no disponible (furrr); ejecutando secuencial.")
    }
    for (iw in weeks_todo) {
      message("Heatmap: columna ", col_index(iw), " / ", n_cols)
      piece <- compute_one_col(iw)
      append_one_col(piece)
      rm(piece)
      gc(verbose = FALSE)
      completed_weeks <- sort(unique(c(completed_weeks, iw)))
      gform_save_heatmap_checkpoint(paths$checkpoint, fingerprint, completed_weeks)
    }
  }

  long_df <- data.table::fread(paths$long)
  finalize_heatmap(long_df)
}

compute_figure4_heatmap <- compute_singleweek_intervention_heatmap

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

bootstrap_one_iter <- function(
    iter,
    model_store,
    person_weeks,
    cox_frame_natural,
    cox_frame_intervention,
    risk_weeks_vec,
    total_births,
    dependent_var,
    risk_entry_week,
    target_week,
    baseline_scenario,
    intervention_scenario,
    boot_seed = GFORM_DEFAULTS$boot_seed) {

  set.seed(boot_seed + iter)
  coefs_list <- simulate_coefs_list(model_store, risk_weeks_vec)
  prob_nat <- predict_weekly_hazards(
    model_store = model_store,
    person_weeks = person_weeks,
    risk_weeks_vec = risk_weeks_vec,
    cox_frame = cox_frame_natural,
    risk_entry_week = risk_entry_week,
    coef_override = coefs_list
  )
  prob_int <- predict_weekly_hazards(
    model_store = model_store,
    person_weeks = person_weeks,
    risk_weeks_vec = risk_weeks_vec,
    cox_frame = cox_frame_intervention,
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

gform_bootstrap_paths <- function(output_stub, dir_bootstrap) {
  base <- file.path(dir_bootstrap, output_stub)
  list(
    dir = base,
    weekly = file.path(base, "weekly_boot.csv"),
    population = file.path(base, "population_boot.csv"),
    checkpoint = file.path(base, "boot_checkpoint.rds")
  )
}

gform_heatmap_paths <- function(output_stub, dir_heatmap) {
  base <- file.path(dir_heatmap, output_stub)
  list(
    dir = base,
    long = file.path(base, "heatmap_long.csv"),
    checkpoint = file.path(base, "heatmap_checkpoint.rds")
  )
}

gform_heatmap_fingerprint <- function(
    output_stub,
    pct,
    total_births,
    sample_frac,
    follow_up_weeks) {
  list(
    output_stub = output_stub,
    pct = pct,
    n_births = as.integer(total_births),
    sample_frac = sample_frac,
    follow_up_min = min(follow_up_weeks),
    follow_up_max = max(follow_up_weeks)
  )
}

gform_heatmap_ck_matches <- function(ck, fingerprint) {
  if (is.null(ck) || is.null(fingerprint)) return(FALSE)
  keys <- c("output_stub", "pct", "n_births", "sample_frac", "follow_up_min", "follow_up_max")
  for (k in keys) {
    a <- ck[[k]]
    b <- fingerprint[[k]]
    if (identical(a, NULL) && identical(b, NULL)) next
    if (is.numeric(a) && is.numeric(b)) {
      if (!isTRUE(all.equal(a, b))) return(FALSE)
    } else if (!identical(a, b)) {
      return(FALSE)
    }
  }
  TRUE
}

gform_save_heatmap_checkpoint <- function(path, fingerprint, completed_weeks) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(c(
    fingerprint,
    list(
      completed_weeks = sort(unique(as.integer(completed_weeks))),
      updated_at = Sys.time()
    )
  ), path)
  invisible(path)
}

gform_bootstrap_is_complete <- function(
    output_stub,
    dir_bootstrap,
    boot_iter,
    boot_seed,
    total_births,
    sample_frac) {
  if (boot_iter <= 0L) return(TRUE)
  paths <- gform_bootstrap_paths(output_stub, dir_bootstrap)
  if (!file.exists(paths$checkpoint)) return(FALSE)
  ck <- tryCatch(readRDS(paths$checkpoint), error = function(e) NULL)
  if (is.null(ck)) return(FALSE)
  fp <- gform_run_fingerprint(
    output_stub = output_stub,
    boot_iter = boot_iter,
    boot_seed = boot_seed,
    total_births = total_births,
    sample_frac = sample_frac
  )
  isTRUE(gform_bootstrap_ck_matches(ck, fp) && ck$last_completed >= boot_iter)
}

compute_bootstrap_ci <- function(weekly_boot, pop_boot) {
  weekly_boot <- data.table::as.data.table(weekly_boot)
  pop_boot <- data.table::as.data.table(pop_boot)

  weekly_ci <- weekly_boot[, .(
    risk_natural_lcl = stats::quantile(risk_natural, 0.025, na.rm = TRUE),
    risk_natural_ucl = stats::quantile(risk_natural, 0.975, na.rm = TRUE),
    risk_intervention_lcl = stats::quantile(risk_intervention, 0.025, na.rm = TRUE),
    risk_intervention_ucl = stats::quantile(risk_intervention, 0.975, na.rm = TRUE),
    risk_ratio_lcl = stats::quantile(risk_ratio, 0.025, na.rm = TRUE),
    risk_ratio_ucl = stats::quantile(risk_ratio, 0.975, na.rm = TRUE),
    risk_difference_lcl = stats::quantile(risk_difference, 0.025, na.rm = TRUE),
    risk_difference_ucl = stats::quantile(risk_difference, 0.975, na.rm = TRUE)
  ), by = "week"]

  pop_ci <- pop_boot[, {
    out <- list(
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
    )
    if ("attributable_fraction" %in% names(pop_boot)) {
      out$attributable_fraction_lcl <- stats::quantile(
        attributable_fraction, 0.025, na.rm = TRUE
      )
      out$attributable_fraction_ucl <- stats::quantile(
        attributable_fraction, 0.975, na.rm = TRUE
      )
    }
    out
  }, by = "scenario"]

  if (!"attributable_fraction" %in% names(pop_boot)) {
    paf_by_iter <- pop_boot[, {
      nat_prev <- prevalence[scenario == "observed"][1L]
      int_prev <- prevalence[scenario == "intervention"][1L]
      paf_val <- if (is.finite(nat_prev) && nat_prev > 0) {
        (nat_prev - int_prev) / nat_prev
      } else {
        NA_real_
      }
      .(attributable_fraction = paf_val)
    }, by = "iter"]

    paf_lcl <- stats::quantile(paf_by_iter$attributable_fraction, 0.025, na.rm = TRUE)
    paf_ucl <- stats::quantile(paf_by_iter$attributable_fraction, 0.975, na.rm = TRUE)

    pop_ci[, attributable_fraction_lcl := NA_real_]
    pop_ci[, attributable_fraction_ucl := NA_real_]
    pop_ci[scenario == "observed", `:=`(
      attributable_fraction_lcl = 0,
      attributable_fraction_ucl = 0
    )]
    pop_ci[scenario == "intervention", `:=`(
      attributable_fraction_lcl = as.numeric(paf_lcl),
      attributable_fraction_ucl = as.numeric(paf_ucl)
    )]
  }

  list(
    weekly_ci = tibble::as_tibble(weekly_ci),
    population_ci = tibble::as_tibble(pop_ci)
  )
}

bootstrap_gformula_effects <- function(
    model_store,
    person_weeks,
    cox_frame_natural,
    cox_frame_intervention,
    risk_weeks_vec,
    total_births,
    output_stub,
    dir_bootstrap,
    dependent_var = GFORM_DEFAULTS$dependent_var,
    risk_entry_week = GFORM_DEFAULTS$risk_entry_week,
    boot_iter = GFORM_DEFAULTS$boot_iter,
    boot_seed = GFORM_DEFAULTS$boot_seed,
    target_week = GFORM_DEFAULTS$population_week,
    baseline_scenario = GFORM_DEFAULTS$baseline_scenario,
    intervention_scenario = "intervention",
    sample_frac = NULL,
    resume = TRUE,
    checkpoint_every = 10L,
    parallel = FALSE,
    bootstrap_batch_size = NULL) {

  paths <- gform_bootstrap_paths(output_stub, dir_bootstrap)
  dir.create(paths$dir, recursive = TRUE, showWarnings = FALSE)
  fingerprint <- gform_run_fingerprint(
    output_stub = output_stub,
    boot_iter = boot_iter,
    boot_seed = boot_seed,
    total_births = total_births,
    sample_frac = sample_frac
  )

  append_one_result <- function(one, iter) {
    wdt <- data.table::as.data.table(one$weekly)
    wdt[, iter := iter]
    pdt <- data.table::as.data.table(one$population)
    pdt[, iter := iter]

    append_mode <- file.exists(paths$weekly) && iter > 1L
    data.table::fwrite(
      wdt, paths$weekly,
      append = append_mode,
      col.names = !append_mode
    )
    data.table::fwrite(
      pdt, paths$population,
      append = append_mode,
      col.names = !append_mode
    )
  }

  run_one_iter <- function(iter) {
    bootstrap_one_iter(
      iter = iter,
      model_store = model_store,
      person_weeks = person_weeks,
      cox_frame_natural = cox_frame_natural,
      cox_frame_intervention = cox_frame_intervention,
      risk_weeks_vec = risk_weeks_vec,
      total_births = total_births,
      dependent_var = dependent_var,
      risk_entry_week = risk_entry_week,
      target_week = target_week,
      baseline_scenario = baseline_scenario,
      intervention_scenario = intervention_scenario,
      boot_seed = boot_seed
    )
  }

  start_iter <- 1L
  if (resume && file.exists(paths$checkpoint)) {
    ck <- readRDS(paths$checkpoint)
    if (gform_bootstrap_ck_matches(ck, fingerprint) &&
        ck$last_completed >= boot_iter) {
      message("Bootstrap ya completo (", boot_iter, " iter); leyendo réplicas desde disco.")
      weekly_boot <- data.table::fread(paths$weekly)
      pop_boot <- data.table::fread(paths$population)
      ci <- compute_bootstrap_ci(weekly_boot, pop_boot)
      return(list(
        weekly_boot = weekly_boot,
        population_boot = pop_boot,
        bootstrap_paths = paths,
        weekly_ci = ci$weekly_ci,
        population_ci = ci$population_ci
      ))
    }
    if (gform_bootstrap_ck_matches(ck, fingerprint)) {
      start_iter <- ck$last_completed + 1L
      message("Bootstrap: reanudando desde iteración ", start_iter, " / ", boot_iter)
    } else {
      message("Checkpoint bootstrap incompatible; reiniciando bootstrap desde 1.")
    }
  }

  if (start_iter == 1L) {
    if (file.exists(paths$weekly)) file.remove(paths$weekly)
    if (file.exists(paths$population)) file.remove(paths$population)
    if (file.exists(paths$checkpoint)) file.remove(paths$checkpoint)
  }

  cfg <- getOption("gform.parallel", gform_parallel_config())
  if (is.null(bootstrap_batch_size)) {
    bootstrap_batch_size <- if (parallel) cfg$bootstrap_batch_size else 1L
  }
  bootstrap_batch_size <- max(1L, as.integer(bootstrap_batch_size))

  if (parallel && requireNamespace("furrr", quietly = TRUE)) {
    gform_setup_parallel(task = "bootstrap", config = cfg)
    globals_limit_gb <- getOption("future.globals.maxSize", 0) / 1024^3
    message(
      "Bootstrap paralelo: lotes de ", bootstrap_batch_size,
      " iteraciones (", cfg$n_workers_bootstrap, " workers)",
      " | globals max: ", round(globals_limit_gb, 1), " GiB",
      if (isTRUE(cfg$use_fork)) " | fork/multicore" else ""
    )
    batch_starts <- seq(start_iter, boot_iter, by = bootstrap_batch_size)
    for (batch_start in batch_starts) {
      batch_end <- min(batch_start + bootstrap_batch_size - 1L, boot_iter)
      iters <- batch_start:batch_end
      batch_results <- furrr::future_map(
        iters,
        run_one_iter,
        .options = gform_furrr_options()
      )
      for (k in seq_along(iters)) {
        append_one_result(batch_results[[k]], iters[k])
      }
      rm(batch_results)
      gc(verbose = FALSE)

      gform_save_bootstrap_checkpoint(paths$checkpoint, fingerprint, batch_end)
      message("Bootstrap checkpoint: ", batch_end, " / ", boot_iter)
    }
    future::plan(future::sequential)
  } else {
    if (isTRUE(parallel)) {
      message("Bootstrap paralelo no disponible (furrr); ejecutando secuencial.")
    } else {
      message("Bootstrap secuencial: ", boot_iter, " iteraciones (checkpoint cada ",
              checkpoint_every, ")")
    }
    for (iter in start_iter:boot_iter) {
      one <- run_one_iter(iter)
      append_one_result(one, iter)
      rm(one)
      gc(verbose = FALSE)

      if (iter %% checkpoint_every == 0L || iter == boot_iter) {
        gform_save_bootstrap_checkpoint(paths$checkpoint, fingerprint, iter)
        message("Bootstrap checkpoint: ", iter, " / ", boot_iter)
      }
    }
  }

  message("Bootstrap: leyendo réplicas y calculando IC...")
  weekly_boot <- data.table::fread(paths$weekly)
  pop_boot <- data.table::fread(paths$population)
  ci <- compute_bootstrap_ci(weekly_boot, pop_boot)

  list(
    weekly_boot = weekly_boot,
    population_boot = pop_boot,
    weekly_ci = ci$weekly_ci,
    population_ci = ci$population_ci,
    bootstrap_paths = paths
  )
}

## Orquestador por intervención ----

run_gform_intervention <- function(
    intervention_spec,
    intervention_path,
    data_base,
    wide_tad_obs,
    model_store,
    person_weeks,
    risk_weeks_vec,
    control_vars,
    total_births,
    cox_frame_natural = NULL,
    data_long = NULL,
    wide_exposicion_natural = NULL,
    raw_wide_pollutant = NULL,
    dependent_var = GFORM_DEFAULTS$dependent_var,
    risk_entry_week = GFORM_DEFAULTS$risk_entry_week,
    follow_up_weeks = GFORM_DEFAULTS$follow_up_weeks,
    boot_iter = GFORM_DEFAULTS$boot_iter,
    boot_seed = GFORM_DEFAULTS$boot_seed,
    target_week = GFORM_DEFAULTS$population_week,
    sample_frac = NULL,
    run_bootstrap = TRUE,
    run_singleweek_heatmap = FALSE,
    parallel_singleweek_heatmap = FALSE,
    parallel_bootstrap = FALSE,
    dir_bootstrap = NULL,
    bootstrap_resume = TRUE,
    dir_heatmap = NULL,
    heatmap_resume = TRUE) {

  timing <- gform_timing_log_init()
  fingerprint <- gform_run_fingerprint(
    output_stub = intervention_spec$output_stub,
    boot_iter = boot_iter,
    boot_seed = boot_seed,
    total_births = total_births,
    sample_frac = sample_frac
  )
  if (run_bootstrap && boot_iter > 0L && !is.null(dir_bootstrap) &&
      gform_bootstrap_is_complete(
        output_stub = intervention_spec$output_stub,
        dir_bootstrap = dir_bootstrap,
        boot_iter = boot_iter,
        boot_seed = boot_seed,
        total_births = total_births,
        sample_frac = sample_frac
      )) {
    message("Bootstrap ya completo (", boot_iter, " iter); lectura desde disco al llegar a bootstrap.")
  }
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
      message(
        "Reanudación: punto estimado en disco (",
        format(pt_ck$saved_at, "%Y-%m-%d %H:%M"), "); omitiendo predicciones."
      )
    }
  }

  block <- gform_time_block("Construir frame Cox intervención", {
    intervention_obj <- readRDS(intervention_path)
    pollutant <- intervention_spec$pollutant

    if (is.null(raw_wide_pollutant)) {
      if (is.null(data_long)) {
        stop("Se requiere data_long o raw_wide_pollutant preconstruido.")
      }
      raw_wide_pollutant <- build_wide_raw_exposure(
        data_long, pollutant, weeks_keep = GFORM_DEFAULTS$weeks_exposure
      )
    }
    if (is.null(wide_exposicion_natural) && is.null(cox_frame_natural)) {
      wide_exposicion_natural <- build_exposicion_wide_from_raw(
        raw_wide = raw_wide_pollutant,
        pollutant = pollutant,
        intervention = list(type = "none")
      )
    }
    if (is.null(cox_frame_natural)) {
      cox_frame_natural <- build_cox_model_frame(
        data_base = data_base,
        wide_exposicion = wide_exposicion_natural,
        wide_tad_obs = wide_tad_obs,
        control_vars = control_vars,
        dependent_var = dependent_var,
        risk_entry_week = risk_entry_week
      )
    }

    wide_exposicion_intervention <- intervention_obj$wide_exposicion
    wide_exposicion_intervention <- wide_exposicion_intervention[id %in% data_base$id]
    cox_frame_intervention <- build_cox_model_frame(
      data_base = data_base,
      wide_exposicion = wide_exposicion_intervention,
      wide_tad_obs = wide_tad_obs,
      control_vars = control_vars,
      dependent_var = dependent_var,
      risk_entry_week = risk_entry_week
    )
    rm(intervention_obj)
    gc()

    list(
      pollutant = pollutant,
      raw_wide_pollutant = raw_wide_pollutant,
      cox_frame_natural = cox_frame_natural,
      cox_frame_intervention = cox_frame_intervention
    )
  })
  timing <- gform_timing_log_add(timing, block$timing)
  pollutant <- block$result$pollutant
  raw_wide_pollutant <- block$result$raw_wide_pollutant
  cox_frame_natural <- block$result$cox_frame_natural
  cox_frame_intervention <- block$result$cox_frame_intervention

  if (isTRUE(skip_point_estimate)) {
    pt_ck <- gform_read_point_checkpoint(point_ck_path, fingerprint)
    nat_mean <- pt_ck$nat_mean
    weekly_effects <- pt_ck$weekly_effects
    cumulative_risk_curves <- pt_ck$cumulative_risk_curves
    population_effects <- pt_ck$population_effects
    timing <- gform_timing_log_add(timing, list(
      label = "Punto estimado (checkpoint)",
      start = pt_ck$saved_at,
      end = pt_ck$saved_at,
      sec = 0
    ))
  } else {
  block <- gform_time_block("Predicción hazards — curso natural", {
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

  block <- gform_time_block("Predicción hazards — intervención", {
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

  block <- gform_time_block("Efectos puntuales (supervivencia, semanal, poblacional, curvas)", {
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
      stop("run_gform_intervention: se requiere dir_bootstrap para bootstrap secuencial.")
    }
    block <- gform_time_block(paste0("Bootstrap paramétrico (", boot_iter, " iteraciones)"), {
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
      weekly_ci <- dplyr::left_join(weekly_effects, boot_out$weekly_ci, by = "week")
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

  singleweek_intervention_heatmap <- NULL
  if (run_singleweek_heatmap) {
    if (is.null(dir_heatmap)) {
      stop("run_gform_intervention: se requiere dir_heatmap para mapa de calor.")
    }
    block <- gform_time_block(
      paste0("Mapa calor semana única (", length(GFORM_DEFAULTS$weeks_exposure), " columnas)"),
      {
        compute_singleweek_intervention_heatmap(
          model_store = model_store,
          person_weeks = person_weeks,
          data_base = data_base,
          raw_wide_pollutant = raw_wide_pollutant,
          pollutant = pollutant,
          wide_tad_obs = wide_tad_obs,
          risk_weeks_vec = risk_weeks_vec,
          control_vars = control_vars,
          nat_mean = nat_mean,
          follow_up_weeks = follow_up_weeks,
          output_stub = intervention_spec$output_stub,
          dir_heatmap = dir_heatmap,
          total_births = total_births,
          sample_frac = sample_frac,
          resume = heatmap_resume,
          parallel = parallel_singleweek_heatmap
        )
      }
    )
    timing <- gform_timing_log_add(timing, block$timing)
    singleweek_intervention_heatmap <- block$result
  }

  list(
    intervention_spec = intervention_spec,
    weekly_effects = weekly_ci,
    population_effects = population_ci,
    cumulative_risk_curves = cumulative_risk_curves,
    singleweek_intervention_heatmap = singleweek_intervention_heatmap,
    figure3 = cumulative_risk_curves,
    figure4 = singleweek_intervention_heatmap,
    bootstrap = boot_out,
    timing = timing
  )
}

gform_intervention_is_complete <- function(
    output_stub,
    dir_weekly,
    dir_population,
    dir_summary,
    dir_bootstrap,
    boot_iter = GFORM_DEFAULTS$boot_iter,
    boot_seed = GFORM_DEFAULTS$boot_seed,
    total_births = NULL,
    sample_frac = NULL) {

  weekly_path <- file.path(dir_weekly, paste0(output_stub, "_weekly_effects.rds"))
  population_path <- file.path(dir_population, paste0(output_stub, "_population_effects.rds"))
  excel_path <- file.path(dir_summary, paste0(output_stub, "_point_estimates.xlsx"))
  boot_ck <- gform_bootstrap_paths(output_stub, dir_bootstrap)$checkpoint

  if (!all(file.exists(c(weekly_path, population_path, excel_path)))) {
    return(FALSE)
  }
  if (boot_iter > 0L && file.exists(boot_ck)) {
    ck <- readRDS(boot_ck)
    if (!is.null(total_births)) {
      fp <- gform_run_fingerprint(
        output_stub, boot_iter, boot_seed, total_births, sample_frac
      )
      if (!gform_bootstrap_ck_matches(ck, fp)) return(FALSE)
    }
    return(isTRUE(ck$last_completed >= boot_iter))
  }
  boot_iter <= 0L
}

gform_load_email_config <- function(secrets_dir = ".secrets") {
  env_path <- file.path(secrets_dir, "gmail.env")
  cfg <- list(
    from = Sys.getenv("GFORM_GMAIL_FROM", "jo.conejeros@gmail.com"),
    to = Sys.getenv("GFORM_GMAIL_TO", "jdconejeros@uc.cl"),
    app_password = Sys.getenv("GFORM_GMAIL_APP_PASSWORD", "")
  )
  if (file.exists(env_path)) {
    lines <- readLines(env_path, warn = FALSE)
    for (ln in lines) {
      if (!nzchar(ln) || grepl("^\\s*#", ln)) next
      parts <- strsplit(ln, "=", fixed = TRUE)[[1L]]
      if (length(parts) < 2L) next
      key <- trimws(parts[[1L]])
      val <- trimws(paste(parts[-1L], collapse = "="))
      switch(key,
        GFORM_GMAIL_FROM = cfg$from <- val,
        GFORM_GMAIL_TO = cfg$to <- val,
        GFORM_GMAIL_APP_PASSWORD = cfg$app_password <- val
      )
    }
  }
  cfg
}

gform_send_email_summary <- function(
    subject,
    body,
    secrets_dir = ".secrets",
    log_path = NULL) {

  cfg <- gform_load_email_config(secrets_dir)
  sent <- FALSE
  err_msg <- character()

  if (nzchar(cfg$app_password)) {
    sent <- tryCatch({
      body_file <- tempfile(fileext = ".txt")
      writeLines(body, body_file, useBytes = TRUE)
      py_code <- sprintf(
        paste(
          "import smtplib, ssl",
          "from email.message import EmailMessage",
          "from pathlib import Path",
          "msg = EmailMessage()",
          "msg['Subject'] = %s",
          "msg['From'] = %s",
          "msg['To'] = %s",
          "msg.set_content(Path(%s).read_text(encoding='utf-8'))",
          "ctx = ssl.create_default_context()",
          "with smtplib.SMTP_SSL('smtp.gmail.com', 465, context=ctx) as s:",
          "    s.login(%s, %s)",
          "    s.send_message(msg)",
          sep = "\n"
        ),
        shQuote(subject, type = "cmd"),
        shQuote(cfg$from, type = "cmd"),
        shQuote(cfg$to, type = "cmd"),
        shQuote(body_file, type = "cmd"),
        shQuote(cfg$from, type = "cmd"),
        shQuote(cfg$app_password, type = "cmd")
      )
      status <- system2("python3", c("-c", py_code), stdout = FALSE, stderr = FALSE)
      unlink(body_file)
      identical(status, 0L)
    }, error = function(e) {
      err_msg <<- c(err_msg, conditionMessage(e))
      FALSE
    })
  }

  if (!sent && .Platform$OS.type == "unix") {
    esc <- function(x) gsub("\\\\", "\\\\\\\\", gsub('"', '\\"', x))
    script <- sprintf(
      paste(
        "tell application \"Mail\"",
        "  set msg to make new outgoing message with properties {subject:\"%s\", content:\"%s\", visible:false}",
        "  tell msg",
        "    make new to recipient at end of to recipients with properties {address:\"%s\"}",
        "    send",
        "  end tell",
        "end tell",
        sep = "\n"
      ),
      esc(subject), esc(body), esc(cfg$to)
    )
    sent <- tryCatch({
      status <- system2("osascript", c("-e", script), stdout = FALSE, stderr = FALSE)
      identical(status, 0L)
    }, error = function(e) {
      err_msg <<- c(err_msg, conditionMessage(e))
      FALSE
    })
  }

  if (is.null(log_path)) {
    log_path <- file.path("02_Output/G-Form/Timing", "email_log.txt")
  }
  dir.create(dirname(log_path), recursive = TRUE, showWarnings = FALSE)
  stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- paste0(
    "[", stamp, "] sent=", sent, " | to=", cfg$to,
    if (length(err_msg)) paste0(" | err=", paste(err_msg, collapse = " ; ")) else "",
    "\n", body, "\n", strrep("-", 72), "\n"
  )
  cat(line, file = log_path, append = TRUE)

  if (!sent) {
    message("Email no enviado (ver ", log_path, "). Configure GFORM_GMAIL_APP_PASSWORD en .secrets/gmail.env")
  } else {
    message("Email de resumen enviado a ", cfg$to)
  }
  invisible(sent)
}

gform_format_intervention_email <- function(result) {
  tl <- result$timing_log
  if (length(tl$phases)) {
    timing_lines <- vapply(tl$phases, function(p) {
      sprintf("  %-42s %8.1f s (%.1f min)", p$label, p$sec, p$sec / 60)
    }, character(1))
    total_sec <- tl$total_sec
  } else {
    timing_lines <- vapply(names(tl), function(nm) {
      if (nm %in% c(
        "phases", "started_at", "n_births", "n_original", "sample_frac",
        "intervention_number", "intervention_id", "total_sec",
        "run_parallel_cox", "run_parallel_bootstrap", "run_singleweek_heatmap",
        "run_parallel_singleweek_heatmap", "n_workers_bootstrap"
      )) {
        return("")
      }
      if (is.numeric(tl[[nm]])) {
        return(sprintf("  %-42s %8.1f s (%.1f min)", nm, tl[[nm]], tl[[nm]] / 60))
      }
      ""
    }, character(1))
    timing_lines <- timing_lines[nzchar(timing_lines)]
    total_sec <- if (is.numeric(tl$total_sec)) tl$total_sec else NA_real_
  }

  paste(
    "G-Formula — intervención completada",
    "",
    sprintf("Intervención %d: %s", result$intervention_number, result$intervention_id),
    sprintf("Descripción: %s", result$description),
    "",
    sprintf("Prevalencia natural (sem. 36): %.6f", result$prevalence_natural),
    sprintf("Prevalencia intervención:      %.6f", result$prevalence_intervention),
    sprintf("Risk difference:               %.6f", result$risk_difference),
    "",
    "Tiempos:",
    paste(timing_lines, collapse = "\n"),
    "",
    sprintf("Tiempo total: %.1f min", total_sec / 60),
    "",
    "Outputs:",
    paste(" ", result$output_files, collapse = "\n"),
    sep = "\n"
  )
}

## Guardado de resultados ----

save_results <- function(
    weekly_effects,
    population_effects,
    weekly_path,
    population_path,
    cumulative_risk_curves = NULL,
    cumulative_risk_curves_path = NULL,
    singleweek_intervention_heatmap = NULL,
    singleweek_intervention_heatmap_path = NULL,
    figure3 = NULL,
    figure3_path = NULL,
    figure4 = NULL,
    figure4_path = NULL,
    weekly_boot = NULL,
    population_boot = NULL,
    metadata = list()) {

  cumulative_risk_curves <- cumulative_risk_curves %||% figure3
  cumulative_risk_curves_path <- cumulative_risk_curves_path %||% figure3_path
  singleweek_intervention_heatmap <- singleweek_intervention_heatmap %||% figure4
  singleweek_intervention_heatmap_path <- singleweek_intervention_heatmap_path %||% figure4_path

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

  if (!is.null(cumulative_risk_curves) && !is.null(cumulative_risk_curves_path)) {
    dir.create(dirname(cumulative_risk_curves_path), recursive = TRUE, showWarnings = FALSE)
    desc <- attr(cumulative_risk_curves, "description")
    saveRDS(c(
      list(
        point = cumulative_risk_curves,
        description = desc
      ),
      metadata
    ), cumulative_risk_curves_path)
  }
  if (!is.null(singleweek_intervention_heatmap) && !is.null(singleweek_intervention_heatmap_path)) {
    dir.create(dirname(singleweek_intervention_heatmap_path), recursive = TRUE, showWarnings = FALSE)
    desc <- attr(singleweek_intervention_heatmap, "description")
    saveRDS(c(
      list(
        point = singleweek_intervention_heatmap,
        description = desc
      ),
      metadata
    ), singleweek_intervention_heatmap_path)
  }

  invisible(list(weekly = weekly_path, population = population_path))
}

save_gform_excel <- function(
    results,
    excel_path,
    intervention_id) {

  curves <- results$cumulative_risk_curves %||% results$figure3
  heatmap <- results$singleweek_intervention_heatmap %||% results$figure4

  sheets <- list(
    weekly_effects = results$weekly_effects,
    population_effects = results$population_effects,
    curvas_riesgo_acumulado_global = curves
  )
  if (!is.null(heatmap)) {
    sheets$mapa_calor_rd_semana_intervencion_long <- heatmap$long
    sheets$mapa_calor_rd_semana_intervencion_wide <- heatmap$wide
  }

  dir.create(dirname(excel_path), recursive = TRUE, showWarnings = FALSE)
  writexl::write_xlsx(sheets, path = excel_path)
  invisible(excel_path)
}

gform_furrr_options <- function() {
  opts <- list(seed = TRUE)
  if ("globals.onReference" %in% names(formals(furrr::furrr_options))) {
    opts$globals.onReference <- "ignore"
  }
  do.call(furrr::furrr_options, opts)
}

gform_finalize_run <- function() {
  if (requireNamespace("future", quietly = TRUE)) {
    future::plan(future::sequential)
  }
  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::setDTthreads(1L)
  }
  gc(verbose = FALSE)
  invisible(NULL)
}
