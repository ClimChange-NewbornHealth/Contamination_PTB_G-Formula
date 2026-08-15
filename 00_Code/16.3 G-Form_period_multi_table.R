# 16.3 G-Formula por períodos multicontaminante — tabla resumen ----
#
# Uso (desde la raíz del proyecto):
#   Rscript "00_Code/16.3 G-Form_period_multi_table.R"
#
# Entrada:
#   02_Output/G-Form-Period-Multi/Summary_results/{stub}_point_estimates.xlsx
# Salida:
#   02_Output/G-Form-Period-Multi/Summary_results/Table_population_effects_summary.xlsx

source("00_Code/0.1 Settings.R")

install_load <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
    suppressPackageStartupMessages(
      library(pkg, character.only = TRUE)
    )
  }
}

install_load(c("readxl", "dplyr", "openxlsx", "writexl"))
source("00_Code/16.0 G-Form_period_multi_functions.R")

data_out_g <- GFORM_PERIOD_MULTI_DATA_OUT
dir_summary <- file.path(data_out_g, "Summary_results")
dir_bootstrap <- file.path(data_out_g, "Bootstrap")
dir_population <- file.path(data_out_g, "PopulationEffects")
path_output <- file.path(dir_summary, "Table_population_effects_summary.xlsx")

time_intervention_order <- c("Trimester", "Overall")

pollutant_specs <- list(
  pm25 = list(
    stub_prefix = "pm25",
    section_label = "PM2.5 (\u00b5g/m\u00b3) [multicontaminante]",
    scenarios = list(
      list(stub = "pm25_pct20", label = "Reduced by 20%")
    )
  ),
  no2 = list(
    stub_prefix = "no2",
    section_label = "NO2 (ppbv) [multicontaminante]",
    scenarios = list(
      list(stub = "no2_pct20", label = "Reduced by 20%")
    )
  ),
  o3 = list(
    stub_prefix = "o3",
    section_label = "O3 (ppbv) [multicontaminante]",
    scenarios = list(
      list(stub = "o3_pct20", label = "Reduced by 20%")
    )
  )
)

table_percent_scale <- 100

table_columns <- c(
  "Exposure and scenario",
  "Time Intervention",
  "Prevalence (95% CI, %)",
  "Cases (95% CI)",
  "Risk Ratio (95% CI)",
  "Risk Difference (95% CI, pp)",
  "Attributable Risk (95% CI, pp)",
  "Population Attributable Fraction (95% CI, %)"
)

metric_specs <- list(
  prevalence = c("prevalence", "prevalence_lcl", "prevalence_ucl"),
  cases = c("cases", "cases_lcl", "cases_ucl"),
  risk_ratio = c("risk_ratio", "risk_ratio_lcl", "risk_ratio_ucl"),
  risk_difference = c("risk_difference", "risk_difference_lcl", "risk_difference_ucl"),
  attributable_risk = c("attributable_risk", "attributable_risk_lcl", "attributable_risk_ucl"),
  attributable_fraction = c(
    "attributable_fraction",
    "attributable_fraction_lcl",
    "attributable_fraction_ucl"
  )
)

format_gform_table_num <- function(x) {
  if (length(x) != 1L || is.na(x) || !is.finite(x)) return(NA_character_)
  if (abs(x) > 2) {
    return(format(round(x, 0), trim = TRUE, decimal.mark = ".", scientific = FALSE))
  }
  formatC(x, format = "f", digits = 4, decimal.mark = ".")
}

format_gform_table_rr_num <- function(x) {
  if (length(x) != 1L || is.na(x) || !is.finite(x)) return(NA_character_)
  formatC(x, format = "f", digits = 3, decimal.mark = ".")
}

format_gform_table_percent_num <- function(x) {
  if (length(x) != 1L || is.na(x) || !is.finite(x)) return(NA_character_)
  formatC(x, format = "f", digits = 2, decimal.mark = ".")
}

format_estimate_ci <- function(
    est,
    lcl,
    ucl,
    percent_display = FALSE,
    num_fmt = NULL) {
  fmt <- if (!is.null(num_fmt)) {
    num_fmt
  } else if (percent_display) {
    format_gform_table_percent_num
  } else {
    format_gform_table_num
  }
  est_chr <- fmt(est)
  lcl_chr <- fmt(lcl)
  ucl_chr <- fmt(ucl)
  if (any(is.na(c(est_chr, lcl_chr, ucl_chr)))) return(NA_character_)
  paste0(est_chr, " (", lcl_chr, "; ", ucl_chr, ")")
}

format_metric_row <- function(row, metric_name, scale = 1) {
  cols <- metric_specs[[metric_name]]
  num_fmt <- if (metric_name == "risk_ratio") format_gform_table_rr_num else NULL
  format_estimate_ci(
    row[[cols[[1]]]] * scale,
    row[[cols[[2]]]] * scale,
    row[[cols[[3]]]] * scale,
    percent_display = scale != 1,
    num_fmt = num_fmt
  )
}

has_dplyr_join_suffix <- function(df) {
  any(grepl("\\.(x|y)$", names(df), perl = TRUE))
}

needs_effects_sheet_repair <- function(df) {
  nms <- names(df)
  has_dplyr_join_suffix(df) ||
    any(grepl("\\.\\.\\.", nms)) ||
    any(grepl("\\.x$|\\.y$", nms, perl = TRUE)) ||
    (sum(nms == "attributable_fraction_lcl") > 1L) ||
    (sum(nms == "attributable_fraction_ucl") > 1L)
}

repair_dplyr_join_suffix <- function(df) {
  df <- as.data.frame(df)
  nms <- names(df)
  if (any(grepl("\\.\\.\\.", nms))) {
    nms <- sub("\\.\\.\\.\\d+$", "", nms)
    names(df) <- nms
  }
  if (any(grepl("\\.y$", nms, perl = TRUE))) {
    df <- df[, !grepl("\\.y$", names(df), perl = TRUE), drop = FALSE]
  }
  if (any(grepl("\\.x$", names(df), perl = TRUE))) {
    names(df) <- sub("\\.x$", "", names(df), perl = TRUE)
  }
  dedupe_df_columns(df)
}

dedupe_df_columns <- function(df) {
  df <- as.data.frame(df)
  if (!ncol(df)) return(df)
  nms <- names(df)
  keep <- !duplicated(nms)
  df[, keep, drop = FALSE]
}

compute_paf_ci_from_boot <- function(boot_path) {
  if (!file.exists(boot_path)) return(NULL)
  boot <- utils::read.csv(boot_path, stringsAsFactors = FALSE)
  if (!all(c("scenario", "prevalence", "iter") %in% names(boot))) return(NULL)

  paf_by_iter <- lapply(split(boot, boot$iter), function(df) {
    nat_prev <- df$prevalence[df$scenario == "observed"][1L]
    int_prev <- df$prevalence[df$scenario == "intervention"][1L]
    if (!is.finite(nat_prev) || nat_prev <= 0) return(NA_real_)
    (nat_prev - int_prev) / nat_prev
  })

  paf_vals <- unlist(paf_by_iter, use.names = FALSE)
  paf_vals <- paf_vals[is.finite(paf_vals)]
  if (!length(paf_vals)) return(NULL)

  c(
    attributable_fraction_lcl = as.numeric(stats::quantile(paf_vals, 0.025, na.rm = TRUE)),
    attributable_fraction_ucl = as.numeric(stats::quantile(paf_vals, 0.975, na.rm = TRUE))
  )
}

augment_population_effects <- function(pop) {
  pop <- dplyr::as_tibble(pop)
  baseline_prev <- pop$prevalence[pop$scenario == "observed"][1L]
  if (!is.finite(baseline_prev) || baseline_prev <= 0) {
    stop("No se pudo determinar la prevalencia del curso natural.")
  }
  if (!"attributable_fraction" %in% names(pop)) {
    pop$attributable_fraction <- ifelse(
      pop$scenario == "observed",
      0,
      (baseline_prev - pop$prevalence) / baseline_prev
    )
  }
  if (!"attributable_fraction_lcl" %in% names(pop)) {
    pop$attributable_fraction_lcl <- ifelse(pop$scenario == "observed", 0, NA_real_)
  }
  if (!"attributable_fraction_ucl" %in% names(pop)) {
    pop$attributable_fraction_ucl <- ifelse(pop$scenario == "observed", 0, NA_real_)
  }
  pop
}

read_population_effects <- function(excel_path) {
  sheets <- readxl::excel_sheets(excel_path)
  hit <- sheets[grepl("^population_effects$", sheets, ignore.case = TRUE)]
  if (!length(hit)) stop("No se encontró hoja 'population_effects' en ", excel_path)
  readxl::read_excel(excel_path, sheet = hit[[1L]])
}

read_scenario_metadata <- function(excel_path) {
  sheets <- readxl::excel_sheets(excel_path)
  if (!"metadata" %in% sheets) return(list(application_scope = NA_character_))
  meta <- readxl::read_excel(excel_path, sheet = "metadata")
  if (!all(c("field", "value") %in% names(meta))) {
    return(list(application_scope = NA_character_))
  }
  scope <- meta$value[meta$field == "application_scope"][1L]
  list(application_scope = as.character(scope))
}

resolve_time_intervention <- function(application_scope, stub) {
  if (identical(application_scope, "trimester") || grepl("_tri_multi$", stub)) {
    return("Trimester")
  }
  if (identical(application_scope, "full") || grepl("_full_multi$", stub)) {
    return("Overall")
  }
  NA_character_
}

load_scenario_population <- function(output_stub, target_scenario = c("observed", "intervention")) {
  target_scenario <- match.arg(target_scenario)
  excel_path <- file.path(dir_summary, paste0(output_stub, "_point_estimates.xlsx"))
  if (!file.exists(excel_path)) return(NULL)

  pop <- read_population_effects(excel_path) |> augment_population_effects()
  boot_path <- file.path(dir_bootstrap, output_stub, "population_boot.csv")
  paf_ci <- compute_paf_ci_from_boot(boot_path)

  if (!is.null(paf_ci) && target_scenario == "intervention") {
    pop$attributable_fraction_lcl[pop$scenario == "intervention"] <- paf_ci[["attributable_fraction_lcl"]]
    pop$attributable_fraction_ucl[pop$scenario == "intervention"] <- paf_ci[["attributable_fraction_ucl"]]
  }

  row <- pop |> dplyr::filter(.data$scenario == target_scenario)
  if (!nrow(row)) stop("No se encontró escenario '", target_scenario, "' en ", excel_path)
  as.list(row[1, , drop = FALSE])
}

format_table_row <- function(scenario_label, pop_row, time_intervention = "") {
  pct <- table_percent_scale
  c(
    scenario_label,
    time_intervention,
    format_metric_row(pop_row, "prevalence", pct),
    format_metric_row(pop_row, "cases"),
    format_metric_row(pop_row, "risk_ratio"),
    format_metric_row(pop_row, "risk_difference", pct),
    format_metric_row(pop_row, "attributable_risk", pct),
    format_metric_row(pop_row, "attributable_fraction", pct)
  )
}

scenario_excel_path <- function(stub) {
  file.path(dir_summary, paste0(stub, "_point_estimates.xlsx")
}

period_multi_stub_candidates <- function(base_stub) {
  unique(c(
    paste0(base_stub, "_tri_multi"),
    paste0(base_stub, "_full_multi")
  ))
}

discover_available_scenarios <- function(scenarios) {
  available <- list()
  seen_stubs <- character()

  add_scenario <- function(stub, base_label, time_intervention) {
    if (stub %in% seen_stubs || !file.exists(scenario_excel_path(stub))) {
      return(invisible(NULL))
    }
    seen_stubs <<- c(seen_stubs, stub)
    available[[length(available) + 1L]] <<- list(
      stub = stub,
      label = base_label,
      time_intervention = time_intervention
    )
    invisible(NULL)
  }

  for (sc in scenarios) {
    tri_stub <- paste0(sc$stub, "_tri_multi")
    full_stub <- paste0(sc$stub, "_full_multi")
    if (file.exists(scenario_excel_path(tri_stub))) {
      add_scenario(tri_stub, sc$label, "Trimester")
    }
    if (file.exists(scenario_excel_path(full_stub))) {
      add_scenario(full_stub, sc$label, "Overall")
    }
  }

  available
}

build_pollutant_block <- function(spec) {
  available <- discover_available_scenarios(spec$scenarios)
  if (!length(available)) {
    warning("Sin resultados para ", spec$stub_prefix, "; se omite bloque.")
    return(NULL)
  }

  rows <- list(c(spec$section_label, rep("", length(metric_specs) + 1L)))

  for (time_int in time_intervention_order) {
    time_scenarios <- Filter(function(x) identical(x$time_intervention, time_int), available)
    if (!length(time_scenarios)) next

    natural_row <- load_scenario_population(time_scenarios[[1L]]$stub, "observed")
    if (is.null(natural_row)) next

    rows[[length(rows) + 1L]] <- format_table_row("Natural Course", natural_row, time_int)

    for (sc in time_scenarios) {
      intervention_row <- load_scenario_population(sc$stub, "intervention")
      rows[[length(rows) + 1L]] <- format_table_row(sc$label, intervention_row, time_int)
    }
  }

  rows
}

build_summary_table <- function() {
  blocks <- lapply(pollutant_specs, build_pollutant_block)
  blocks <- Filter(Negate(is.null), blocks)
  if (!length(blocks)) {
    stop("No hay resultados disponibles para construir la tabla.")
  }
  table_matrix <- do.call(rbind, unlist(blocks, recursive = FALSE))
  colnames(table_matrix) <- table_columns
  as.data.frame(table_matrix, stringsAsFactors = FALSE)
}

write_summary_workbook <- function(table_df, path_output) {
  n_cols <- ncol(table_df)
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "Table")

  header_style <- openxlsx::createStyle(
    textDecoration = "bold",
    halign = "center",
    valign = "center",
    border = "Bottom",
    wrapText = TRUE
  )
  section_style <- openxlsx::createStyle(
    textDecoration = "bold",
    fgFill = "#D9D9D9",
    border = c("top", "bottom"),
    halign = "left"
  )
  body_style <- openxlsx::createStyle(
    halign = "center",
    valign = "center",
    wrapText = TRUE
  )
  label_style <- openxlsx::createStyle(
    halign = "left",
    valign = "center"
  )

  openxlsx::writeData(wb, sheet = "Table", x = table_df, startRow = 1, colNames = TRUE)
  openxlsx::addStyle(
    wb, sheet = "Table", style = header_style,
    rows = 1, cols = seq_len(n_cols), gridExpand = TRUE
  )

  for (i in seq_len(nrow(table_df))) {
    is_section <- all(table_df[i, 2:n_cols, drop = TRUE] == "")
    openxlsx::addStyle(
      wb, sheet = "Table",
      style = if (is_section) section_style else body_style,
      rows = i + 1L, cols = 2:n_cols, gridExpand = TRUE, stack = TRUE
    )
    openxlsx::addStyle(
      wb, sheet = "Table",
      style = if (is_section) section_style else label_style,
      rows = i + 1L, cols = 1, gridExpand = TRUE, stack = TRUE
    )
    if (!is_section) {
      openxlsx::addStyle(
        wb, sheet = "Table", style = label_style,
        rows = i + 1L, cols = 2, gridExpand = TRUE, stack = TRUE
      )
    }
  }

  openxlsx::setColWidths(
    wb, sheet = "Table", cols = seq_len(n_cols),
    widths = c(28, 18, rep(24, n_cols - 2L))
  )

  dir.create(dirname(path_output), recursive = TRUE, showWarnings = FALSE)
  openxlsx::saveWorkbook(wb, path_output, overwrite = TRUE)
}

message("Construyendo tabla resumen de efectos poblacionales (período multi)...")
summary_table <- build_summary_table()
write_summary_workbook(summary_table, path_output)

message("Tabla guardada: ", path_output)
print(summary_table)
