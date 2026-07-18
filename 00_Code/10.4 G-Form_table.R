# 10.4 G-Formula — tabla resumen de efectos poblacionales ----
#
# Uso (desde la raíz del proyecto):
#   Rscript "00_Code/10.4 G-Form_table.R"
#
# Entrada:
#   02_Output/G-Form/Summary_results/{stub}_point_estimates.xlsx
#     stubs: pm25/no2/o3 _pct20, _lt20, _lt15, _lt10, _lt5 (según disponibilidad)
#   02_Output/G-Form/Bootstrap/{stub}/population_boot.csv  (IC de PAF si existe)
# Salida:
#   02_Output/G-Form/Summary_results/Table_population_effects_summary.xlsx
#
# Nota: augmenta métricas derivadas (PAF) en Excel existentes sin re-estimación.
# Repara columnas .x/.y (join bootstrap duplicado tras GFORM_HEATMAP_ONLY) sin re-estimar.
# La tabla resumen presenta prevalencia, RD, AR y PAF en escala % (×100);
# Cases y Risk Ratio permanecen en unidades originales (conteo y ratio).

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

## ===== Configuración =====
data_out_g <- "02_Output/G-Form/"
dir_summary <- file.path(data_out_g, "Summary_results")
dir_bootstrap <- file.path(data_out_g, "Bootstrap")
dir_weekly <- file.path(data_out_g, "WeeklyEffects")
dir_population <- file.path(data_out_g, "PopulationEffects")
path_output <- file.path(dir_summary, "Table_population_effects_summary.xlsx")

pollutant_specs <- list(
  pm25 = list(
    stub_prefix = "pm25",
    section_label = "PM2.5 (\u00b5g/m\u00b3)",
    scenarios = list(
      list(stub = "pm25_lt20", label = "< 20 \u00b5g/m\u00b3"),
      list(stub = "pm25_lt15", label = "< 15 \u00b5g/m\u00b3"),
      list(stub = "pm25_lt10", label = "< 10 \u00b5g/m\u00b3"),
      list(stub = "pm25_lt5", label = "< 5 \u00b5g/m\u00b3"),
      list(stub = "pm25_pct20", label = "Reduced by 20%")
    )
  ),
  no2 = list(
    stub_prefix = "no2",
    section_label = "NO2 (ppbv)",
    scenarios = list(
      list(stub = "no2_lt20", label = "< 20 ppbv"),
      list(stub = "no2_lt15", label = "< 15 ppbv"),
      list(stub = "no2_lt10", label = "< 10 ppbv"),
      list(stub = "no2_lt5", label = "< 5 ppbv"),
      list(stub = "no2_pct20", label = "Reduced by 20%")
    )
  ),
  o3 = list(
    stub_prefix = "o3",
    section_label = "O3 (ppbv)",
    scenarios = list(
      list(stub = "o3_pct20", label = "Reduced by 20%")
    )
  )
)

table_percent_scale <- 100

table_columns <- c(
  "Exposure and scenario",
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

## ===== Formato numérico =====
format_gform_table_num <- function(x) {
  if (length(x) != 1L || is.na(x) || !is.finite(x)) {
    return(NA_character_)
  }
  if (abs(x) > 2) {
    return(format(round(x, 0), trim = TRUE, decimal.mark = ".", scientific = FALSE))
  }
  formatC(x, format = "f", digits = 4, decimal.mark = ".")
}

format_gform_table_percent_num <- function(x) {
  if (length(x) != 1L || is.na(x) || !is.finite(x)) {
    return(NA_character_)
  }
  formatC(x, format = "f", digits = 2, decimal.mark = ".")
}

format_estimate_ci <- function(est, lcl, ucl, percent_display = FALSE) {
  fmt <- if (percent_display) format_gform_table_percent_num else format_gform_table_num
  est_chr <- fmt(est)
  lcl_chr <- fmt(lcl)
  ucl_chr <- fmt(ucl)
  if (any(is.na(c(est_chr, lcl_chr, ucl_chr)))) {
    return(NA_character_)
  }
  paste0(est_chr, " (", lcl_chr, "; ", ucl_chr, ")")
}

format_metric_row <- function(row, metric_name, scale = 1) {
  cols <- metric_specs[[metric_name]]
  format_estimate_ci(
    row[[cols[[1]]]] * scale,
    row[[cols[[2]]]] * scale,
    row[[cols[[3]]]] * scale,
    percent_display = scale != 1
  )
}

## ===== Reparar columnas .x/.y (join bootstrap duplicado) =====
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
  if (!ncol(df)) {
    return(df)
  }
  nms <- names(df)
  keep <- !duplicated(nms)
  df[, keep, drop = FALSE]
}

repair_point_estimates_stub <- function(output_stub) {
  repaired <- FALSE
  weekly_path <- file.path(dir_weekly, paste0(output_stub, "_weekly_effects.rds"))
  population_path <- file.path(
    dir_population, paste0(output_stub, "_population_effects.rds")
  )
  excel_path <- file.path(dir_summary, paste0(output_stub, "_point_estimates.xlsx"))

  if (file.exists(weekly_path)) {
    weekly_obj <- readRDS(weekly_path)
    fixed <- repair_dplyr_join_suffix(weekly_obj$point)
    if (!identical(names(fixed), names(weekly_obj$point))) {
      weekly_obj$point <- fixed
      saveRDS(weekly_obj, weekly_path)
      repaired <- TRUE
    }
  }

  if (file.exists(population_path)) {
    population_obj <- readRDS(population_path)
    fixed <- repair_dplyr_join_suffix(population_obj$point)
    if (!identical(names(fixed), names(population_obj$point)) ||
        ncol(fixed) != ncol(population_obj$point)) {
      population_obj$point <- fixed
      saveRDS(population_obj, population_path)
      repaired <- TRUE
    }
  }

  if (file.exists(excel_path)) {
    sheets <- readxl::excel_sheets(excel_path)
    repair_sheets <- c("weekly_effects", "population_effects")
    excel_needs_repair <- any(vapply(
      intersect(sheets, repair_sheets),
      function(sheet_name) {
        df <- as.data.frame(readxl::read_excel(excel_path, sheet = sheet_name))
        needs_effects_sheet_repair(df)
      },
      logical(1L)
    ))
    if (excel_needs_repair) {
      sheets_data <- stats::setNames(
        lapply(sheets, function(sheet_name) {
          df <- as.data.frame(readxl::read_excel(excel_path, sheet = sheet_name))
          if (sheet_name %in% repair_sheets) {
            df <- repair_dplyr_join_suffix(df)
          }
          df
        }),
        sheets
      )
      writexl::write_xlsx(sheets_data, path = excel_path)
      repaired <- TRUE
    }
  }

  if (repaired) {
    message("Columnas .x/.y reparadas: ", output_stub)
  }
  invisible(repaired)
}

repair_all_point_estimates <- function(stubs = NULL) {
  if (is.null(stubs)) {
    excel_files <- list.files(
      dir_summary,
      pattern = "_point_estimates\\.xlsx$",
      full.names = FALSE
    )
    stubs <- sub("_point_estimates\\.xlsx$", "", excel_files)
  }
  if (!length(stubs)) {
    message("No hay stubs para reparar columnas .x/.y.")
    return(invisible(NULL))
  }
  repaired_any <- vapply(stubs, repair_point_estimates_stub, logical(1L))
  if (!any(repaired_any)) {
    message("No se encontraron columnas .x/.y ni duplicados en RDS/Excel.")
  }
  invisible(stubs[repaired_any])
}

## ===== Augmentar métricas derivadas sin re-estimación =====
compute_paf_ci_from_boot <- function(boot_path) {
  if (!file.exists(boot_path)) {
    return(NULL)
  }

  boot <- utils::read.csv(boot_path, stringsAsFactors = FALSE)
  if (!all(c("scenario", "prevalence", "iter") %in% names(boot))) {
    return(NULL)
  }

  paf_by_iter <- lapply(split(boot, boot$iter), function(df) {
    nat_prev <- df$prevalence[df$scenario == "observed"][1L]
    int_prev <- df$prevalence[df$scenario == "intervention"][1L]
    if (!is.finite(nat_prev) || nat_prev <= 0) {
      return(NA_real_)
    }
    (nat_prev - int_prev) / nat_prev
  })

  paf_vals <- unlist(paf_by_iter, use.names = FALSE)
  paf_vals <- paf_vals[is.finite(paf_vals)]
  if (!length(paf_vals)) {
    return(NULL)
  }

  c(
    attributable_fraction_lcl = as.numeric(stats::quantile(paf_vals, 0.025, na.rm = TRUE)),
    attributable_fraction_ucl = as.numeric(stats::quantile(paf_vals, 0.975, na.rm = TRUE))
  )
}

augment_population_effects <- function(pop) {
  pop <- dplyr::as_tibble(pop)
  if (!all(c("scenario", "prevalence") %in% names(pop))) {
    stop("population_effects requiere columnas scenario y prevalence.")
  }

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

augment_point_estimates_excel <- function(excel_path, boot_path = NULL) {
  sheets <- readxl::excel_sheets(excel_path)
  pop_sheet <- sheets[grepl("^population_effects$", sheets, ignore.case = TRUE)][1L]
  if (is.na(pop_sheet)) {
    warning("Sin hoja population_effects en ", excel_path)
    return(invisible(FALSE))
  }

  pop <- readxl::read_excel(excel_path, sheet = pop_sheet)
  pop <- augment_population_effects(pop)

  paf_ci <- compute_paf_ci_from_boot(boot_path)
  if (!is.null(paf_ci)) {
    pop$attributable_fraction_lcl[pop$scenario == "intervention"] <- paf_ci[["attributable_fraction_lcl"]]
    pop$attributable_fraction_ucl[pop$scenario == "intervention"] <- paf_ci[["attributable_fraction_ucl"]]
    pop$attributable_fraction_lcl[pop$scenario == "observed"] <- 0
    pop$attributable_fraction_ucl[pop$scenario == "observed"] <- 0
  }

  sheets_data <- stats::setNames(
    lapply(sheets, function(sheet_name) {
      if (sheet_name == pop_sheet) {
        return(pop)
      }
      readxl::read_excel(excel_path, sheet = sheet_name)
    }),
    sheets
  )

  writexl::write_xlsx(sheets_data, path = excel_path)
  invisible(TRUE)
}

augment_all_summary_excels <- function() {
  excel_files <- list.files(
    dir_summary,
    pattern = "_point_estimates\\.xlsx$",
    full.names = TRUE
  )
  if (!length(excel_files)) {
    message("No hay archivos *_point_estimates.xlsx para augmentar.")
    return(invisible(NULL))
  }

  for (excel_path in excel_files) {
    stub <- sub("_point_estimates\\.xlsx$", "", basename(excel_path))
    boot_path <- file.path(dir_bootstrap, stub, "population_boot.csv")
    augment_point_estimates_excel(excel_path, boot_path)
    message("Métricas augmentadas en: ", excel_path)
  }
}

## ===== Lectura de resultados =====
read_population_effects <- function(excel_path) {
  sheets <- readxl::excel_sheets(excel_path)
  hit <- sheets[grepl("^population_effects$", sheets, ignore.case = TRUE)]
  if (!length(hit)) {
    stop("No se encontró hoja 'population_effects' en ", excel_path)
  }
  readxl::read_excel(excel_path, sheet = hit[[1L]])
}

load_scenario_population <- function(output_stub, target_scenario = c("observed", "intervention")) {
  target_scenario <- match.arg(target_scenario)
  excel_path <- file.path(dir_summary, paste0(output_stub, "_point_estimates.xlsx"))
  if (!file.exists(excel_path)) {
    return(NULL)
  }

  pop <- read_population_effects(excel_path) |> augment_population_effects()
  boot_path <- file.path(dir_bootstrap, output_stub, "population_boot.csv")
  paf_ci <- compute_paf_ci_from_boot(boot_path)

  if (!is.null(paf_ci) && target_scenario == "intervention") {
    if (!"attributable_fraction_lcl" %in% names(pop) || is.na(pop$attributable_fraction_lcl[pop$scenario == "intervention"])) {
      pop$attributable_fraction_lcl[pop$scenario == "intervention"] <- paf_ci[["attributable_fraction_lcl"]]
      pop$attributable_fraction_ucl[pop$scenario == "intervention"] <- paf_ci[["attributable_fraction_ucl"]]
    }
  }

  required_cols <- unique(unlist(metric_specs))
  missing_cols <- setdiff(required_cols, names(pop))
  if (length(missing_cols)) {
    stop(
      "Hoja population_effects incompleta en ", excel_path,
      ". Faltan: ", paste(missing_cols, collapse = ", ")
    )
  }

  row <- pop |> dplyr::filter(.data$scenario == target_scenario)
  if (!nrow(row)) {
    stop("No se encontró escenario '", target_scenario, "' en ", excel_path)
  }

  as.list(row[1, , drop = FALSE])
}

format_table_row <- function(scenario_label, pop_row) {
  pct <- table_percent_scale
  c(
    scenario_label,
    format_metric_row(pop_row, "prevalence", pct),
    format_metric_row(pop_row, "cases"),
    format_metric_row(pop_row, "risk_ratio"),
    format_metric_row(pop_row, "risk_difference", pct),
    format_metric_row(pop_row, "attributable_risk", pct),
    format_metric_row(pop_row, "attributable_fraction", pct)
  )
}

scenario_excel_path <- function(stub) {
  file.path(dir_summary, paste0(stub, "_point_estimates.xlsx"))
}

build_pollutant_block <- function(spec) {
  section_label <- spec$section_label
  scenarios <- spec$scenarios

  available <- Filter(
    function(sc) file.exists(scenario_excel_path(sc$stub)),
    scenarios
  )

  if (!length(available)) {
    for (sc in scenarios) {
      message("  [omitido] ", sc$stub, ": no se encontró ", sc$stub, "_point_estimates.xlsx")
    }
    warning("Sin resultados para ", spec$stub_prefix, "; se omite bloque.")
    return(NULL)
  }

  for (sc in scenarios) {
    if (!file.exists(scenario_excel_path(sc$stub))) {
      message("  [omitido] ", sc$stub, ": no se encontró ", sc$stub, "_point_estimates.xlsx")
    }
  }

  natural_row <- load_scenario_population(available[[1L]]$stub, "observed")
  if (is.null(natural_row)) {
    warning("No se pudo leer curso natural para ", spec$stub_prefix)
    return(NULL)
  }

  rows <- list(format_table_row("Natural Course", natural_row))

  for (sc in available) {
    intervention_row <- load_scenario_population(sc$stub, "intervention")
    rows[[length(rows) + 1L]] <- format_table_row(sc$label, intervention_row)
  }

  c(
    list(c(section_label, rep("", length(metric_specs)))),
    rows
  )
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

  openxlsx::writeData(
    wb,
    sheet = "Table",
    x = table_df,
    startRow = 1,
    colNames = TRUE
  )

  openxlsx::addStyle(
    wb,
    sheet = "Table",
    style = header_style,
    rows = 1,
    cols = seq_len(n_cols),
    gridExpand = TRUE
  )

  for (i in seq_len(nrow(table_df))) {
    is_section <- all(table_df[i, 2:n_cols, drop = TRUE] == "")
    openxlsx::addStyle(
      wb,
      sheet = "Table",
      style = if (is_section) section_style else body_style,
      rows = i + 1L,
      cols = 2:n_cols,
      gridExpand = TRUE,
      stack = TRUE
    )
    openxlsx::addStyle(
      wb,
      sheet = "Table",
      style = if (is_section) section_style else label_style,
      rows = i + 1L,
      cols = 1,
      gridExpand = TRUE,
      stack = TRUE
    )
  }

  openxlsx::setColWidths(
    wb,
    sheet = "Table",
    cols = seq_len(n_cols),
    widths = c(28, rep(24, n_cols - 1L))
  )

  dir.create(dirname(path_output), recursive = TRUE, showWarnings = FALSE)
  openxlsx::saveWorkbook(wb, path_output, overwrite = TRUE)
}

## ===== Ejecución =====
message("Reparando columnas .x/.y en RDS/Excel (sin re-estimación)...")
repair_all_point_estimates()

message("Augmentando métricas derivadas en Summary_results (sin re-estimación)...")
augment_all_summary_excels()

message("Construyendo tabla resumen de efectos poblacionales...")
summary_table <- build_summary_table()
write_summary_workbook(summary_table, path_output)

message("Tabla guardada: ", path_output)
print(summary_table)
