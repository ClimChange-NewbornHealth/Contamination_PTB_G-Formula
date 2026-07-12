# 10.4 G-Formula — tabla resumen de efectos poblacionales ----
#
# Uso (desde la raíz del proyecto):
#   Rscript "00_Code/10.4 G-Form_table.R"
#
# Entrada:
#   02_Output/G-Form/Summary_results/{pm25,no2,o3}_pct{p}_point_estimates.xlsx
#   02_Output/G-Form/Bootstrap/{stub}/population_boot.csv  (IC de PAF si existe)
# Salida:
#   02_Output/G-Form/Summary_results/Table_population_effects_summary.xlsx
#
# Nota: augmenta métricas derivadas (PAF) en Excel existentes sin re-estimación.
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
path_output <- file.path(dir_summary, "Table_population_effects_summary.xlsx")

reduction_pcts <- c(5, 10, 15, 20)

pollutant_specs <- list(
  pm25 = list(
    stub_prefix = "pm25",
    section_label = "PM2.5 (\u00b5g/m\u00b3)"
  ),
  no2 = list(
    stub_prefix = "no2",
    section_label = "NO2 (ppbv)"
  ),
  o3 = list(
    stub_prefix = "o3",
    section_label = "O3 (ppbv)"
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

build_pollutant_block <- function(spec, reduction_pcts) {
  stub_prefix <- spec$stub_prefix
  section_label <- spec$section_label

  available_stubs <- Filter(
    function(pct) {
      file.exists(file.path(dir_summary, paste0(stub_prefix, "_pct", pct, "_point_estimates.xlsx")))
    },
    reduction_pcts
  )

  if (!length(available_stubs)) {
    missing <- paste0(stub_prefix, "_pct", reduction_pcts)
    for (stub in missing) {
      message("  [omitido] ", stub, ": no se encontró ", stub, "_point_estimates.xlsx")
    }
    warning("Sin resultados para ", stub_prefix, "; se omite bloque.")
    return(NULL)
  }

  for (pct in setdiff(reduction_pcts, available_stubs)) {
    stub <- paste0(stub_prefix, "_pct", pct)
    message("  [omitido] ", stub, ": no se encontró ", stub, "_point_estimates.xlsx")
  }

  first_stub <- paste0(stub_prefix, "_pct", available_stubs[[1L]])
  natural_row <- load_scenario_population(first_stub, "observed")
  if (is.null(natural_row)) {
    warning("No se pudo leer curso natural para ", stub_prefix)
    return(NULL)
  }

  rows <- list(format_table_row("Natural Course", natural_row))

  for (pct in available_stubs) {
    stub <- paste0(stub_prefix, "_pct", pct)
    intervention_row <- load_scenario_population(stub, "intervention")
    rows[[length(rows) + 1L]] <- format_table_row(
      paste0("Reduced by ", pct, "%"),
      intervention_row
    )
  }

  c(
    list(c(section_label, rep("", length(metric_specs)))),
    rows
  )
}

build_summary_table <- function() {
  blocks <- lapply(pollutant_specs, build_pollutant_block, reduction_pcts = reduction_pcts)
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
message("Augmentando métricas derivadas en Summary_results (sin re-estimación)...")
augment_all_summary_excels()

message("Construyendo tabla resumen de efectos poblacionales...")
summary_table <- build_summary_table()
write_summary_workbook(summary_table, path_output)

message("Tabla guardada: ", path_output)
print(summary_table)
