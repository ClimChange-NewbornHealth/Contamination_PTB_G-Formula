# 10.3 G-Formula — figuras a partir de resultados exportados (Excel / Heatmap CSV) ----
#
# Uso (desde la raíz del proyecto, Mac / local):
#   Rscript "00_Code/10.3 G-Form_plots.R"
#
# Salida:
#   02_Output/G-Form/Figures/Figure_cumulative_risk_interventions.png
#   02_Output/G-Form/Figures/Figure_heatmap_rd_interventions.png
#   02_Output/G-Form/Summary_results/heatmap_rd_min_max_by_scenario.xlsx

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

install_load(c("readxl", "dplyr", "tidyr", "ggplot2", "scales", "ggpubr", "grid", "writexl"))

## ===== Configuración =====
data_out_g <- "02_Output/G-Form/"
dir_figures <- file.path(data_out_g, "Figures")
dir_summary <- file.path(data_out_g, "Summary_results")
dir_heatmap <- file.path(data_out_g, "Heatmap")

pollutants <- c("pm25", "no2", "o3")

panel_titles <- list(
  pm25 = expression("A. PM"[2.5]),
  no2 = expression("B. NO"[2]),
  o3 = expression("C. O"[3])
)

cap_labels <- list(
  pm25 = list(lt20 = "< 20 \u00b5g/m\u00b3", lt5 = "< 5 \u00b5g/m\u00b3"),
  no2 = list(lt20 = "< 20 ppbv", lt5 = "< 5 ppbv")
)

pollutant_scenarios <- list(
  pm25 = list(
    list(stub = "pm25_lt20", series = "lt20"),
    list(stub = "pm25_lt5", series = "lt5"),
    list(stub = "pm25_pct20", series = "pct20")
  ),
  no2 = list(
    list(stub = "no2_lt20", series = "lt20"),
    list(stub = "no2_lt5", series = "lt5"),
    list(stub = "no2_pct20", series = "pct20")
  ),
  o3 = list(
    list(stub = "o3_pct20", series = "pct20")
  )
)

heatmap_rows <- list(
  pct20 = list(
    pm25 = "pm25_pct20",
    no2 = "no2_pct20",
    o3 = "o3_pct20"
  ),
  lt20 = list(
    pm25 = "pm25_lt20",
    no2 = "no2_lt20",
    o3 = NA_character_
  ),
  lt5 = list(
    pm25 = "pm25_lt5",
    no2 = "no2_lt5",
    o3 = NA_character_
  )
)


follow_up_weeks <- 28:36
intervention_week_range <- c(0L, 44L)
fig_dpi <- 300
label_5dec <- scales::label_number(accuracy = 0.00001)
label_rd_legend <- scales::label_number(accuracy = 0.00001)
axis_title_size <- 9
fig_width <- 30
fig_height <- 12
heatmap_fig_height <- 34
heatmap_legend_barwidth <- 2.6
heatmap_legend_barheight <- 0.35

heatmap_rd_limits <- c(
  pm25 = 0.00025,
  no2 = 0.00020,
  o3 = 0.00045
)

legend_series_levels <- c(
  "Natural Course",
  "< 20 limit",
  "< 5 limit",
  "Reduced by 20%"
)

series_colours <- c(
  "Natural Course" = "#8C8C8C",
  "< 20 limit" = "#E41A1C",
  "< 5 limit" = "#377EB8",
  "Reduced by 20%" = "#4DAF4A"
)
series_linetypes <- c(
  "Natural Course" = "solid",
  "< 20 limit" = "22",
  "< 5 limit" = "44",
  "Reduced by 20%" = "1313"
)
series_linewidths <- c(
  "Natural Course" = 0.65,
  "< 20 limit" = 0.45,
  "< 5 limit" = 0.45,
  "Reduced by 20%" = 0.45
)

series_key_to_label <- c(
  lt20 = "< 20 limit",
  lt5 = "< 5 limit",
  pct20 = "Reduced by 20%"
)

gform_panel_theme <- function() {
  ggplot2::theme_light(base_size = 12) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = 11, hjust = 0),
      plot.subtitle = ggplot2::element_text(size = 9, hjust = 0, colour = "grey30"),
      axis.title.x = ggplot2::element_text(size = axis_title_size),
      axis.title.y = ggplot2::element_text(size = axis_title_size)
    )
}

gform_excel_sheet <- function(path, pattern) {
  sheets <- readxl::excel_sheets(path)
  hit <- sheets[grepl(pattern, sheets, ignore.case = TRUE)]
  if (!length(hit)) {
    return(NA_character_)
  }
  hit[[1L]]
}

scenario_excel_path <- function(stub) {
  file.path(dir_summary, paste0(stub, "_point_estimates.xlsx"))
}

load_curves_stub <- function(stub) {
  excel_path <- scenario_excel_path(stub)
  if (!file.exists(excel_path)) {
    return(NULL)
  }

  sheet_curves <- gform_excel_sheet(excel_path, "^curvas_riesgo_acumulado")
  if (is.na(sheet_curves)) {
    warning("Sin hoja de curvas en ", excel_path)
    return(NULL)
  }

  curves <- readxl::read_excel(excel_path, sheet = sheet_curves) |>
    dplyr::filter(.data$follow_up_week %in% follow_up_weeks)

  if (!nrow(curves)) {
    warning("Curvas vacías en ", excel_path)
    return(NULL)
  }

  curves
}

load_heatmap_stub <- function(stub) {
  if (is.na(stub) || !nzchar(stub)) {
    return(NULL)
  }

  excel_path <- scenario_excel_path(stub)
  if (file.exists(excel_path)) {
    sheet_heatmap <- gform_excel_sheet(excel_path, "^mapa_calor_rd_semana_interve$")
    if (!is.na(sheet_heatmap)) {
      heatmap_long <- readxl::read_excel(excel_path, sheet = sheet_heatmap)
      if (nrow(heatmap_long) && "risk_difference" %in% names(heatmap_long)) {
        return(
          heatmap_long |>
            dplyr::filter(is.finite(.data$risk_difference))
        )
      }
    }
  }

  csv_path <- file.path(dir_heatmap, stub, "heatmap_long.csv")
  if (file.exists(csv_path)) {
    heatmap_long <- utils::read.csv(csv_path, stringsAsFactors = FALSE)
    if (nrow(heatmap_long) && "risk_difference" %in% names(heatmap_long)) {
      return(
        heatmap_long |>
          dplyr::filter(is.finite(.data$risk_difference))
      )
    }
  }

  warning("No se encontró heatmap para ", stub)
  NULL
}

build_pollutant_curves_long <- function(pollutant) {
  scenarios <- pollutant_scenarios[[pollutant]]
  available <- Filter(
    function(sc) file.exists(scenario_excel_path(sc$stub)),
    scenarios
  )
  if (!length(available)) {
    stop("No hay curvas disponibles para ", pollutant)
  }

  first_curves <- load_curves_stub(available[[1L]]$stub)
  out <- tibble::tibble(
    follow_up_week = first_curves$follow_up_week,
    cumulative_risk = first_curves$cumulative_risk_observed,
    series = "Natural Course"
  )

  for (sc in available) {
    curves <- load_curves_stub(sc$stub)
    if (is.null(curves)) {
      next
    }
    out <- dplyr::bind_rows(
      out,
      tibble::tibble(
        follow_up_week = curves$follow_up_week,
        cumulative_risk = curves$cumulative_risk_intervention,
        series = series_key_to_label[[sc$series]]
      )
    )
  }

  out |>
    dplyr::mutate(
      series = factor(.data$series, levels = legend_series_levels)
    )
}

## ===== Datos =====
dir.create(dir_figures, recursive = TRUE, showWarnings = FALSE)

curves_by_pollutant <- lapply(pollutants, build_pollutant_curves_long)
names(curves_by_pollutant) <- pollutants

## ===== Panel: riesgo acumulado (todas las intervenciones por contaminante) =====
y_vals <- unlist(lapply(curves_by_pollutant, function(df) df$cumulative_risk))
y_max <- max(y_vals, na.rm = TRUE)
y_upper <- ceiling(y_max * 1e5 * 1.02) / 1e5
y_breaks <- scales::pretty_breaks(n = 8)(c(0, y_upper))

plot_cumulative_risk_panel <- function(df, panel_title, y_upper, y_breaks) {
  ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data$follow_up_week,
      y = .data$cumulative_risk,
      colour = .data$series,
      linetype = .data$series,
      linewidth = .data$series
    )
  ) +
    ggplot2::geom_line() +
    ggplot2::scale_x_continuous(
      breaks = follow_up_weeks,
      limits = range(follow_up_weeks)
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, y_upper),
      breaks = y_breaks,
      labels = label_5dec
    ) +
    ggplot2::scale_colour_manual(
      values = series_colours,
      breaks = legend_series_levels,
      drop = FALSE
    ) +
    ggplot2::scale_linetype_manual(
      values = series_linetypes,
      breaks = legend_series_levels,
      drop = FALSE
    ) +
    ggplot2::scale_linewidth_manual(
      values = series_linewidths,
      breaks = legend_series_levels,
      drop = FALSE
    ) +
    ggplot2::labs(
      title = panel_title,
      x = "Gestational Weeks",
      y = "Risk of Birth",
      colour = NULL,
      linetype = NULL,
      linewidth = NULL
    ) +
    gform_panel_theme() +
    ggplot2::theme(
      legend.position = "top",
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 8)
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        order = 1,
        nrow = 1,
        override.aes = list(
          linetype = series_linetypes,
          linewidth = series_linewidths
        )
      ),
      linetype = "none",
      linewidth = "none"
    )
}

p_cumulative_pm25 <- plot_cumulative_risk_panel(
  curves_by_pollutant$pm25,
  panel_titles$pm25,
  y_upper,
  y_breaks
)
p_cumulative_no2 <- plot_cumulative_risk_panel(
  curves_by_pollutant$no2,
  panel_titles$no2,
  y_upper,
  y_breaks
)
p_cumulative_o3 <- plot_cumulative_risk_panel(
  curves_by_pollutant$o3,
  panel_titles$o3,
  y_upper,
  y_breaks
)

p_cumulative_panel <- ggpubr::ggarrange(
  p_cumulative_pm25,
  p_cumulative_no2,
  p_cumulative_o3,
  ncol = 3,
  common.legend = TRUE,
  legend = "top"
)

path_cumulative <- file.path(dir_figures, "Figure_cumulative_risk_interventions.png")

ggplot2::ggsave(
  path_cumulative,
  p_cumulative_panel,
  width = fig_width,
  height = fig_height,
  units = "cm",
  dpi = fig_dpi,
  bg = "white"
)
message("Panel de riesgo acumulado guardado: ", path_cumulative)

## ===== Tabla: min / max RD del heatmap por escenario =====
heatmap_stubs <- unique(unlist(heatmap_rows, use.names = FALSE))
heatmap_stubs <- heatmap_stubs[!is.na(heatmap_stubs) & nzchar(heatmap_stubs)]

heatmap_rd_summary <- dplyr::bind_rows(lapply(heatmap_stubs, function(stub) {
  hm <- load_heatmap_stub(stub)
  if (is.null(hm)) {
    return(NULL)
  }
  hm |>
    dplyr::summarise(
      scenario = stub,
      min_risk_difference = min(.data$risk_difference, na.rm = TRUE),
      max_risk_difference = max(.data$risk_difference, na.rm = TRUE),
      .groups = "drop"
    )
}))

path_heatmap_summary <- file.path(dir_summary, "heatmap_rd_min_max_by_scenario.xlsx")
writexl::write_xlsx(heatmap_rd_summary, path = path_heatmap_summary)
message("Tabla min/max RD guardada: ", path_heatmap_summary)
print(heatmap_rd_summary)

## ===== Panel: mapa de calor RD acumulado (3 filas × 3 columnas) =====
panel_rd_lim <- function(hm, pollutant) {
  cap <- heatmap_rd_limits[[pollutant]]
  if (is.null(hm) || !nrow(hm)) {
    return(cap)
  }
  val <- suppressWarnings(max(abs(hm$risk_difference), na.rm = TRUE))
  if (!is.finite(val) || val <= 0) {
    return(cap)
  }
  val <- min(cap, ceiling(val * 1e5 * 1.02) / 1e5)
  max(val, 1e-5)
}

heatmap_fill_guide <- function(rd_lim) {
  ggplot2::guide_colorbar(
    title.position = "left",
    title.hjust = 1,
    title.vjust = 0.92,
    label.position = "bottom",
    barwidth = grid::unit(heatmap_legend_barwidth, "cm"),
    barheight = grid::unit(heatmap_legend_barheight, "cm"),
    frame.colour = NA,
    ticks.colour = NA,
    title.theme = ggplot2::element_text(
      size = 8,
      hjust = 1,
      margin = ggplot2::margin(0, 3, 0, 0)
    )
  )
}

heatmap_axis_scales <- function(
    intervention_weeks = intervention_week_range,
    follow_up_weeks_hm = follow_up_weeks) {
  list(
    ggplot2::scale_x_continuous(
      limits = c(intervention_week_range[[1L]] - 0.5, intervention_week_range[[2L]] + 0.5),
      breaks = seq(intervention_week_range[[1L]], intervention_week_range[[2L]], by = 5),
      expand = c(0, 0)
    ),
    ggplot2::scale_y_reverse(
      limits = c(max(follow_up_weeks_hm) + 0.5, min(follow_up_weeks_hm) - 0.5),
      breaks = follow_up_weeks_hm,
      expand = c(0, 0)
    )
  )
}

heatmap_fill_scale <- function(rd_lim) {
  ggplot2::scale_fill_gradient2(
    name = "Risk Difference",
    low = "#B2182B",
    mid = "white",
    high = "#542788",
    midpoint = 0,
    limits = c(-rd_lim, rd_lim),
    breaks = c(-rd_lim, rd_lim),
    labels = label_rd_legend(c(-rd_lim, rd_lim)),
    oob = scales::squish,
    guide = heatmap_fill_guide(rd_lim)
  )
}

heatmap_panel_theme <- function(
    show_legend,
    show_x_title = TRUE,
    show_x_text = FALSE,
    show_y = TRUE) {
  gform_panel_theme() +
    ggplot2::theme(
      legend.position = if (show_legend) "top" else "none",
      legend.justification = "center",
      legend.box = "horizontal",
      legend.box.just = "center",
      legend.title = ggplot2::element_text(size = 8, hjust = 1),
      legend.text = ggplot2::element_text(size = 7, vjust = 0),
      legend.background = ggplot2::element_blank(),
      legend.margin = ggplot2::margin(b = 2, t = 0),
      legend.title.align = 0,
      plot.margin = ggplot2::margin(4, 4, 4, 4),
      axis.title.x = if (show_x_title) {
        ggplot2::element_text(size = axis_title_size)
      } else {
        ggplot2::element_blank()
      },
      axis.text.x = if (show_x_text) {
        ggplot2::element_text()
      } else {
        ggplot2::element_blank()
      },
      axis.title.y = if (show_y) {
        ggplot2::element_text(size = axis_title_size)
      } else {
        ggplot2::element_blank()
      },
      axis.text.y = if (show_y) {
        ggplot2::element_text()
      } else {
        ggplot2::element_blank()
      },
      axis.ticks = if (show_y || show_x_text) {
        ggplot2::element_line()
      } else {
        ggplot2::element_blank()
      },
      axis.line = if (show_y || show_x_text) {
        ggplot2::element_line()
      } else {
        ggplot2::element_blank()
      }
    )
}

heatmap_subtitle <- function(pollutant, row_id) {
  if (row_id == "pct20") {
    return("Reduced by 20%")
  }
  cap_labels[[pollutant]][[row_id]]
}

plot_heatmap_rd_panel <- function(
    data,
    main_title,
    subtitle,
    rd_lim,
    show_x_title = TRUE,
    show_x_text = TRUE) {

  if (is.null(data) || !nrow(data)) {
    stop("plot_heatmap_rd_panel requiere datos con filas.")
  }

  follow_up_weeks_hm <- sort(unique(data$follow_up_week))

  ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = .data$intervention_week,
      y = .data$follow_up_week,
      fill = .data$risk_difference
    )
  ) +
    ggplot2::geom_tile(colour = NA) +
    heatmap_axis_scales(follow_up_weeks_hm = follow_up_weeks_hm) +
    heatmap_fill_scale(rd_lim) +
    ggplot2::labs(
      title = main_title,
      subtitle = subtitle,
      x = "Gestational Week of Intervention",
      y = "Follow-up Week"
    ) +
    heatmap_panel_theme(
      show_legend = TRUE,
      show_x_title = show_x_title,
      show_x_text = show_x_text
    )
}

plot_heatmap_no_intervention_panel <- function(
    main_title,
    subtitle,
    show_x_title = TRUE,
    show_x_text = FALSE) {

  x_mid <- mean(intervention_week_range)
  y_mid <- mean(follow_up_weeks)

  ggplot2::ggplot() +
    heatmap_axis_scales() +
    ggplot2::annotate(
      "text",
      x = x_mid,
      y = y_mid,
      label = "No Intervention",
      size = 4.2,
      colour = "grey35",
      fontface = "italic"
    ) +
    ggplot2::labs(
      title = main_title,
      subtitle = subtitle,
      x = "Gestational Week of Intervention",
      y = "Follow-up Week"
    ) +
    heatmap_panel_theme(
      show_legend = FALSE,
      show_x_title = show_x_title,
      show_x_text = show_x_text
    )
}

plot_heatmap_no_intervention_blank_panel <- function(main_title, subtitle) {
  ggplot2::ggplot() +
    ggplot2::annotate(
      "text",
      x = 0.5,
      y = 0.5,
      label = "No Intervention",
      size = 4.2,
      colour = "grey35",
      fontface = "italic"
    ) +
    ggplot2::labs(
      title = main_title,
      subtitle = subtitle
    ) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE, clip = "off") +
    gform_panel_theme() +
    ggplot2::theme(
      legend.position = "none",
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.margin = ggplot2::margin(4, 4, 4, 4)
    )
}

plot_heatmap_missing_panel <- function(
    main_title,
    subtitle,
    show_x_title = TRUE,
    show_x_text = FALSE) {

  x_mid <- mean(intervention_week_range)
  y_mid <- mean(follow_up_weeks)

  ggplot2::ggplot() +
    heatmap_axis_scales() +
    ggplot2::annotate(
      "text",
      x = x_mid,
      y = y_mid,
      label = "Data not available",
      size = 3.8,
      colour = "grey50"
    ) +
    ggplot2::labs(
      title = main_title,
      subtitle = subtitle,
      x = "Gestational Week of Intervention",
      y = "Follow-up Week"
    ) +
    heatmap_panel_theme(
      show_legend = FALSE,
      show_x_title = show_x_title,
      show_x_text = show_x_text
    )
}

build_heatmap_grid <- function() {
  row_ids <- names(heatmap_rows)
  plots <- vector("list", length(row_ids) * length(pollutants))
  idx <- 1L

  for (row_i in seq_along(row_ids)) {
    row_id <- row_ids[[row_i]]
    row_cfg <- heatmap_rows[[row_id]]
    show_x_text <- row_i == length(row_ids)
    show_x_title <- TRUE

    for (col_i in seq_along(pollutants)) {
      pollutant <- pollutants[[col_i]]
      stub <- row_cfg[[pollutant]]
      main_title <- if (row_i == 1L) panel_titles[[pollutant]] else NULL
      subtitle <- heatmap_subtitle(pollutant, row_id)
      is_o3_blank <- pollutant == "o3" && row_id %in% c("lt20", "lt5")

      if (is_o3_blank) {
        p <- plot_heatmap_no_intervention_blank_panel(
          main_title = main_title,
          subtitle = subtitle
        )
      } else if (is.na(stub) || !nzchar(stub)) {
        p <- plot_heatmap_no_intervention_panel(
          main_title = main_title,
          subtitle = subtitle,
          show_x_title = show_x_title,
          show_x_text = show_x_text
        )
      } else {
        hm <- load_heatmap_stub(stub)
        if (is.null(hm)) {
          p <- plot_heatmap_missing_panel(
            main_title = main_title,
            subtitle = subtitle,
            show_x_title = show_x_title,
            show_x_text = show_x_text
          )
        } else {
          rd_lim <- panel_rd_lim(hm, pollutant)
          p <- plot_heatmap_rd_panel(
            data = hm,
            main_title = main_title,
            subtitle = subtitle,
            rd_lim = rd_lim,
            show_x_title = show_x_title,
            show_x_text = TRUE
          )
        }
      }

      plots[[idx]] <- p
      idx <- idx + 1L
    }
  }

  do.call(
    ggpubr::ggarrange,
    c(
      plots,
      list(ncol = length(pollutants), nrow = length(row_ids), common.legend = FALSE)
    )
  )
}

p_heatmap_panel <- build_heatmap_grid()

path_heatmap <- file.path(dir_figures, "Figure_heatmap_rd_interventions.png")

ggplot2::ggsave(
  path_heatmap,
  p_heatmap_panel,
  width = fig_width,
  height = heatmap_fig_height,
  units = "cm",
  dpi = fig_dpi,
  bg = "white"
)
message("Panel de mapa de calor RD guardado: ", path_heatmap)
message(
  "  Filas: pct20, lt20, lt5 | Columnas: PM2.5, NO2, O3 (O3 vacío en lt20/lt5)"
)
message(
  "  Escalas RD por contaminante: PM2.5 [-0.00025, 0.00025], ",
  "NO2 [-0.00020, 0.00020], O3 [-0.00045, 0.00045]"
)

## ===== Eliminar figuras antiguas reemplazadas =====
old_figures <- c(
  file.path(dir_figures, "Figure3_pm25_pct20_cumulative_risk.png"),
  file.path(dir_figures, "Figure4_pm25_pct20_heatmap_rd.png")
)
for (old_path in old_figures) {
  if (file.exists(old_path)) {
    unlink(old_path)
    message("Figura antigua eliminada: ", old_path)
  }
}

message("\nListo. Figuras en: ", dir_figures)
