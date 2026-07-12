# 10.3 G-Formula — figuras a partir de resultados exportados (Excel) ----
#
# Uso (desde la raíz del proyecto, Mac / local):
#   Rscript "00_Code/10.3 G-Form_plots.R"
#
# Entrada:
#   02_Output/G-Form/Summary_results/{pm25,no2,o3}_pct20_point_estimates.xlsx
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

intervention_stubs <- c(
  pm25 = "pm25_pct20",
  no2 = "no2_pct20",
  o3 = "o3_pct20"
)

panel_titles <- list(
  pm25 = expression("A. PM"[2.5]),
  no2 = expression("B. NO"[2]),
  o3 = expression("C. O"[3])
)

follow_up_weeks <- 28:36
fig_dpi <- 300
label_5dec <- scales::label_number(accuracy = 0.00001)
label_rd_legend <- scales::label_number(accuracy = 0.00001)
intervention_label <- "\u2193 20%"
axis_title_size <- 9
fig_width <- 30
fig_height <- 12
annotation_margin_right <- 56

heatmap_rd_limits <- c(
  pm25 = 0.00025,
  no2 = 0.00020,
  o3 = 0.00045
)

series_colours <- c(
  "Observed" = "#8C8C8C",
  "Intervention" = "#377EB8"
)
series_linetypes <- c(
  "Observed" = "solid",
  "Intervention" = "22"
)
series_linewidths <- c(
  "Observed" = 0.65,
  "Intervention" = 0.45
)

add_intervention_annotation <- function(plot, y_pos) {
  plot +
    ggplot2::annotate(
      "text",
      x = Inf,
      y = y_pos,
      label = intervention_label,
      fontface = "bold",
      size = 3.8,
      hjust = -0.28,
      vjust = 0.5
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme(
      plot.margin = ggplot2::margin(5.5, annotation_margin_right, 5.5, 5.5)
    )
}

gform_panel_theme <- function() {
  ggplot2::theme_light(base_size = 12) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = 11, hjust = 0),
      axis.title.x = ggplot2::element_text(size = axis_title_size),
      axis.title.y = ggplot2::element_text(size = axis_title_size)
    )
}

gform_excel_sheet <- function(path, pattern) {
  sheets <- readxl::excel_sheets(path)
  hit <- sheets[grepl(pattern, sheets, ignore.case = TRUE)]
  if (!length(hit)) {
    stop(
      "No se encontró hoja que coincida con '", pattern,
      "' en ", path, ". Hojas: ", paste(sheets, collapse = ", ")
    )
  }
  hit[[1L]]
}

load_gform_excel <- function(stub) {
  excel_path <- file.path(dir_summary, paste0(stub, "_point_estimates.xlsx"))
  if (!file.exists(excel_path)) {
    stop("No se encontró el Excel de resultados: ", excel_path)
  }

  sheet_curves <- gform_excel_sheet(excel_path, "^curvas_riesgo_acumulado")
  sheet_heatmap <- gform_excel_sheet(excel_path, "^mapa_calor_rd_semana_interve$")

  curves <- readxl::read_excel(excel_path, sheet = sheet_curves) |>
    dplyr::filter(.data$follow_up_week %in% follow_up_weeks)

  heatmap_long <- readxl::read_excel(excel_path, sheet = sheet_heatmap)

  if (!nrow(curves)) {
    stop(
      "La hoja de curvas en ", excel_path,
      " no tiene filas para semanas ", min(follow_up_weeks), "-", max(follow_up_weeks)
    )
  }
  if (!nrow(heatmap_long)) {
    stop("La hoja del mapa de calor en ", excel_path, " está vacía.")
  }
  if (!"risk_difference" %in% names(heatmap_long)) {
    stop("La hoja del mapa de calor en ", excel_path, " no contiene risk_difference.")
  }

  heatmap_long <- heatmap_long |>
    dplyr::filter(is.finite(.data$risk_difference))

  list(curves = curves, heatmap = heatmap_long)
}

## ===== Datos =====
gform_data <- lapply(intervention_stubs, load_gform_excel)
names(gform_data) <- names(intervention_stubs)

dir.create(dir_figures, recursive = TRUE, showWarnings = FALSE)

## ===== Panel: riesgo acumulado (Observed + Intervention por contaminante) =====
y_vals <- unlist(lapply(gform_data, function(x) {
  c(x$curves$cumulative_risk_observed, x$curves$cumulative_risk_intervention)
}))
y_max <- max(y_vals, na.rm = TRUE)
y_upper <- ceiling(y_max * 1e5 * 1.02) / 1e5
y_breaks <- scales::pretty_breaks(n = 8)(c(0, y_upper))

plot_cumulative_risk_panel <- function(data, panel_title, y_upper, y_breaks) {
  df <- data |>
    tidyr::pivot_longer(
      cols = c("cumulative_risk_observed", "cumulative_risk_intervention"),
      names_to = "series",
      values_to = "cumulative_risk"
    ) |>
    dplyr::mutate(
      series = factor(
        dplyr::recode(
          .data$series,
          cumulative_risk_observed = "Observed",
          cumulative_risk_intervention = "Intervention"
        ),
        levels = c("Observed", "Intervention")
      )
    )

  p <- ggplot2::ggplot(
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
      breaks = c("Observed", "Intervention")
    ) +
    ggplot2::scale_linetype_manual(
      values = series_linetypes,
      breaks = c("Observed", "Intervention")
    ) +
    ggplot2::scale_linewidth_manual(
      values = series_linewidths,
      breaks = c("Observed", "Intervention")
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
      legend.title = ggplot2::element_blank()
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        order = 1,
        override.aes = list(
          linetype = series_linetypes,
          linewidth = series_linewidths
        )
      ),
      linetype = "none",
      linewidth = "none"
    )

  add_intervention_annotation(
    p,
    y_pos = y_upper * 0.5
  )
}

p_cumulative_pm25 <- plot_cumulative_risk_panel(
  gform_data$pm25$curves,
  panel_titles$pm25,
  y_upper,
  y_breaks
)
p_cumulative_no2 <- plot_cumulative_risk_panel(
  gform_data$no2$curves,
  panel_titles$no2,
  y_upper,
  y_breaks
)
p_cumulative_o3 <- plot_cumulative_risk_panel(
  gform_data$o3$curves,
  panel_titles$o3,
  y_upper,
  y_breaks
)

p_cumulative_panel <- ggpubr::ggarrange(
  p_cumulative_pm25,
  p_cumulative_no2,
  p_cumulative_o3,
  ncol = 3,
  common.legend = TRUE
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
heatmap_rd_summary <- dplyr::bind_rows(lapply(names(gform_data), function(pollutant) {
  gform_data[[pollutant]]$heatmap |>
    dplyr::summarise(
      pollutant = pollutant,
      scenario = intervention_stubs[[pollutant]],
      min_risk_difference = min(.data$risk_difference, na.rm = TRUE),
      max_risk_difference = max(.data$risk_difference, na.rm = TRUE),
      .groups = "drop"
    )
}))

path_heatmap_summary <- file.path(dir_summary, "heatmap_rd_min_max_by_scenario.xlsx")
writexl::write_xlsx(heatmap_rd_summary, path = path_heatmap_summary)
message("Tabla min/max RD guardada: ", path_heatmap_summary)
print(heatmap_rd_summary)

## ===== Panel: mapa de calor RD acumulado =====
plot_heatmap_rd_panel <- function(data, panel_title, rd_lim) {
  intervention_weeks <- sort(unique(data$intervention_week))
  follow_up_weeks_hm <- sort(unique(data$follow_up_week))
  y_mid <- mean(range(follow_up_weeks_hm))

  p <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = .data$intervention_week,
      y = .data$follow_up_week,
      fill = .data$risk_difference
    )
  ) +
    ggplot2::geom_tile(colour = NA) +
    ggplot2::scale_x_continuous(
      limits = c(min(intervention_weeks) - 0.5, max(intervention_weeks) + 0.5),
      breaks = seq(0, max(intervention_weeks), by = 5),
      expand = c(0, 0)
    ) +
    ggplot2::scale_y_reverse(
      limits = c(max(follow_up_weeks_hm) + 0.5, min(follow_up_weeks_hm) - 0.5),
      breaks = follow_up_weeks_hm,
      expand = c(0, 0)
    ) +
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
      guide = ggplot2::guide_colorbar(
        title.position = "left",
        title.hjust = 1,
        title.vjust = 0.92,
        label.position = "bottom",
        barwidth = grid::unit(4.5, "cm"),
        barheight = grid::unit(0.45, "cm"),
        frame.colour = NA,
        ticks.colour = NA,
        title.theme = ggplot2::element_text(
          size = 10,
          hjust = 1,
          margin = ggplot2::margin(0, 4, 0, 0)
        )
      )
    ) +
    ggplot2::labs(
      title = panel_title,
      x = "Gestational Week of Intervention",
      y = "Follow-up Week"
    ) +
    gform_panel_theme() +
    ggplot2::theme(
      legend.position = "top",
      legend.justification = "center",
      legend.box = "horizontal",
      legend.box.just = "center",
      legend.title = ggplot2::element_text(size = 10, hjust = 1),
      legend.text = ggplot2::element_text(size = 9, vjust = 0),
      legend.background = ggplot2::element_blank(),
      legend.margin = ggplot2::margin(b = 6),
      legend.title.align = 0
    )

  add_intervention_annotation(
    p,
    y_pos = y_mid
  )
}

p_heatmap_pm25 <- plot_heatmap_rd_panel(
  gform_data$pm25$heatmap,
  panel_titles$pm25,
  heatmap_rd_limits[["pm25"]]
)
p_heatmap_no2 <- plot_heatmap_rd_panel(
  gform_data$no2$heatmap,
  panel_titles$no2,
  heatmap_rd_limits[["no2"]]
)
p_heatmap_o3 <- plot_heatmap_rd_panel(
  gform_data$o3$heatmap,
  panel_titles$o3,
  heatmap_rd_limits[["o3"]]
)

p_heatmap_panel <- ggpubr::ggarrange(
  p_heatmap_pm25,
  p_heatmap_no2,
  p_heatmap_o3,
  ncol = 3,
  common.legend = FALSE
)

path_heatmap <- file.path(dir_figures, "Figure_heatmap_rd_interventions.png")

ggplot2::ggsave(
  path_heatmap,
  p_heatmap_panel,
  width = fig_width,
  height = fig_height,
  units = "cm",
  dpi = fig_dpi,
  bg = "white"
)
message("Panel de mapa de calor RD guardado: ", path_heatmap)
message("  Escalas RD por panel: PM2.5 [-0.00025, 0.00025], NO2 [-0.00020, 0.00020], O3 [-0.00045, 0.00045]")

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
