# 10.3 G-Formula — figuras a partir de resultados exportados (Excel) ----
#
# Uso (desde la raíz del proyecto, Mac / local):
#   Rscript "00_Code/10.3 G-Form_plots.R"
#
# Entrada por defecto:
#   02_Output/G-Form/Summary_results/pm25_pct20_point_estimates.xlsx
# Salida:
#   02_Output/G-Form/Figures/

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

install_load(c("readxl", "dplyr", "ggplot2", "scales", "patchwork", "grid"))

## ===== Configuración =====
output_stub <- Sys.getenv("GFORM_PLOT_STUB", unset = "pm25_pct20")
data_out_g <- "02_Output/G-Form/"
dir_figures <- file.path(data_out_g, "Figures")
excel_path <- file.path(
  data_out_g,
  "Summary_results",
  paste0(output_stub, "_point_estimates.xlsx")
)

follow_up_weeks <- 28:36
fig_dpi <- 300
label_5dec <- scales::label_number(accuracy = 0.00001)

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

## ===== Datos =====
if (!file.exists(excel_path)) {
  stop("No se encontró el Excel de resultados: ", excel_path)
}

sheet_curves <- gform_excel_sheet(excel_path, "^curvas_riesgo_acumulado")
sheet_heatmap <- gform_excel_sheet(excel_path, "^mapa_calor_rd_semana_interve$")

curves <- readxl::read_excel(excel_path, sheet = sheet_curves) |>
  dplyr::filter(.data$follow_up_week %in% follow_up_weeks)

heatmap_long <- readxl::read_excel(excel_path, sheet = sheet_heatmap)

if (!nrow(curves)) {
  stop("La hoja de curvas no tiene filas para semanas ", min(follow_up_weeks), "-", max(follow_up_weeks))
}
if (!nrow(heatmap_long)) {
  stop("La hoja del mapa de calor está vacía.")
}
if (!"risk_difference" %in% names(heatmap_long)) {
  stop("La hoja del mapa de calor no contiene la columna risk_difference.")
}

heatmap_long <- heatmap_long |>
  dplyr::filter(is.finite(.data$risk_difference))

intervention_weeks <- sort(unique(heatmap_long$intervention_week))
follow_up_weeks_hm <- sort(unique(heatmap_long$follow_up_week))

dir.create(dir_figures, recursive = TRUE, showWarnings = FALSE)

## ===== Figura 3: riesgo acumulado (Observed e Intervention por separado) =====
y_vals <- c(curves$cumulative_risk_observed, curves$cumulative_risk_intervention)
y_max <- max(y_vals, na.rm = TRUE)
y_upper <- ceiling(y_max * 1e5 * 1.02) / 1e5
y_breaks <- scales::pretty_breaks(n = 8)(c(0, y_upper))

plot_cumulative_risk <- function(data, y_var, series_label, linetype_val) {
  df <- data |>
    dplyr::transmute(
      follow_up_week = .data$follow_up_week,
      cumulative_risk = .data[[y_var]],
      series = series_label
    )

  ggplot2::ggplot(df, ggplot2::aes(x = follow_up_week, y = cumulative_risk, linetype = series)) +
    ggplot2::geom_line(linewidth = 0.7, colour = "black") +
    ggplot2::scale_x_continuous(
      breaks = follow_up_weeks,
      limits = range(follow_up_weeks)
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, y_upper),
      breaks = y_breaks,
      labels = label_5dec
    ) +
    ggplot2::scale_linetype_manual(values = stats::setNames(linetype_val, series_label)) +
    ggplot2::labs(
      x = "Gestational Weeks",
      y = "Risk of Birth",
      linetype = NULL
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      legend.position = "top",
      legend.justification = "left",
      legend.background = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )
}

p_observed <- plot_cumulative_risk(
  curves,
  "cumulative_risk_observed",
  "Observed",
  "solid"
)
p_intervention <- plot_cumulative_risk(
  curves,
  "cumulative_risk_intervention",
  "Intervention",
  "22"
)

p_fig3_panel <- p_observed + p_intervention +
  patchwork::plot_layout(ncol = 2)

path_fig3_panel <- file.path(
  dir_figures,
  paste0("Figure3_", output_stub, "_cumulative_risk.png")
)

ggplot2::ggsave(
  path_fig3_panel,
  p_fig3_panel,
  width = 14,
  height = 5,
  dpi = fig_dpi
)
message("Figura 3 (panel Observed + Intervention) guardada: ", path_fig3_panel)

## ===== Figura 4: mapa de calor RD acumulado =====
rd_lim <- 0.00021
label_rd_legend <- scales::label_number(accuracy = 0.00001)

fig4_width <- 10
fig4_height <- 6

p_fig4 <- ggplot2::ggplot(
  heatmap_long,
  ggplot2::aes(
    x = intervention_week,
    y = follow_up_week,
    fill = risk_difference
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
    x = "Gestational Week of Intervention",
    y = "Follow-up Week"
  ) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
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

path_fig4_png <- file.path(
  dir_figures,
  paste0("Figure4_", output_stub, "_heatmap_rd.png")
)

ggplot2::ggsave(
  path_fig4_png,
  p_fig4,
  width = fig4_width,
  height = fig4_height,
  dpi = fig_dpi
)
message("Figura 4 guardada: ", path_fig4_png)
message("  Escala RD (leyenda): [-0.00021, 0.00021]")

message("\nListo. Figuras en: ", dir_figures)
