# 11.2 Exposure–perinatal outcome models — figuras Cox ----
#
# Uso:
#   Rscript "00_Code/11.2 Exposure_PO_Plots_Models.R"
#
# Requiere: 02_Output/Exposure_PO/Models/Exposure_models_PO_cox.xlsx
# Salida:   02_Output/Exposure_PO/Figures/
#           Figure_HR_preterm_pollutants_raw.png
#           Figure_HR_preterm_pollutants_iqr.png

source("00_Code/0.1 Settings.R")

install_load <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
}

install_load(c("readxl", "dplyr", "ggplot2", "scales", "ggpubr"))
source("00_Code/11.0 Exposure_PO_functions.R")

PO_EXPOSURE_PLOT_LEVELS <- c(
  "Trimester 1", "Trimester 2", "Trimester 3", "Overall"
)

pick_pretty_breaks <- function(limits, n_min = 5L, n_max = 6L, center = 1) {
  breaks <- pretty(limits, n = n_min - 1L)
  if (length(breaks) < n_min) {
    breaks <- pretty(limits, n = n_max - 1L)
  }
  if (length(breaks) > n_max) {
    idx <- unique(round(seq(1, length(breaks), length.out = n_max)))
    breaks <- breaks[idx]
  }
  if (center >= min(breaks) && center <= max(breaks) &&
      !any(abs(breaks - center) < 1e-9)) {
    breaks <- sort(unique(c(breaks, center)))
  }
  if (length(breaks) > n_max) {
    center_idx <- which.min(abs(breaks - center))
    other_idx <- setdiff(seq_along(breaks), center_idx)
    keep_other <- unique(round(seq(1, length(other_idx), length.out = n_max - 1L)))
    breaks <- sort(c(breaks[center_idx], breaks[other_idx[keep_other]]))
  }
  breaks
}

auto_panel_y_axis <- function(
    data_panel,
    n_breaks_min = 5L,
    n_breaks_max = 6L,
    pad_frac = 0.10,
    center = 1,
    limit_margin_frac = 0.15) {
  y_vals <- c(data_panel$estimate, data_panel$conf.low, data_panel$conf.high)
  y_vals <- y_vals[is.finite(y_vals)]
  if (!length(y_vals)) return(NULL)

  max_dist <- max(center - min(y_vals), max(y_vals) - center)
  pad <- max(max_dist * pad_frac, 0.02)
  half_span <- max_dist + pad

  limits <- c(center - half_span, center + half_span)
  breaks <- pick_pretty_breaks(limits, n_breaks_min, n_breaks_max, center)
  half_span <- max(center - min(breaks), max(breaks) - center)

  if (length(breaks) > 1L) {
    step <- min(diff(breaks))
    half_span <- half_span + step * limit_margin_frac
  }

  limits <- c(center - half_span, center + half_span)
  breaks <- breaks[breaks >= limits[[1L]] - 1e-9 & breaks <= limits[[2L]] + 1e-9]
  labels <- sprintf("%.2f", breaks)

  list(limits = limits, breaks = breaks, labels = labels)
}

data_out <- "02_Output/Exposure_PO/"
dir_models <- file.path(data_out, "Models")
dir_figures <- file.path(data_out, "Figures")
dir.create(dir_figures, recursive = TRUE, showWarnings = FALSE)

path_results <- file.path(dir_models, "Exposure_models_PO_cox.xlsx")
if (!file.exists(path_results)) {
  stop("Ejecute primero 11.1 Exposure_PO_Models.R")
}

models_cox <- readxl::read_excel(path_results, sheet = "cox_models")

contaminant_labels <- list(
  pm25 = expression("PM"[2.5]),
  no2 = expression("NO"[2]),
  o3 = expression("O"[3])
)

plot_preterm_panel <- function(
    data_preterm,
    panel_title,
    y_axis,
    show_legend = FALSE) {

  if (!nrow(data_preterm) || is.null(y_axis)) return(NULL)

  pd <- ggplot2::position_dodge(width = 0.6)
  rect_data <- data.frame(xmin = 3.5, xmax = 4.5, ymin = -Inf, ymax = Inf)

  ggplot2::ggplot(
    data_preterm,
    ggplot2::aes(
      y = .data$estimate,
      x = .data$exposure,
      colour = .data$adjustment,
      shape = .data$adjustment
    )
  ) +
    ggplot2::geom_rect(
      data = rect_data,
      ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax, ymin = .data$ymin, ymax = .data$ymax),
      inherit.aes = FALSE,
      fill = "grey95",
      alpha = 0.7
    ) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$conf.low, ymax = .data$conf.high),
      width = 0.25,
      position = pd
    ) +
    ggplot2::geom_point(size = 2, position = pd) +
    ggplot2::scale_color_manual(values = c("Unadjusted" = "grey50", "Adjusted" = "black")) +
    ggplot2::scale_shape_manual(values = c("Unadjusted" = 16, "Adjusted" = 15)) +
    ggplot2::scale_y_continuous(
      limits = y_axis$limits,
      breaks = y_axis$breaks,
      labels = y_axis$labels,
      expand = c(0, 0)
    ) +
    ggplot2::scale_x_discrete(expand = c(0.05, 0)) +
    ggplot2::labs(y = "HR (95% CI)", x = NULL, title = panel_title) +
    ggplot2::theme_light(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 10, face = "bold", hjust = 0.5),
      legend.position = if (show_legend) "top" else "none",
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 9),
      panel.grid = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 9),
      axis.text.x = ggplot2::element_text(size = 8),
      plot.margin = ggplot2::margin(4, 4, 4, 4, "pt")
    )
}

build_preterm_figure <- function(plot_data_cox, scale_key, out_fname) {
  plot_data_cox <- plot_data_cox |>
    dplyr::filter(.data$dependent_var == "birth_preterm") |>
    dplyr::mutate(
      exposure = factor(.data$exposure, levels = PO_EXPOSURE_PLOT_LEVELS),
      adjustment = factor(.data$adjustment, levels = c("Unadjusted", "Adjusted"))
    )

  plots <- lapply(seq_along(PO_EXPOSURE_POLLUTANTS), function(i) {
    cont <- PO_EXPOSURE_POLLUTANTS[[i]]
    data_panel <- plot_data_cox |> dplyr::filter(.data$contaminante == cont)
    y_axis <- if (identical(scale_key, "iqr")) {
      auto_panel_y_axis(data_panel, n_breaks_min = 5L, n_breaks_max = 6L)
    } else {
      auto_panel_y_axis(data_panel, n_breaks_min = 5L, n_breaks_max = 7L)
    }
    plot_preterm_panel(
      data_panel,
      panel_title = contaminant_labels[[cont]],
      y_axis = y_axis,
      show_legend = (i == 1L)
    )
  })
  plots <- plots[!vapply(plots, is.null, logical(1L))]
  if (!length(plots)) return(invisible(NULL))

  fig <- ggpubr::ggarrange(
    plotlist = plots,
    ncol = 3,
    nrow = 1,
    common.legend = TRUE,
    legend = "top",
    align = "none"
  )

  out <- file.path(dir_figures, out_fname)
  ggplot2::ggsave(
    out, fig, width = 22, height = 10, units = "cm", dpi = 300, bg = "white"
  )
  message("Guardado: ", out)
  invisible(fig)
}

old_figs <- list.files(dir_figures, pattern = "\\.png$", full.names = TRUE)
if (length(old_figs)) {
  unlink(old_figs)
  message("Eliminadas ", length(old_figs), " figuras anteriores.")
}

plot_data_all <- prepare_po_plot_table_data(models_cox)
compiled_figures <- list()

for (scale in c("raw", "iqr")) {
  plot_data_scale <- plot_data_all |> dplyr::filter(.data$exposure_scale == scale)
  fname <- paste0("Figure_HR_preterm_pollutants_", scale, ".png")
  fig <- build_preterm_figure(plot_data_scale, scale, fname)
  if (!is.null(fig)) compiled_figures[[paste0("preterm_", scale)]] <- fig
}

save(compiled_figures, file = file.path(dir_models, "Plots_compiled_figures.RData"))
message("Figuras en: ", dir_figures)
