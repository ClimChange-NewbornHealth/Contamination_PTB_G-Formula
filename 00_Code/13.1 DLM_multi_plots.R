# 13.1 Plot DLM multicontaminante -----
#
# Extensión de 9.1 DLM_plots.R para resultados de 13.0 DLM_multi_pollution.R
# Requisito previo: 13.0 DLM_multi_pollution.R

rm(list = ls())

## Settings ----
source("00_Code/0.1 Settings.R")
source("00_Code/0.2 Packages.R")
source("00_Code/0.3 Functions.R")

data_out_model <- "02_Output/Models/"

krg_path <- paste0(data_out_model, "DLM_multi_cox_krg_results.RData")
idw_path <- paste0(data_out_model, "DLM_multi_cox_idw_results.RData")

if (!file.exists(krg_path)) {
  stop("No se encontró ", krg_path, ". Ejecute primero 00_Code/13.0 DLM_multi_pollution.R")
}
if (!file.exists(idw_path)) {
  stop("No se encontró ", idw_path, ". Ejecute primero 00_Code/13.0 DLM_multi_pollution.R")
}

## Load models ----

idw <- rio::import(idw_path)$results_cox |>
  as.data.frame() |>
  dplyr::select(dplyr::matches("\\.(term|estimate|conf\\.low|conf\\.high)$")) |>
  tidyr::pivot_longer(
    cols = dplyr::matches("^(pm25|o3|no2)_idw\\.(term|estimate|conf\\.low|conf\\.high)$"),
    names_to = c("pollutant", ".value"),
    names_pattern = "^(pm25|o3|no2)_idw\\.(term|estimate|conf\\.low|conf\\.high)$"
  ) |>
  dplyr::filter(grepl("^exposicion_[0-9]+$", term)) |>
  tidyr::drop_na() |>
  dplyr::mutate(
    term = stringr::str_extract(term, "[:digit:]+$") |> as.numeric(),
    method = "IDW",
    risk = dplyr::if_else(conf.low > 1, 1, 0),
    protect = dplyr::if_else(conf.high < 1, 1, 0)
  ) |>
  dplyr::rename(week = term)

glimpse(idw)

krg <- rio::import(krg_path)$results_cox |>
  as.data.frame() |>
  dplyr::select(dplyr::matches("\\.(term|estimate|conf\\.low|conf\\.high)$")) |>
  tidyr::pivot_longer(
    cols = dplyr::matches("^(pm25|o3|no2)_krg\\.(term|estimate|conf\\.low|conf\\.high)$"),
    names_to = c("pollutant", ".value"),
    names_pattern = "^(pm25|o3|no2)_krg\\.(term|estimate|conf\\.low|conf\\.high)$"
  ) |>
  dplyr::filter(grepl("^exposicion_[0-9]+$", term)) |>
  tidyr::drop_na() |>
  dplyr::mutate(
    term = stringr::str_extract(term, "[:digit:]+$") |> as.numeric(),
    method = "Kriging",
    risk = dplyr::if_else(conf.low > 1, 1, 0),
    protect = dplyr::if_else(conf.high < 1, 1, 0)
  ) |>
  dplyr::rename(week = term)

glimpse(krg)

data_models <- krg |>
  dplyr::bind_rows(idw) |>
  dplyr::mutate(
    pollutant = factor(
      pollutant,
      levels = c("pm25", "no2", "o3"),
      labels = c("PM<sub>2.5</sub>", "NO<sub>2</sub>", "O<sub>3</sub>")
    ),
    pollutant_panel = factor(
      dplyr::recode(
        as.character(pollutant),
        "PM<sub>2.5</sub>" = "A. PM<sub>2.5</sub>",
        "NO<sub>2</sub>" = "B. NO<sub>2</sub>",
        "O<sub>3</sub>" = "C. O<sub>3</sub>"
      ),
      levels = c("A. PM<sub>2.5</sub>", "B. NO<sub>2</sub>", "C. O<sub>3</sub>")
    ),
    hr_color = dplyr::case_when(
      risk == 1 ~ "Increased risk",
      protect == 1 ~ "Protective",
      TRUE ~ "Null"
    )
  )

glimpse(data_models)

## Figure with the models ----

y_delta <- max(abs(data_models$conf.low - 1), abs(data_models$conf.high - 1), na.rm = TRUE)
y_step <- y_delta / 2
y_breaks <- 1 + (-3:3) * y_step
y_limits <- range(y_breaks)

trimester_bands <- tibble::tribble(
  ~xmin, ~xmax, ~fill_col,
  -Inf, 12, "gray70",
  12, 24, "white",
  24, Inf, "gray70"
)

trimester_labels <- tibble::tribble(
  ~x, ~label,
  6, "T1",
  18, "T2",
  30.5, "T3"
)

plot_dlm_single <- function(data, method_filter, pollutant_filter, panel_title) {
  data |>
    dplyr::filter(method == method_filter, pollutant == pollutant_filter) |>
    ggplot2::ggplot(ggplot2::aes(x = week, y = estimate)) +
    ggplot2::geom_rect(
      data = trimester_bands,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill_col),
      inherit.aes = FALSE,
      alpha = 0.15
    ) +
    ggplot2::geom_text(
      data = trimester_labels,
      ggplot2::aes(x = x, y = Inf, label = label),
      inherit.aes = FALSE,
      vjust = 1.2,
      size = 4,
      fontface = "bold"
    ) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = conf.low, ymax = conf.high, color = hr_color),
      width = 0.3
    ) +
    ggplot2::geom_point(ggplot2::aes(color = hr_color), size = 2) +
    ggplot2::scale_color_manual(
      values = c(
        "Increased risk" = "#E41A1C",
        "Protective" = "#377EB8",
        "Null" = "black"
      ),
      breaks = c("Increased risk", "Protective", "Null"),
      name = NULL
    ) +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_y_continuous(
      limits = y_limits,
      breaks = y_breaks,
      labels = scales::label_number(accuracy = 0.01, decimal.mark = ".")
    ) +
    ggplot2::scale_x_continuous(breaks = seq(1, 37, by = 3)) +
    ggplot2::labs(
      title = panel_title,
      subtitle = "Multicontaminant-adjusted (weekly co-exposure at w)",
      y = "HR (95% CI)",
      x = "Gestational week"
    ) +
    ggplot2::theme_light(base_size = 10) +
    ggplot2::theme(
      plot.title = ggtext::element_markdown(size = 14, hjust = 0),
      plot.subtitle = ggplot2::element_text(size = 9, hjust = 0, color = "grey30"),
      legend.position = "none",
      panel.grid = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 9),
      axis.text.x = ggplot2::element_text(size = 8),
      plot.margin = ggplot2::margin(4, 4, 4, 4, "pt")
    )
}

plot_dlm_by_method <- function(data, method_filter) {
  p_pm25 <- plot_dlm_single(data, method_filter, "PM<sub>2.5</sub>", "A. PM<sub>2.5</sub>")
  p_no2 <- plot_dlm_single(data, method_filter, "NO<sub>2</sub>", "B. NO<sub>2</sub>")
  p_o3 <- plot_dlm_single(data, method_filter, "O<sub>3</sub>", "C. O<sub>3</sub>")

  ggpubr::ggarrange(
    p_pm25, p_no2, p_o3,
    ncol = 3, nrow = 1,
    align = "hv"
  )
}

p_krg <- plot_dlm_by_method(data_models, "Kriging")
p_idw <- plot_dlm_by_method(data_models, "IDW")

p_krg
p_idw

ggplot2::ggsave(
  paste0(data_out_model, "DLM_multi_models_krg.png"),
  plot = p_krg,
  res = 300,
  width = 30,
  height = 13,
  scale = 1,
  units = "cm",
  device = ragg::agg_png
)

ggplot2::ggsave(
  paste0(data_out_model, "DLM_multi_models_idw.png"),
  plot = p_idw,
  res = 300,
  width = 30,
  height = 13,
  scale = 1,
  units = "cm",
  device = ragg::agg_png
)
