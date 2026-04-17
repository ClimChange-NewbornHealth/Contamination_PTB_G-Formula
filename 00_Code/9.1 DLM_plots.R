# 9.1 Plot DLM results -----
rm(list = ls())

## Settings ----
source("00_Code/0.1 Settings.R")
source("00_Code/0.2 Packages.R")
source("00_Code/0.3 Functions.R")

data_inp <- "01_Data/Output/"
data_out <- "02_Output/Descriptives/"
data_out_model <- "02_Output/Models/"

## Load models ----

idw <- rio::import(paste0(data_out_model, "DLM_cox_idw_results.RData"))$results_cox |> 
  as.data.frame() |> 
  select(matches("\\.(term|estimate|conf\\.low|conf\\.high)$")) |> 
  pivot_longer(
    cols = matches("^(pm25|o3|no2)_idw\\.(term|estimate|conf\\.low|conf\\.high)$"),
    names_to = c("pollutant", ".value"),
    names_pattern = "^(pm25|o3|no2)_idw\\.(term|estimate|conf\\.low|conf\\.high)$"
  ) |>
  filter(grepl("^exposicion_[0-9]+$", term)) |> 
  drop_na() |> 
  mutate(
    term = str_extract(term, "[:digit:]+$") |> as.numeric(),
    method = "IDW",
    risk = if_else(conf.low > 1, 1, 0),
    protect = if_else(conf.high < 1, 1, 0)
    ) |> 
  rename(week = term)

glimpse(idw)

krg <- rio::import(paste0(data_out_model, "DLM_cox_krg_results.RData"))$results_cox |> 
  as.data.frame() |> 
  select(matches("\\.(term|estimate|conf\\.low|conf\\.high)$")) |> 
  pivot_longer(
    cols = matches("^(pm25|o3|no2)_krg\\.(term|estimate|conf\\.low|conf\\.high)$"),
    names_to = c("pollutant", ".value"),
    names_pattern = "^(pm25|o3|no2)_krg\\.(term|estimate|conf\\.low|conf\\.high)$"
  ) |>
  filter(grepl("^exposicion_[0-9]+$", term)) |> 
  drop_na() |> 
  mutate(
    term = str_extract(term, "[:digit:]+$") |> as.numeric(),
    method = "Kriging",
    risk = if_else(conf.low > 1, 1, 0),
    protect = if_else(conf.high < 1, 1, 0)
    ) |> 
  rename(week = term)

glimpse(krg)

data_models <- krg |> 
  bind_rows(idw) |> 
  mutate(
    pollutant = factor(
      pollutant, 
      levels = c("pm25", "no2", "o3"), 
      labels = c("PM<sub>2.5</sub>", "NO<sub>2</sub>", "O<sub>3</sub>")),
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
y_step <- y_delta / 3
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
    ggplot(aes(x = week, y = estimate)) +
    geom_rect(
      data = trimester_bands,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill_col),
      inherit.aes = FALSE,
      alpha = 0.15
    ) +
    geom_text(
      data = trimester_labels,
      aes(x = x, y = Inf, label = label),
      inherit.aes = FALSE,
      vjust = 1.2,
      size = 4, 
      fontface = "bold"
    ) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
    geom_errorbar(
      aes(ymin = conf.low, ymax = conf.high, color = hr_color),
      width = 0.3
    ) +
    geom_point(aes(color = hr_color), size = 2) +
    scale_color_manual(
      values = c(
        "Increased risk" = "#E41A1C",
        "Protective" = "#377EB8",
        "Null" = "black"
      ),
      breaks = c("Increased risk", "Protective", "Null"),
      name = NULL
    ) +
    scale_fill_identity() +
    scale_y_continuous(
      limits = y_limits,
      breaks = y_breaks,
      labels = scales::label_number(accuracy = 0.01, decimal.mark = ".")
    ) +
    scale_x_continuous(breaks = seq(1, 37, by = 3)) +
    labs(title = panel_title, y = "HR (95% CI)", x = "Gestational week") +
    theme_light(base_size = 10) +
    theme(
      plot.title = element_markdown(size = 14, hjust = 0),
      legend.position = "none",
      panel.grid = element_blank(),
      axis.text.y = element_text(size = 9),
      axis.text.x = element_text(size = 8),
      plot.margin = margin(4, 4, 4, 4, "pt")
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
  paste0(data_out_model, "DLM_models_krg.png"),
  plot = p_krg,
  res = 300,
  width = 30,
  height = 13,
  scale = 1,
  units = "cm",
  device = ragg::agg_png
)

ggplot2::ggsave(
  paste0(data_out_model, "DLM_models_idw.png"),
  plot = p_idw,
  res = 300,
  width = 30,
  height = 13,
  scale = 1,
  units = "cm",
  device = ragg::agg_png
)
