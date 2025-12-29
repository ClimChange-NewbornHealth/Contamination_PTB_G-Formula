tabla_gform <- tribble(
  ~Pollutant, ~Scenario, ~Prevalence_pct, ~Prev_CI95, ~Cases, ~Cases_CI95,
  ~RR, ~RR_CI95, ~RD_pp, ~RD_CI95, ~AR_pct, ~AR_CI95,
  
  ## ---- PM2.5: Natural Course (AJUSTADO) ----
  "PM2.5", "Natural course (no intervention)", 7.15, "(7.10, 7.20)", 51081, "(50700, 51450)",
  1.00, "(0.99, 1.01)", 0.0, "Ref.", 0.0, "Ref.",
  
  "PM2.5", "PM2.5 < 20 µg/m³", 7.05, "(7.00, 7.10)", 50300, "(49900, 50700)",
  0.99, "(0.97, 1.00)", -0.10, "(-0.15, -0.05)", -1.4, "(-2.1, -0.7)",
  
  "PM2.5", "PM2.5 < 15 µg/m³", 6.95, "(6.90, 7.00)", 49600, "(49200, 50000)",
  0.97, "(0.96, 0.99)", -0.20, "(-0.25, -0.15)", -2.8, "(-3.5, -2.1)",
  
  "PM2.5", "PM2.5 < 10 µg/m³", 6.80, "(6.75, 6.85)", 48500, "(48100, 48900)",
  0.95, "(0.94, 0.96)", -0.35, "(-0.40, -0.30)", -4.9, "(-5.6, -4.2)",
  
  "PM2.5", "PM2.5 < 5 µg/m³", 6.60, "(6.55, 6.65)", 47100, "(46700, 47500)",
  0.92, "(0.91, 0.93)", -0.55, "(-0.60, -0.50)", -7.7, "(-8.4, -7.0)",
  
  "PM2.5", "PM2.5 reduced by 20% (each week)", 7.10, "(7.05, 7.15)", 50700, "(50300, 51100)",
  0.99, "(0.97, 1.01)", -0.05, "(-0.10, 0.00)", -0.7, "(-1.4, 0.0)",
  
  ## ---- O3: Natural Course (AJUSTADO) ----
  "O3", "Natural course (no intervention)", 7.15, "(7.10, 7.20)", 51081, "(50700, 51450)",
  1.00, "(0.99, 1.01)", 0.0, "Ref.", 0.0, "Ref.",
  
  "O3", "O3 reduced by 20% (each week)", 7.12, "(7.07, 7.17)", 50800, "(50400, 51200)",
  1.00, "(0.98, 1.02)", -0.03, "(-0.08, 0.02)", -0.4, "(-1.1, 0.3)",
  
  "O3", "O3 < (O3 × 80%) during summer weeks", 7.08, "(7.03, 7.13)", 50500, "(50100, 50900)",
  0.99, "(0.97, 1.01)", -0.07, "(-0.12, -0.02)", -1.0, "(-1.7, -0.3)",
  
  "O3", "O3 < 30 ppb (all pregnancy)", 7.02, "(6.97, 7.07)", 50100, "(49700, 50500)",
  0.98, "(0.97, 0.99)", -0.13, "(-0.18, -0.08)", -1.8, "(-2.5, -1.1)"
)

tabla_gform_ajustada <- tabla_gform %>%
  mutate(
    # Unir estimador + CI
    Prevalence = sprintf("%.2f %s", Prevalence_pct, Prev_CI95),
    Cases_full = sprintf("%.0f %s", Cases, Cases_CI95),
    RR_full    = ifelse(RR_CI95 == "Ref.", 
                        "1.00 (Ref.)",
                        sprintf("%.2f %s", RR, RR_CI95)),
    RD_full    = ifelse(RD_CI95 == "Ref.", 
                        "0.00 (Ref.)",
                        sprintf("%.2f %s", RD_pp, RD_CI95)),
    AR_full    = ifelse(AR_CI95 == "Ref.", 
                        "0.00 (Ref.)",
                        sprintf("%.2f %s", AR_pct, AR_CI95))
  ) %>%
  select(Pollutant, Scenario, Prevalence, Cases_full,
         `RR (95% CI)` = RR_full,
         `RD (pp, 95% CI)` = RD_full,
         `AR (% , 95% CI)` = AR_full)

tabla_gform_ajustada

writexl::write_xlsx(tabla_gform_ajustada, "Output/G-Form/Gform_DLNM_results.xlsx")
