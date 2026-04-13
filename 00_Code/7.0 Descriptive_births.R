# Code 6: Descriptive exposition — exposure summary table and plots ----

rm(list = ls())

## Settings ----
source("00_Code/0.1 Settings.R")
source("00_Code/0.2 Packages.R")

data_inp <- "01_Data/Output/"
data_out <- "02_Output/Descriptives/"

## Load exposure (contaminación + metadatos comunales: lat, long, sup) ----
births <- rio::import(paste0(data_inp, "births_2010_2020_exposure.RData"))
births_weeks <- rio::import(paste0(data_inp, "births_2010_2020_exposure_weeks.RData"))
glimpse(births)
glimpse(births_weeks)

## Descriptive table births ----

# Fixed cohort bias 

weeks <- births |> 
  group_by(date_start_week_gest, date_ends_week_gest, weeks) |> 
  summarise(n_gestantes=n(),
            min_semana_gestacion=min(weeks),
            max_semana_gestacion=max(weeks), 
            ultimo_nacimiento=max(date_nac) 
          )

write.xlsx(weeks, paste0(data_out, "Start_ends_gestational_weeks.xlsx"))

# Figure trends preterms

table <- births |> 
  group_by(year_nac) |> 
  summarise(
    tasa_vpt=mean(birth_very_preterm, na.rm=TRUE)*1000,
    tasa_mpt=mean(birth_moderately_preterm, na.rm=TRUE)*1000,
    tasa_lpt=mean(birth_late_preterm, na.rm=TRUE)*1000,
    tasa_pt=mean(birth_preterm, na.rm=TRUE)*1000,
    tasa_t=mean(birth_term, na.rm=TRUE)*1000,
    tasa_post=mean(birth_posterm, na.rm=TRUE)*1000,
  ) |> 
  pivot_longer(
    cols=!year_nac, 
    names_to="preterm",
    values_to="prev"
  ) |> 
  mutate(preterm=case_when(
    preterm=="tasa_vpt" ~ "Very Preterm (28-31 weeks)",
    preterm=="tasa_mpt" ~ "Moderately Preterm (32-33 weeks)",
    preterm=="tasa_lpt" ~ "Late Preterm (34-36 weeks)",
    preterm=="tasa_pt" ~ "Preterm (<37 weeks)",
    preterm=="tasa_t" ~  "Term (38-42 weeks)",
    preterm=="tasa_post" ~ "Post-term (>42 weeks)"
  )) |> 
  mutate(preterm=factor(preterm, levels=c(
    "Very Preterm (28-31 weeks)",
    "Moderately Preterm (32-33 weeks)",
    "Late Preterm (34-36 weeks)",
    "Preterm (<37 weeks)",
    "Term (38-42 weeks)",
    "Post-term (>42 weeks)"
  )))

table |> 
  ggplot(aes(y=prev, x=year_nac)) +
  geom_line(color="#08519c") +
  geom_point(color="#08519c") +
  facet_wrap(~preterm, ncol = 2, scales = "free") +
  scale_x_continuous(breaks = seq(2010, 2020, by=1)) +
  labs(y ="Prevalence (per 1.000)", x=NULL) +
  theme_light() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid       = element_blank(),
    strip.background = element_rect(fill = "white"),
    strip.text       = element_text(size = 11, color = "black", hjust=0),
    axis.text.y      = element_text(size = 9),
    axis.ticks.y     = element_blank()
  )

ggsave(filename = paste0(data_out, "Preterm_trends_2010_2020", ".png"), # "Preterm_trendsrm1991"
       res = 300,
       width = 20,
       height = 12,
       units = 'cm',
       scaling = 0.90,
       device = ragg::agg_png
      )

# Descriptive stats table 

tab1 <-  births |> 
   select(
          birth_preterm,
          tbw, weeks, sex,  
          age_group_mom, educ_group_mom, job_group_mom,
          age_group_dad, educ_group_dad, job_group_dad, 
          year_nac, month_nac, 
          sovi, vulnerability
          
   ) |> 
   mutate(
    birth_preterm=factor(birth_preterm)
   ) |> 
   st(,
   digits = 1, 
   out="return", 
   add.median = TRUE,
   fixed.digits = TRUE, 
   simple.kable = FALSE,
   title="",
   numformat = NA) |> 
   data.frame() 

tab1

tab2 <-  births |> 
   filter(birth_preterm==1) |> 
   select(
          tbw, weeks, sex,  
          age_group_mom, educ_group_mom, job_group_mom,
          age_group_dad, educ_group_dad, job_group_dad, 
          year_nac, month_nac, 
          sovi, vulnerability
          
   ) |>  
   st(,
   digits = 1, 
   out="return", 
   add.median = TRUE,
   fixed.digits = TRUE, 
   simple.kable = FALSE,
   title="",
   numformat = NA) |> 
   data.frame() 

tab2

writexl::write_xlsx(list(tab1, tab2), path =  paste0(data_out,  "Descriptives_stats", ".xlsx"))

## Desccriptive exposure -----

glimpse(births_weeks)

births_weeks_short <- births_weeks |> 
  dplyr::select(
    id, com, weeks, date_start_week_gest, date_ends_week_gest, week_gest_num,
    pm25_krg:no2_idw 
  )

table_weeks <- births_weeks_short |> 
  group_by(week_gest_num) |> 
  summarise(
    pm25_krg_mean=mean(pm25_krg, na.rm=TRUE),
    pm25_krg_min=min(pm25_krg, na.rm=TRUE),
    pm25_krg_max=max(pm25_krg, na.rm=TRUE),
    no2_krg_mean=mean(no2_krg, na.rm=TRUE),
    no2_krg_min=min(no2_krg, na.rm=TRUE),
    no2_krg_max=max(no2_krg, na.rm=TRUE),
    o3_krg_mean=mean(o3_krg, na.rm=TRUE),
    o3_krg_min=min(o3_krg, na.rm=TRUE),
    o3_krg_max=max(o3_krg, na.rm=TRUE),

    pm25_idw_mean=mean(pm25_idw, na.rm=TRUE),
    pm25_idw_min=min(pm25_idw, na.rm=TRUE),
    pm25_idw_max=max(pm25_idw, na.rm=TRUE),
    no2_idw_mean=mean(no2_idw, na.rm=TRUE),
    no2_idw_min=min(no2_idw, na.rm=TRUE),
    no2_idw_max=max(no2_idw, na.rm=TRUE),
    o3_idw_mean=mean(o3_idw, na.rm=TRUE),
    o3_idw_min=min(o3_idw, na.rm=TRUE),
    o3_idw_max=max(o3_idw, na.rm=TRUE)
  ) |> 
  ungroup()

vars <- colnames(table_weeks)[-1]

table_weeks <- table_weeks |>
  mutate(across(all_of(vars), ~ round(.x, 2))) |> 
  mutate(across(all_of(vars), ~ formatC(.x, format = "f", digits = 1, decimal.mark = ".")))

glimpse(table_weeks)

writexl::write_xlsx(table_weeks, path =  paste0(data_out,  "Descriptives_exposure_stats_time", ".xlsx"))
