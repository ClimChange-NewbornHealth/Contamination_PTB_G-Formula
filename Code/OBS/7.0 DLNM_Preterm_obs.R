# Code 7.0 DLNM PretermBirth ----

## Settings ----
source("Code/0.1 Settings.R")
source("Code/0.2 Functions.R")
source("Code/0.3 Packages.R")

# Data path
data_out <- "Data/Output/"

## Data ---- 

bw <- rio::import(paste0(data_out, "series_births_exposition_pm25_o3_kriging_idw_long_obs", ".RData")) |> 
  drop_na() 

# rs
ids <- bw |>
  slice_sample(n = 50000) |>
  pull(id)

bw <- bw |> 
  filter(id %in% ids) |> 
  arrange(id, week_gest_num) |> 
  group_by(id) |> 
  mutate(tstart = week_gest_num - 1,
         tstop  = week_gest_num,
         event  = if_else((week_gest_num == weeks & birth_preterm == 1), 1, 0)) |> # Event in week 
  ungroup()

glimpse(bw) #  27.509.363 obs

## Specification data ---- 

# Variables 
time <- c("id", "week_gest_num")
vd <- c("event")
vi1 <- c("pm25_krg_week_iqr") 
vi2 <- c("pm25_idw_week_iqr") 
vi3 <- c("o3_krg_week_iqr") 
vi4 <- c("o3_idw_week_iqr") 
vc <- c("sex", "age_group_mom", "educ_group_mom", "job_group_mom", 
        "age_group_dad", "educ_group_dad", "job_group_dad",  "covid", "vulnerability")
trend <- c("year_week1")
        
# Select and transform to wide data 
mat1 <- bw |>
  select(id, week_gest_num, pm25_week) |> 
  #filter(week_gest_num<=37) |>
  #dplyr::select(all_of(c(time, vd, vi1, vc, trend, "weeks"))) |>
  pivot_wider(
    id_cols = id, 
    names_from = "week_gest_num", 
    values_from = "pm25_week",
    names_prefix = "w_") |> 
  arrange(id) |>          # asegúrate de que el orden de filas coincida con bw
  select(starts_with("w_")) 

mat1 <- mat1 |> as.matrix()

mat2 <- bw |>
  select(id, week_gest_num, o3_week) |> 
  #filter(week_gest_num<=37) |>
  #dplyr::select(all_of(c(time, vd, vi1, vc, trend, "weeks"))) |>
  pivot_wider(
    id_cols = id, 
    names_from = "week_gest_num", 
    values_from = "o3_week",
    names_prefix = "w_") |> 
  arrange(id) |>          # asegúrate de que el orden de filas coincida con bw
  select(starts_with("w_")) 

####################################################/
# DLNM -------- 
####################################################/

options(future.globals.maxSize = 3000 * 1024^2)  

lagknots <- logknots(x = 36, nk=1, fun = "ns")

lagknots <- equalknots(x = c(2, 36), # same as the one that we define earlier
                         nk = 1, 
                         fun = "ns")
  
cb1 <- crossbasis(bw$pm25_week, 
                   lag = c(0, 37),
                   argvar = list(fun = "bs", degree = 3, knots =  2), # Cubic
                   arglag = list(fun = "ns", knots = lagknots)) # the form of the curve across all the lags, that is, the lag constrain 

cb2 <- crossbasis(bw$o3_week, 
                   lag = c(0, 37),
                   argvar = list(fun = "bs", degree = 3, knots =  2), # Cubic
                   arglag = list(fun = "ns", knots = lagknots)) # the form of the curve across all the lags, that is, the lag constrain 



form1 <- as.formula(paste("Surv(tstart, tstop, ", "event", ") ~ ", "cb1", 
                              "+ sex + age_group_mom + educ_group_mom + job_group_mom +",
                              "age_group_dad + educ_group_dad + job_group_dad +",
                              "vulnerability + factor(month_week1) + factor(year_week1) + factor(covid)"))

form2 <- as.formula(paste("Surv(tstart, tstop, ", "event", ") ~ ", "cb2", 
                              "+ sex + age_group_mom + educ_group_mom + job_group_mom +",
                              "age_group_dad + educ_group_dad + job_group_dad +",
                              "vulnerability + factor(month_week1) + factor(year_week1) + factor(covid)"))


tic()
mod1 <- coxph(form1, data = bw, ties = "efron", cluster = id)
mod2 <- coxph(form2, data = bw, ties = "efron", cluster = id)
toc()

pred1 <- crosspred(cb1, mod1, cen=0) # pm25_krg_week_iqr
pred2 <- crosspred(cb2, mod2, cen=0) # pm25_idw_week_iqr

# Plots: pm25_week_iqr
plot(pred1, ptype = "3d")
plot(pred1, ptype = "contour")
plot(pred1, ptype = "overall", main="Cummulative effect")

# HR lag response 
plot(pred1, ptype = "slice", var = 1)
plot(pred1, ptype = "slice", var = 2)
plot(pred1, ptype = "slice", var = 3)
plot(pred1, ptype = "slice", var = 4)

# Plots: o3_week_iqr
plot(pred2, ptype = "3d")
plot(pred2, ptype = "contour")
plot(pred2, ptype = "overall", main="Cummulative effect")

# HR lag response 
plot(pred2, ptype = "slice", var = 1)
plot(pred2, ptype = "slice", var = 2)




