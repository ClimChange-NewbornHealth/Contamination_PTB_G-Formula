# Code 8.0 DLNM PretermBirth ----

## Settings ----
source("Code/0.1 Settings.R")
source("Code/0.2 Packages.R")
source("Code/0.3 Functions.R")

# Data path
data_out <- "Data/Output/"

## Data ---- 

bw <- rio::import(paste0(data_out, "series_births_exposition_pm25_o3_kriging_idw_long", ".RData")) |> 
  drop_na() 

bw <- bw |> 
  arrange(id, week_gest_num) |> 
  group_by(id) |> 
  mutate(tstart = week_gest_num - 1,
         tstop  = week_gest_num,
         event  = if_else((week_gest_num == weeks & birth_preterm == 1), 1, 0)) |> # Event in week 
  ungroup()

glimpse(bw) #  27.509.363 obs

## Specification data ---- 

# rs
set.seed(1234)
ids <- bw |>
  slice_sample(n=50000) |>
  pull(id)

bw_sample <- bw |> filter(id %in% ids)

bw_sample <- bw_sample |> 
  select(id:name_com, week_gest_num,
         pm25_krg_week_iqr,
         o3_krg_week_iqr,
         tstart, tstop, event,
         weeks, sex,  
         age_group_mom, educ_group_mom, job_group_mom,
         age_group_dad, educ_group_dad,job_group_dad,
         year_week1, month_week1, vulnerability) 

glimpse(bw_sample)
setDT(bw_sample)
# Remove original data 
#rm(bw)

####################################################/
# DLNM -------- 
####################################################/

options(future.globals.maxSize = 3000 * 1024^2)  

range_pm25 <- range(bw_sample$pm25_krg_week_iqr, na.rm = TRUE); range_pm25
range_o3   <- range(bw_sample$o3_krg_week_iqr, na.rm = TRUE); range_o3
q_pm25 <- equalknots(bw_sample$pm25_krg_week_iqr, nk = 2); q_pm25
q_o3 <- equalknots(bw_sample$o3_krg_week_iqr, nk = 2); q_o3 

# Equal efect
lagknots <- equalknots(x = c(0, 36), nk = 6, fun = "ns"); lagknots

cb1 <- crossbasis(bw_sample$pm25_krg_week_iqr, 
                   lag = c(0, 36),
                   #argvar = list(fun = "lin"),
                   argvar = list(fun = "bs", degree = 2,
                                 knots = q_pm25, #knots = q_pm25,
                                Boundary.knots = range_pm25
                                ), # Cubic Boundary.knots = range_pm25
                   arglag = list(fun = "ns", knots = lagknots)) #knots = lagknots the form of the curve across all the lags, that is, the lag constrain 

form1 <- as.formula(paste("Surv(tstart, tstop, ", "event", ") ~ ", "cb1", 
                              "+ sex + age_group_mom + educ_group_mom + job_group_mom +",
                              "age_group_dad + educ_group_dad + job_group_dad +",
                              "factor(month_week1) + factor(year_week1) + vulnerability")) # strata(year_week1, month_week1)

tic()
mod1 <- coxph(form1, data = bw_sample, ties = "breslow", cluster = id)
toc() # 103,647 sec elapsed

saveRDS(mod1, file = "Output/DLNM/Model_PM25_PTB.rds")

pred1 <- crosspred(cb1, mod1, cen=0) # pm25_krg_week_iqr
plot(pred1, ptype = "contour")

plot(pred1, ptype = "slice", var = 1,
    ylab="HR (95% CI)", xlab="Gestational Weeks", main="")
mtext(expression(bold("A. Overall - PM" [2.5])),  side = 3, adj = 0, font = 2, line = 1)


cb2 <- crossbasis(bw_sample$o3_krg_week_iqr, 
                   lag = c(0, 36),
                   #argvar = list(fun = "lin"),
                   argvar = list(fun = "bs", degree = 3, 
                                 knots = q_o3, #knots = q_o3,
                                 Boundary.knots = range_o3
                                 ), # Cubic 
                   arglag = list(fun = "ns", knots = lagknots)) # the form of the curve across all the lags, that is, the lag constrain 


form2 <- as.formula(paste("Surv(tstart, tstop, ", "event", ") ~ ", "cb2", 
                              "+ sex + age_group_mom + educ_group_mom + job_group_mom +",
                              "age_group_dad + educ_group_dad + job_group_dad +",
                              "factor(month_week1) + factor(year_week1) + vulnerability"))


tic()
mod2 <- coxph(form2, data = bw_sample, ties = "efron", cluster = id)
toc()

pred2 <- crosspred(cb2, mod2, cen=0) # o3_krg_week_iqr

plot(pred2, ptype = "slice", var = 1,
    ylab="HR (95% CI)", xlab="Gestational Weeks", main="")
mtext(expression(bold("B. Overall - O" [3])),  side = 3, adj = 0, font = 2, line = 1)
