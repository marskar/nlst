library(dplyr)
library(here)
library(readr)
library(survival)
# Read in the Kovalchik prediction function
source(here("R/kovalchik.R"))

# Control arm of PLCO who had no chest xray
plco_control <- readRDS(here('data/plco.rds')) %>%
    subset(control.group == 1) %>%
    identity()
    # In the PLCO dataset,
    # impute missing family history values to 0
    # mutate(fam_lung_trend = ifelse(is.na(fam.lung.trend),
    #                                 0,
    #                                 fam.lung.trend))

write_rds(plco_control, here("data/plco_control.rds"))

# To later calculate pre-screening risk, we must first fit the incidence model and other-cause death models in PLCO.
coxph(
    Surv(incidence.years, case)
    ~ female
    + race
    + edu6
    + fam.lung.trend
    + emp
    + I(bmi <= 18.5)
    + I(cpd > 20)
    + as.factor(pkyears.cat)
    + I(log(age))
    + I(log(bmi))
    + I(log(qtyears + 1))
    + smkyears,
    data = plco_control
) %>%
    write_rds('models/lcrat.rds')

coxph(
    Surv(years.followed, other.cause.death)
    ~ female
    + race
    + edu6
    + emp
    + I(bmi <= 18.5)
    + I(cpd > 20)
    + as.factor(pkyears.cat)
    + I((age) ^ 2) + I((bmi - 25) ^ 2)
    + I(log(qtyears + 1))
    + smkyears,
    data = plco_control
) %>%
    write_rds('models/cox_death.rds')
