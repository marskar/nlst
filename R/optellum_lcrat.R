library(here)
library(survival)
library(dplyr)
source(here("R/kovalchik.R"))

nlst <- readRDS(here("data/nlst.rds"))
plco <- readRDS("data/plco.rds")

plco.control <-
    subset(plco, control.group == 1) # control arm of PLCO who had no chest xray

# TODO double check against lcmodels
LCRAT <- coxph(
    Surv(incidence.years, case) ~
        female +
        race +
        edu6 +
        fam.lung.trend +
        emp +
        I(bmi <= 18.5) +
        I(cpd > 20) +
        as.factor(pkyears.cat) +
        I(log(age)) +
        I(log(bmi)) +
        I(log(qtyears + 1)) +
        smkyears,
    data = plco.control
)

cox_death <- coxph(
    Surv(years.followed, other.cause.death) ~
        female +
        race +
        edu6 +
        emp +
        I(bmi <= 18.5) +
        I(cpd > 20) +
        as.factor(pkyears.cat) +
        I((age) ^ 2) +
        I((bmi - 25) ^ 2) +
        I(log(qtyears + 1)) +
        smkyears,
    data = plco.control
)

optellum <- nlst %>%
    filter(screen_group == "CT") %>%
    mutate(prescr_1yrisk = risk.kovalchik(0, 1, ., LCRAT, cox_death)) %>%
    select(pid, prescr_1yrisk)

readr::write_csv(optellum, path = "data/optellum_matrix.csv")
