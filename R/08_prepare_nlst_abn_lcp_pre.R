# Read in the Kovalchik prediction function
library(survival)
library(glmnet)
library(dplyr)
library(stats)
library(here)
library(readr)
library(tidyr)
source(here("R/kovalchik.R"))
plco <- readRDS(here("data/plco.rds"))
lcp <- read_rds(here("data/nlst_abn_lcp.rds"))

t0_nlst_emp <- read_csv(here('data/T0_data.csv'))
# t1_nlst_emp <- read_csv(here('data/T1_data.csv')) %>% 
#     mutate(pid = as.numeric(pid))
# t2_nlst_emp <- read_csv(here('data/T2_data.csv')) %>% 
#     mutate(pid = as.numeric(pid))

# Control arm of PLCO who had no chest xray
plco_control <- subset(plco, control.group == 1)

LCRAT <- coxph(
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
)

cox_death <- coxph(
    Surv(years.followed, other.cause.death)
    ~ female
    + race
    + edu6
    + emp
    + I(bmi <= 18.5)
    + I(cpd > 20)
    + as.factor(pkyears.cat)
    + I((age) ^ 2)
    + I((bmi - 25) ^ 2)
    + I(log(qtyears + 1))
    + smkyears,
    data = plco_control
)


lcp_pre <- lcp %>%
    # mutate(log_diam = log(longest_diam + 1)) %>%
    # Calculate pre-screening risk inside this dataset
    mutate(prescr_1yrisk = risk.kovalchik(0,
                                          1,
                                          .,
                                          LCRAT,
                                          cox_death)) %>%
    mutate(log1yrisk = log(prescr_1yrisk),
           logit1yrisk = log(prescr_1yrisk / (1 - prescr_1yrisk))) %>%
    mutate(
        diam.cat = case_when(
            longest.diam == 0 ~ 1,
            longest.diam > 0 & longest.diam <= 5 ~ 2,
            longest.diam > 5 & longest.diam <= 7 ~ 3,
            longest.diam > 7 & longest.diam <= 10 ~ 4,
            longest.diam > 10 & longest.diam <= 13 ~ 5,
            longest.diam > 13 & longest.diam < 100 ~ 6
        )
    )
    

# lcp_pre_emph <- lcp_pre %>% 
#     left_join(t0_nlst_emp, by = "pid") %>%
#     identity()
    # full_join(t1_nlst_emp, by = "pid") %>% 
    # full_join(t2_nlst_emp, by = "pid")

lcp_pre %>% write_rds(here("data/nlst_abn_lcp_pre.rds"))
lcp_pre %>% write_csv(here("data/nlst_abn_lcp_pre.csv"))
