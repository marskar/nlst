# List packages to be loaded (and installed if needed)
packages <-
    c(
        "lmtest",
        "here",
        "dplyr",
        "ggplot2",
        "survival",
        "gmodels",
        "devtools",
        "geepack",
        "MESS",
        "psych",
        "Hmisc",
        "glmnet",
        "boot"
    )

# List packages that are not installed
not_installed <-
    packages[!(packages %in% installed.packages()[, "Package"])]

# Install packages that are not installed
if (length(not_installed))
    install.packages(not_installed)

# Load all packages
lapply(packages, require, character.only = TRUE)

# Read in the Kovalchik prediction function
source(here("R/kovalchik.R"))

devtools::install_github('marskar/coxph_risk')
devtools::install_github('marskar/lcmodels')

# Function to create "not in" operator
'%!in%' <- function(x, y)
    ! ('%in%'(x, y))

# Control arm of PLCO who had no chest xray
plco_control <- readRDS('data/plco.rds') %>% 
    subset(control.group == 1) %>%
    # In the PLCO dataset, impute missing family history values to 0
    mutate(fam.lung.trend = ifelse(is.na(fam.lung.trend),
                                   0,
                                   plco$fam.lung.trend))

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
    saveRDS(file = 'models/lcrat.rds')

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
    saveRDS(file = 'models/cox_death.rds')
