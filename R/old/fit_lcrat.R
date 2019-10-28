#### Setup ####

packages <- c(
    "dplyr",
    "tidyr",
    "ggplot2",
    "readr",
    "here",
    "survival",
    "gmodels",
    "coxph.risk",
    "geepack",
    "MESS",
    "psych",
    "pROC",
    "Hmisc",
    "glmnet",
    "boot"
)

not_installed <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(not_installed)) install.packages(not_installed)
lapply(packages, require, character.only = TRUE)

# Read in the Kovalchik prediction function
source(here::here("R/kovalchik.R"))

# Install coxph.risk package
devtools::install_github('marskar/coxph_risk')
# Install lcmodels package
devtools::install_github('marskar/lcmodels')

# Function to create "not in" operator
'%!in%' <- function(x,y)!('%in%'(x,y))

# Load PLCO and NLST data sets.
# New versions provided by Li Cheung
# On 11 July 2016, 14 July 2016, 20 July 2016, 8 Aug 2016.
# load("hilary.RData")  # Load NLST and PLCO data
nlst = readRDS('data/nlst.rds')
plco = readRDS('data/plco.rds')
# In the PLCO dataset, impute missing family history values to 0
plco$fam_lung_trend <- ifelse(is.na(plco$fam.lung.trend), 0, plco$fam.lung.trend)
# Control arm of PLCO who had no chest xray
plco_control <- subset(plco, control.group == 1)

# Remove people with <30 pack-years and age<55 or age>74 from NLST
nlst <- nlst %>%
    filter(pkyears.cat != "[0,30)" & age >= 55 & age <= 74) %>%
    # Make a new pack-years variable to get rid of the [0,30) level
    mutate(pkyears.cat.clone = ifelse(
        pkyears.cat == "[30,40)",
        "[30,40)",
        ifelse(
            pkyears.cat == "[40,50)",
            "[40,50)",
            ifelse(pkyears.cat == "[50,Inf]", "[50,Inf]", NA)
        )
    )) %>%
    mutate(pkyears.cat = as.factor(pkyears.cat.clone)) %>%
    # Make a variable for days to diagnosis
    mutate(days_to_dx = ifelse(case == 1, 365 * incidence.years, NA)) %>%
    identity()

# Make a subset of NLST data with the LCRAT variables that we will need later to merge back with at-risk datasets
varlist <-
    c(
        "age",
        "female",
        "smkyears",
        "qtyears",
        "cpd",
        "race",
        "emp",
        "fam.lung.trend",
        "bmi",
        "edu6"
    )

lcratvars <-
    c(
        "age",
        "female",
        "smkyears",
        "qtyears",
        "cpd",
        "race",
        "emp",
        "fam.lung.trend",
        "bmi",
        "edu6",
        "pkyears.cat"
    )

nlst_lcrat <- as.data.frame(cbind(nlst[, lcratvars], pid=nlst$pid, lss=as.numeric(nlst$lss)))

varlist <-
    c(
        "age",
        "female",
        "smkyears",
        "qtyears",
        "cpd",
        "race",
        "emp"
    )
nlst_sub <- nlst %>% select(varlist)
nlst.sub <- as.data.frame(cbind(nlst[, varlist], pid = nlst$pid, lss = as.numeric(nlst$lss)))
lcmodels::lcmodels(nlst_sub)
View()
lcmodels::lcmodels(lcmodels::LCRAT)

names(lcmodels::LCRAT$coefficients)
