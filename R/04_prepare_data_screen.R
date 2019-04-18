#### Setup ####

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


# Load NLST data
nlst_sub = readRDS('data/nlst_sub.rds')
nlst_ct_t1_scrisk = readRDS('data/nlst_ct_t1_scrisk.rds')
nlst_ct_t2_scrisk = readRDS('data/nlst_ct_t2_scrisk.rds')

# Construct dataset to model risk of ALL screen-detected cancers (at T1 and T2)
# Interval variable in data.screen datasets is 1 for risk at T1 and 2 for risk at T2
data_screen <-
    data.frame(
        pid = c(nlst_ct_t1_scrisk$pid, nlst_ct_t2_scrisk$pid),
        case = c(
            nlst_ct_t1_scrisk$case_T1_screen,
            nlst_ct_t2_scrisk$case_T2_screen
        ),
        screen_result = c(nlst_ct_t1_scrisk$T0posneg,
                          nlst_ct_t2_scrisk$T1posneg),
        interval = c(rep(1, times = nrow(nlst_ct_t1_scrisk)),
                     rep(2, times = nrow(nlst_ct_t2_scrisk)))
    ) %>%
    # Merge this back with covariates from NLST
    merge(nlst_sub, by = "pid", all.x = TRUE) %>%
    # Add a variable for lagged screen result & a 6-level variable for all combinations
    group_by(pid) %>%
    mutate(lag_screen = dplyr::lag(screen_result, order_by = interval)) %>%
    mutate(
        screen_comb = case_when(
            interval == 1 & screen_result == 0 ~ 1,
            interval == 1 & screen_result == 1 ~ 2,
            interval == 2 &
                lag_screen == 0 & screen_result == 0 ~ 3,
            interval == 2 &
                lag_screen == 0 & screen_result == 1 ~ 4,
            interval == 2 &
                lag_screen == 1 & screen_result == 0 ~ 5,
            interval == 2 & lag_screen == 1 & screen_result == 1 ~ 6
        )
    ) %>%
    mutate(screen_comb = factor(
        screen_comb,
        levels = c(1, 2, 3, 4, 5, 6),
        labels = c("Neg", "Pos", "Neg-Neg", "Neg-Pos", "Pos-Neg", "Pos-Pos")
    )) %>%
    # Update age, quit-years, and smoke-years by adding a year for T1
    mutate_at(c("age", "smkyears", "qtyears"), funs(as.numeric)) %>%
    mutate(
        age = ifelse(interval == 2, age + 1, age),
        smkyears = ifelse(interval == 2 &
                              qtyears == 0, smkyears + 1, smkyears),
        qtyears = ifelse(interval == 2 &
                             qtyears > 0, qtyears + 1, qtyears)
    ) %>%
    # Using new smoke-years, update pack-years, then re-categorize
    mutate(pkyears_cont = cpd * smkyears / 20) %>%
    mutate(
        pkyears_cat = case_when(
            pkyears_cont >= 30 & pkyears_cont < 40 ~ 1,
            pkyears_cont >= 40 & pkyears_cont < 50 ~ 2,
            pkyears_cont >= 50 & pkyears_cont < 999 ~ 3
        )
    ) %>%
    mutate(pkyears_cat = factor(
        pkyears_cat,
        levels = 1:3,
        labels = c("[30,40)", "[40,50)", "[50,Inf]")
    )) %>%
    saveRDS(file = "data/data_screen.rds")
