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


nlst <-
    # Load NLST data
    readRDS('data/nlst.rds') %>%
    # Remove people with <30 pack-years and age<55 or age>74 from NLST
    filter(pkyears.cat != "[0,30)"
           & age >= 55
           & age <= 74) %>%
    # Drop levels to remove pack-years variable [0,30) empty level
    droplevels() %>%
    # Make a variable for days to diagnosis
    mutate(days_to_dx = ifelse(case == 1, 365 * incidence.years, NA))

# Make a subset of NLST data with the LCRAT variables that we will need later
varlist <-
    c(
        "female",
        "race",
        "edu6",
        "fam.lung.trend",
        "emp",
        "bmi",
        "cpd",
        "pkyears.cat",
        "age",
        "qtyears",
        "smkyears"
    )

nlst %>%
    select(varlist, pid, lss) %>%
    saveRDS(file = 'data/nlst_sub.rds')


# Subset to CT arm in NLST and make a pos/neg variable for the first, second, and third screens
nlst %>%
    mutate(T0posneg = case_when(
        truefalse_scrnres_ly0 %in% c(4, 5, 6) ~ 0,
        truefalse_scrnres_ly0 %in% c(1, 2, 3) ~ 1
    )) %>%
    mutate(T1posneg = case_when(
        truefalse_scrnres_ly1 %in% c(4, 5, 6) ~ 0,
        truefalse_scrnres_ly1 %in% c(1, 2, 3) ~ 1
    )) %>%
    mutate(T2posneg = case_when(
        truefalse_scrnres_ly2 %in% c(4, 5, 6) ~ 0,
        truefalse_scrnres_ly2 %in% c(1, 2, 3) ~ 1
    )) %>%
    # Subset to CT arm and create screening history variables
    subset(screen_group == "CT") %>%
    mutate(
        hist.T0.T1 =
            case_when(
                T0posneg == 0 & T1posneg == 0 ~ 1,
                T0posneg == 0 & T1posneg == 1 ~ 2,
                T0posneg == 1 & T1posneg == 0 ~ 3,
                T0posneg == 1 & T1posneg == 1 ~ 4
            )
    ) %>%
    mutate(
        hist.T1.T2 =
            case_when(
                T1posneg == 0 & T2posneg == 0 ~ 1,
                T1posneg == 0 & T2posneg == 1 ~ 2,
                T1posneg == 1 & T2posneg == 0 ~ 3,
                T1posneg == 1 & T2posneg == 1 ~ 4
            )
    ) %>%
    mutate_at(vars(hist.T0.T1, hist.T1.T2),
              funs(factor(
                  .,
                  levels = c(1, 2, 3, 4),
                  labels = c("Neg-Neg", "Neg-Pos", "Pos-Neg", "Pos-Pos")
              ))) %>%
    saveRDS(file = 'data/nlst_ct.rds')

nlst_ct <- readRDS(file = 'data/nlst_ct.rds')

### Create dataset for risk from T0 to T1.
# 0 inadeq, 1 true-pos, 2 poss true-pos, 3 false-pos, 4 true-neg, 5 poss false-neg, 6 false-neg
# At risk for screen-detected at T1: either false-positive or true-negative at T0, and did not have any of the following at T1:
# inadequate image, left study, refused, wrong screen, erroneous report of LC, form not submitted (no missing values of scr_res0).
# Case status: case=1 AND either of (true-pos at T1 or T1 is coded as "not expected: cancer/death in screening window")

nlst_ct %>%
    filter(truefalse_scrnres_ly0 %in% c(2, 3, 4, 5)
           & scr_res1 %!in% c(10, 11, 15, 17, 95, 97)) %>%
    mutate(case_T1_screen = if_else((
        truefalse_scrnres_ly1 == 1 |
            scr_res1 %in% c(23, 24)
    )
    & case == 1,
    1,
    0)) %>%
    saveRDS(file = 'data/nlst_ct_t1_scrisk.rds')

### Create dataset for risk from T1 to T2.
# At risk for screen-detected at T2: either false-positive or true-negative at T1, and did not have any of the following at T2:
# inadequate image, left study, refused, wrong screen, erroneous report of LC, form not submitted (no missing values of scr_res0).
# Case status: case=1 AND either of (true-pos at T2 or T2 is coded as "not expected: cancer/death in screening window")

nlst_ct %>%
    filter(truefalse_scrnres_ly1 %in% c(2, 3, 4, 5)
           & scr_res2 %!in% c(10, 11, 15, 17, 95, 97)) %>%
    mutate(case_T2_screen = if_else((
        truefalse_scrnres_ly2 == 1
        | scr_res2 %in% c(23, 24)
    ) & case == 1,
    1,
    0)) %>%
    saveRDS(file = 'data/nlst_ct_t2_scrisk.rds')

### Create datasets for risk from T0 to T1. One for interval, one for screen-detected.
## .neg intended for analysis of interval cancer risk; .scrisk intended for analysis of screen-detected cancers
# 0 inadeq, 1 true-pos, 2 poss true-pos, 3 false-pos, 4 true-neg, 5 poss false-neg, 6 false-neg
# At risk in interval: T0 negatives only. Case status: false-negative at T0.
nlst_ct %>%
    filter(truefalse_scrnres_ly0 %in% c(4, 5, 6)) %>%
    mutate(case_T1_interval = ifelse(truefalse_scrnres_ly0 == 6, 1, 0)) %>%
    saveRDS(file = 'data/nlst_ct_t1_neg.rds')

# At risk in interval: T1 negatives only. Case status: false-negative at T1.
nlst_ct %>%
    filter(truefalse_scrnres_ly1 %in% c(4, 5, 6)) %>%
    mutate(case_T2_interval = ifelse(truefalse_scrnres_ly1 == 6, 1, 0)) %>%
    saveRDS(file = 'data/nlst_ct_t2_neg.rds')
### Create a dataset for risk during "interval" after T2 (within 1 year)
nlst_ct %>%
    filter(truefalse_scrnres_ly2 %in% c(4, 5, 6)) %>%
    mutate(case_T3_interval = ifelse(truefalse_scrnres_ly2 == 6, 1, 0)) %>%
    saveRDS(file = 'data/nlst_ct_t3_neg.rds')
