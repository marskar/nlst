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


# Load PLCO and NLST data sets. New versions provided by Li Cheung on 11 July 2016, 14 July 2016, 20 July 2016, 8 Aug 2016.
# Load NLST and PLCO data
nlst = readRDS('data/nlst.rds')
plco = readRDS('data/plco.rds')

# In the PLCO dataset, impute missing family history values to 0
plco$fam.lung.trend <-
    ifelse(is.na(plco$fam.lung.trend),
           0,
           plco$fam.lung.trend)

# Control arm of PLCO who had no chest xray
plco_control <- subset(plco, control.group == 1)

# Remove people with <30 pack-years and age<55 or age>74 from NLST
nlst <-
    nlst %>%
    filter(pkyears.cat != "[0,30)"
           & age >= 55
           & age <= 74) %>%
    # Drop levels to remove pack-years variable [0,30) empty level
    droplevels() %>%
    # Make a variable for days to diagnosis
    mutate(days_to_dx = ifelse(case == 1, 365 * incidence.years, NA)) %>%
    identity()

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

nlst_sub <- nlst %>% select(varlist, pid, lss)

# To later calculate pre-screening risk, we must first fit the incidence model and other-cause death models in PLCO.
LCRAT <-
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
    )

cox.death <-
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
    )

# Subset to CT arm in NLST and make a pos/neg variable for the first, second, and third screens
nlst <-
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
    ))

# Subset to CT arm and create screening history variables
nlst_CT <-
    nlst %>%
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
    identity()

### Create dataset for risk from T0 to T1.
# 0 inadeq, 1 true-pos, 2 poss true-pos, 3 false-pos, 4 true-neg, 5 poss false-neg, 6 false-neg
# At risk for screen-detected at T1: either false-positive or true-negative at T0, and did not have any of the following at T1:
# inadequate image, left study, refused, wrong screen, erroneous report of LC, form not submitted (no missing values of scr_res0).
# Case status: case=1 AND either of (true-pos at T1 or T1 is coded as "not expected: cancer/death in screening window")

nlst_CT_T1_scrisk <-
    nlst_CT %>%
    filter(truefalse_scrnres_ly0 %in% c(2, 3, 4, 5)
           & scr_res1 %!in% c(10, 11, 15, 17, 95, 97)) %>%
    mutate(case_T1_screen = if_else((
        truefalse_scrnres_ly1 == 1 |
            scr_res1 %in% c(23, 24)
    )
    & case == 1,
    1,
    0))

### Create dataset for risk from T1 to T2.
# At risk for screen-detected at T2: either false-positive or true-negative at T1, and did not have any of the following at T2:
# inadequate image, left study, refused, wrong screen, erroneous report of LC, form not submitted (no missing values of scr_res0).
# Case status: case=1 AND either of (true-pos at T2 or T2 is coded as "not expected: cancer/death in screening window")

nlst_CT_T2_scrisk <-
    nlst_CT %>%
    filter(truefalse_scrnres_ly1 %in% c(2, 3, 4, 5)
           & scr_res2 %!in% c(10, 11, 15, 17, 95, 97)) %>%
    mutate(case_T2_screen = if_else((
        truefalse_scrnres_ly2 == 1
        | scr_res2 %in% c(23, 24)
    ) & case == 1,
    1,
    0))

# Construct dataset to model risk of ALL screen-detected cancers (at T1 and T2)
data_screen <-
    data.frame(
        pid = c(nlst_CT_T1_scrisk$pid, nlst_CT_T2_scrisk$pid),
        case = c(
            nlst_CT_T1_scrisk$case_T1_screen,
            nlst_CT_T2_scrisk$case_T2_screen
        ),
        screen_result = c(nlst_CT_T1_scrisk$T0posneg,
                          nlst_CT_T2_scrisk$T1posneg),
        interval = c(rep(1, times = nrow(nlst_CT_T1_scrisk)),
                     rep(2, times = nrow(nlst_CT_T2_scrisk)))
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
