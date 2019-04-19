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
nlst_ct_t1_neg = readRDS('data/nlst_ct_t1_neg.rds')
nlst_ct_t2_neg = readRDS('data/nlst_ct_t2_neg.rds')
nlst_ct_t3_neg = readRDS('data/nlst_ct_t3_neg.rds')
nlst_sub = readRDS('data/nlst_sub.rds')
nlst_ct = readRDS('data/nlst_ct.rds')

# Load abnormalities data (person-screen level) and merge with data_screen
# This dataset was created by 02_prepare_abn_data.R script
abn <- readRDS(here("data/abn_lrads_merged.rds"))
plco_control <- readRDS(here("data/plco_control.rds"))
lcrat <- readRDS(here("models/lcrat.rds"))
cox_death <- readRDS(here("models/cox_death.rds"))


### Construct dataset to model risk of ALL interval cancers across all 3 screens
# Interval variable in data.interval datasets is 1 for T0-T1, 2 for T1-T2, and 3 for post-T2 intervals
data_interval <- data.frame(
    pid = c(nlst_ct_t1_neg$pid, nlst_ct_t2_neg$pid, nlst_ct_t3_neg$pid),
    case = c(
        nlst_ct_t1_neg$case_T1_interval,
        nlst_ct_t2_neg$case_T2_interval,
        nlst_ct_t3_neg$case_T3_interval
    ),
    interval = c(rep(1, times = nrow(nlst_ct_t1_neg)),
                 rep(2, times = nrow(nlst_ct_t2_neg)),
                 rep(3, times = nrow(nlst_ct_t3_neg)))
) %>%
    # Merge this back with covariates from NLST. Add screening history variable.
    merge(nlst_sub, by = "pid", all.x = T) %>%
    merge(
        select(nlst_ct, pid, hist.T0.T1, hist.T1.T2),
        by = "pid",
        all.x = T,
        all.y = F
    ) %>%
    mutate(screen.hist = ifelse(interval == 2, hist.T0.T1, ifelse(interval ==
                                                                      3, hist.T1.T2, NA))) %>%
    mutate(screen.hist = factor(
        screen.hist,
        levels = c(1, 2, 3, 4),
        labels = c("Neg-Neg", "Neg-Pos", "Pos-Neg", "Pos-Pos")
    )) %>%
    select(-c(hist.T0.T1, hist.T1.T2)) %>%
    arrange(pid, interval) %>%
    # Update age, quit-years, and smoke-years by adding a year for T1 and T2
    mutate(
        age = ifelse(interval == 2, age + 1, ifelse(interval == 3, age + 2, age)),
        smkyears = ifelse(
            interval == 2 &
                qtyears == 0,
            smkyears + 1,
            ifelse(interval == 3 &
                       qtyears == 0, smkyears + 2, smkyears)
        ),
        qtyears = ifelse(
            interval == 2 &
                qtyears > 0,
            qtyears + 1,
            ifelse(interval == 3 &
                       qtyears > 0, qtyears + 2, qtyears)
        )
    ) %>%
    # Using new smoke-years, update pack-years, then re-categorize
    mutate(pkyears.cont = cpd * smkyears / 20) %>%
    mutate(pkyears.cat = as.factor(ifelse(
        pkyears.cont >= 30 & pkyears.cont < 40,
        "[30,40)",
        ifelse(
            pkyears.cont >= 40 &
                pkyears.cont < 50,
            "[40,50)",
            ifelse(pkyears.cont >= 50 &
                       pkyears.cont < 999, "[50,Inf]", NA)
        )
    )))

# Merge abnormalities data (person-screen level) with data.interval
merge(
    data_interval,
    abn,
    by = c("pid", "interval"),
    all.x = T,
    all.y = F
) %>%
    # Replace NAs with 0 (not present) for appropriate variables
    mutate(
        pid = ifelse(is.na(pid), 0, pid),
        interval = ifelse(is.na(interval), 0, interval),
        LRcat  = ifelse(is.na(LRcat), 0, LRcat)
    ) %>%
    mutate(longest.diam = ifelse(is.infinite(longest.diam), 0, longest.diam)) %>%
    # Make a variable for including observations in Lung-RADS analysis
    mutate(LR_include = (LRcat %in% c("1", "2", "3", "4A", "4B", "4X"))) %>%
    # Create variable for log(diameter)
    mutate(log_diam = log(longest.diam + 1)) %>%
    # Calculate pre-screening risk inside this dataset
    mutate(prescr_1yrisk = risk.kovalchik(0, 1, ., lcrat, cox_death)) %>%
    mutate(log1yrisk = log(prescr_1yrisk),
           logit1yrisk = log(prescr_1yrisk / (1 - prescr_1yrisk))) %>%
    saveRDS(here('data/data_interval_abn.rds'))
