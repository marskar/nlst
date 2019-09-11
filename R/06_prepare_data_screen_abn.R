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


# Load abnormalities data (person-screen level) and merge with data_screen
# This dataset was created by 02_prepare_abn_data.R script
abn <- readRDS(here("data/abn_lrads_merged.rds"))
plco_control <- readRDS(here("data/plco_control.rds"))
# This dataset was created by 03_prepare_data_screen.R script
data_screen <- readRDS(here("data/data_screen.rds"))
lcrat <- readRDS(here("models/lcrat.rds"))
cox_death <- readRDS(here("models/cox_death.rds"))

data_screen_abn <- merge(
    data_screen,
    abn,
    by = c("pid", "interval"),
    all.x = TRUE,
    all.y = FALSE
) %>%
    # Replace NAs with 0 (not present) for appropriate variables
    mutate(
        pid = ifelse(is.na(pid), 0, pid),
        interval = ifelse(is.na(interval), 0, interval),
        LRcat  = ifelse(is.na(LRcat), 0, LRcat)
    ) %>%
    mutate(longest_diam = ifelse(is.infinite(longest_diam), 0, longest_diam)) %>%
    # Make a variable for including observations in Lung-RADS analysis
    mutate(LR_include = (LRcat %in% c("1", "2", "3", "4A", "4B", "4X"))) %>%
    # Create variable for log(diameter)
    mutate(log_diam = log(longest_diam + 1)) %>%
    # Calculate pre-screening risk inside this dataset
    mutate(prescr_1yrisk = risk.kovalchik(0, 1, ., lcrat, cox_death)) %>%
    mutate(log1yrisk = log(prescr_1yrisk),
           logit1yrisk = log(prescr_1yrisk / (1 - prescr_1yrisk))) %>%
    saveRDS(here('data/data_screen_abn.rds'))

# This datasets is needed to separately model screen-detected cancers incorporating abnormalities for false-positives
filter(data_screen_abn, screen_result == 0) %>%
    saveRDS(here('data/data_screen_abn_neg.rds'))

filter(data_screen_abn, screen_result == 1) %>%
    saveRDS(here('data/data_screen_abn_pos.rds'))
