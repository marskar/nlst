library(dplyr)
library(here)
library(readr)
library(survival)
library(purrr)
source(here("R/kovalchik.R"))
plco_control <- read_rds(here("data/plco_control.rds"))
lcp <- read_rds(here("data/nlst_abn_lcp.rds"))
emph <- read_rds(here("data/nlst_abn_emph.rds"))
lcrat <- read_rds(here("models/lcrat.rds"))
cox_death <- read_rds(here("models/cox_death.rds"))

lcp_pre <- lcp %>%
    # mutate(log_diam = log(longest_diam + 1)) %>%
    # Calculate pre-screening risk inside this dataset
    mutate(prescr_1yrisk = risk.kovalchik(0, 1, ., lcrat, cox_death)) %>%
    mutate(log1yrisk = log(prescr_1yrisk),
           logit1yrisk = log(prescr_1yrisk / (1 - prescr_1yrisk)))

emph_pre <- emph %>%
    mutate(prescr_1yrisk = risk.kovalchik(0, 1, ., lcrat, cox_death)) %>%
    mutate(log1yrisk = log(prescr_1yrisk),
           logit1yrisk = log(prescr_1yrisk / (1 - prescr_1yrisk)))

lcp_pre %>% write_rds(here("data/nlst_abn_lcp_pre.rds"))
lcp_pre %>% write_csv(here("data/nlst_abn_lcp_pre.csv"))

emph_pre %>% write_rds(here("data/nlst_abn_emph_pre.rds"))
emph_pre %>% write_csv(here("data/nlst_abn_emph_pre.csv"))
