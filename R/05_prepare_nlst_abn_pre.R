# Libraries ####
library(dplyr)
library(here)
library(readr)
library(survival)
library(purrr)
library(tidyr)
# Read in the Kovalchik prediction function
source(here("R/kovalchik.R"))

# Data Files ####
plco_control <- read_rds(here("data/plco_control.rds"))
lcrat <- read_rds(here("models/lcrat.rds"))
cox_death <- read_rds(here("models/cox_death.rds"))
# Provided by Optellum
lcp <- read_csv(here("data/nov17_max_lcp_for_nih.csv"))
nlst_abn_lrads <- read_rds(here("data/nlst_abn_lrads.rds"))

# Add prescreen risk ####
nlst_abn_lrads_pre <- nlst_abn_lrads %>%
    # mutate(log_diam = log(longest_diam + 1)) %>%
    # Calculate pre-screening risk inside this dataset
    mutate(prescr_1yrisk = risk.kovalchik(0, 1, ., lcrat, cox_death)) %>%
    mutate(log1yrisk = log(prescr_1yrisk),
           logit1yrisk = log(prescr_1yrisk / (1 - prescr_1yrisk)))

write_csv(nlst_abn_lrads_pre, here("data/nlst_abn_pre.csv"))
write_rds(nlst_abn_lrads_pre, here("data/nlst_abn_pre.rds"))

# Merge Optellum LCP score with NLST #### 
lcp <- lcp %>%
    drop_na(max_lcp_score) %>%
    rename(interval = yr) %>%
    inner_join(nlst_abn_lrads_pre, by = c("pid", "interval"))

write_csv(lcp, here("data/nlst_abn_pre_lcp.csv"))
write_rds(lcp, here("data/nlst_abn_pre_lcp.rds"))

# Emphysema #### 
split_pid <- function(filename) {
    readr::read_csv(here(filename)) %>%
        dplyr::distinct() %>%
        tidyr::separate(pid,
                        into = c("pid", "interval"),
                        sep = "/T") %>%
        dplyr::mutate(pid = as.numeric(pid),
                      interval = as.numeric(interval))
}

emph0 <- readr::read_csv(here("data/T0_data.csv")) %>%
    tibble::add_column(interval = 0) %>%
    dplyr::distinct()

emph <- bind_rows(emph0,
                  map_dfr(c("data/T1_data.csv",
                            "data/T2_data.csv"),
                          split_pid)) %>%
    group_by(pid, interval) %>%
    mutate(p_emph = mean(p_emph),
           p_no_emph = mean(p_no_emph)) %>%
    distinct()  %>%
    drop_na(p_emph) %>%
    left_join(nlst_abn_lrads,
              by = c("pid", "interval"))

write_csv(emph, here("data/nlst_abn_emph.csv"))
write_rds(emph, here("data/nlst_abn_emph.rds"))

# emph_pre <- emph %>%
#     mutate(prescr_1yrisk = risk.kovalchik(0, 1, ., lcrat, cox_death)) %>%
#     mutate(log1yrisk = log(prescr_1yrisk),
#            logit1yrisk = log(prescr_1yrisk / (1 - prescr_1yrisk)))

# emph_pre %>% write_rds(here("data/nlst_abn_emph_pre.rds")
# emph_pre %>% write_csv(here("data/nlst_abn_emph_pre.csv"))
