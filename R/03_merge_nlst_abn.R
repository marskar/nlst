# Libraries ####
library(dplyr)
library(here)
library(readr)
library(tidyr)

# Data Files ####

# Created by script 00
nlst <- read_rds(here("data/nlst.rds"))
# created by script 02
abn_lrads <- read_rds(here("data/abn_lrads_merged.rds"))

# Pivot longer NLST and join with lrads ####

nlst_abn_lrads <- nlst %>%
    filter(screen_group == "CT") %>%
    pivot_longer(
        cols = starts_with("truefalse_scrnres_ly"),
        names_to = "interval",
        values_to = "screen_result",
        names_prefix = "truefalse_scrnres_ly",
        values_drop_na = TRUE
    ) %>%
    mutate(
        # add binary column cancer detected 1 year later (for optellum)
        case_at_next_screen = if_else(
            (interval != 2) &
            (lead(pid) == pid)
            & (lead(screen_result) %in% c(1, 6)),
                                      1, 0),
        # add binary column cancer detected at this screen (for immediate risk)
        case_at_current_screen = if_else(screen_result %in% c(1, 6), 1, 0)
    ) %>%
    mutate(interval = as.numeric(interval)) %>%
    left_join(abn_lrads, by = c("pid", "interval"))

# Data exploration ####

# total rows with case == 1: 2370
sum(nlst_abn_lrads$case)
sum(nlst_abn_lrads %>% select(case) %>% filter(case == 1))
# total cases: 1074 pids
nlst_abn_lrads %>% 
    group_by(pid) %>% 
    tally(case) %>% 
    filter(n > 0) %>% 
    tally()

# total cases at next screen: 401 pids
sum(nlst_abn_lrads$case_at_next_screen)
nlst_abn_lrads %>% 
    group_by(pid) %>% 
    tally(case_at_next_screen) %>% 
    filter(n > 0) %>% 
    tally()

# total cases at current screen: 692 pids
sum(nlst_abn_lrads$case_at_current_screen)
nlst_abn_lrads %>% 
    group_by(pid) %>% 
    tally(case_at_current_screen) %>% 
    filter(n > 0) %>% 
    tally()

write_rds(nlst_abn_lrads, here("data/nlst_abn_lrads.rds"))
dim(nlst_abn_lrads)

# total pids = 26408
nlst_abn_lrads %>% 
    group_by(pid) %>% 
    tally() %>% tally()
