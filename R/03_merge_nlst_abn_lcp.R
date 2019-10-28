library(dplyr)
library(here)
library(readr)
library(tidyr)

# Created by script 00
nlst <- read_rds(here("data/nlst.rds"))
# created by script 02
abn_lrads <- read_rds(here("data/abn_lrads_merged.rds"))
# Provided by Optellum
lcp <- read_csv(here("data/nov17_max_lcp_for_nih.csv"))

nlst_abn_lrads <- nlst %>%
    filter(screen_group == "CT") %>%
    pivot_longer(
        cols = starts_with("truefalse_scrnres_ly"),
        names_to = "interval",
        values_to = "screen_result",
        names_prefix = "truefalse_scrnres_ly",
        values_drop_na = TRUE
    ) %>%
    mutate(interval = is.numeric(interval)) %>%
    mutate(screen_result = if_else(screen_result %in% c(1, 2, 5, 6), 1, 0)) %>%
    left_join(abn_lrads, by = c("pid", "interval")) %>% 
    mutate(
        diam.cat = case_when(
            longest_diam == 0 ~ 1,
            longest_diam > 0 & longest_diam <= 5 ~ 2,
            longest_diam > 5 & longest_diam <= 7 ~ 3,
            longest_diam > 7 & longest_diam <= 10 ~ 4,
            longest_diam > 10 & longest_diam <= 13 ~ 5,
            longest_diam > 13 & longest_diam < 100 ~ 6
        )
    )
    

nlst_abn_lrads %>% write_rds(here("data/nlst_abn_lrads.rds"))

lcp <- lcp %>%
    drop_na(max_lcp_score) %>% 
    mutate(interval = yr) %>%
    inner_join(nlst_abn_lrads, by = c("pid", "interval"))

write_csv(lcp, "data/nlst_abn_lcp.csv")
write_rds(lcp, "data/nlst_abn_lcp.rds")
