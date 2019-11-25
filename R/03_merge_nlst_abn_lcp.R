library(dplyr)
library(here)
library(readr)
library(tidyr)
library(purrr)

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
    mutate(interval = as.numeric(interval)) %>%
    mutate(screen_result = if_else(screen_result %in% c(1, 2, 5, 6), 1, 0)) %>%
    left_join(abn_lrads, by = c("pid", "interval"))

# Optellum LCP score
lcp <- lcp %>%
    drop_na(max_lcp_score) %>% 
    rename(interval = yr) %>%
    inner_join(nlst_abn_lrads, by = c("pid", "interval"))

# Emphysema
split_pid <- function(filename) {
    readr::read_csv(here(filename)) %>%
    dplyr::distinct() %>% 
    tidyr::separate(pid,
                    into=c("pid", "interval"),
                    sep = "/T") %>% 
    dplyr::mutate(pid = as.numeric(pid),
                  interval = as.numeric(interval))
}

emph0 <- readr::read_csv(here("data/T0_data.csv")) %>%
    tibble::add_column(interval=0) %>%
    dplyr::distinct()

emph <- bind_rows(
    emph0,
    map_dfr(c("data/T1_data.csv",
              "data/T2_data.csv"),
            split_pid)
    ) %>% 
    group_by(pid, interval) %>% 
    mutate(p_emph = mean(p_emph),
           p_no_emph = mean(p_no_emph)) %>% 
    distinct()  %>%
    drop_na(p_emph) %>%
    left_join(nlst_abn_lrads,
               by = c("pid", "interval"))

write_rds(nlst_abn_lrads, here("data/nlst_abn_lrads.rds"))
write_csv(lcp, here("data/nlst_abn_lcp.csv"))
write_rds(lcp, here("data/nlst_abn_lcp.rds"))
write_csv(emph, here("data/nlst_abn_emph.csv"))
write_rds(emph, here("data/nlst_abn_emph.rds"))
