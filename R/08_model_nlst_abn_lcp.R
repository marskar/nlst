library(here)
library(dplyr)
library(tidyr)

nlst_abn_lcp <- read_rds(here("data/nlst_abn_lcp.rds"))

nlst_abn_lag %>%
    filter(screen_group == "CT") %>%
    merge(abn, by = "pid") %>%
    merge(lcp, by = "pid")

lag_vars <- function(.data, .group_var, .vars) {
    .data %>%
        group_by(!!enquo(.group_var)) %>%
        mutate_at(.vars, .funs = list(lag = ~ lag))
}

colname_vector <- names(nlst_abn)[35:68]
nlst_abn %>%
    lag_vars(pid, colname_vector) %>%
    lag_vars(pid, paste0(colname_vector, "_lag")) %>%
    saveRDS(file = "data/nlst_abn_lag.rds")
names(nlst_abn)
