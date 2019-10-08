library(here)
library(h2o)

h2o.init()

nlst_abn_lcp <- h2o.importFile(here("data/nlst_abn_lcp.csv"))
nlst_abn_lcp$case <- as.factor(nlst_abn_lcp$case)
aml <- h2o.automl(y = "case", training_frame = nlst_abn_lcp)

lb <- aml@leaderboard

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
