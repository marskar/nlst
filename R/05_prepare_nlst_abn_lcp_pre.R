library(dplyr)
library(here)
library(readr)
plco <- read_rds(here("data/plco.rds"))
lcp <- read_rds(here("data/nlst_abn_lcp.rds"))
lcrat <- read_rds(here("models/lcrat.rds"))
cox_death <- read_rds(here("models/cox_death.rds"))

# t0_nlst_emp <- read_csv(here('data/T0_data.csv'))
# t1_nlst_emp <- read_csv(here('data/T1_data.csv')) %>%
#     mutate(pid = as.numeric(pid))
# t2_nlst_emp <- read_csv(here('data/T2_data.csv')) %>%
#     mutate(pid = as.numeric(pid))

lcp_pre <- lcp %>%
    # mutate(log_diam = log(longest_diam + 1)) %>%
    # Calculate pre-screening risk inside this dataset
    mutate(prescr_1yrisk = risk.kovalchik(0, 1, ., lcrat, cox_death)) %>%
    mutate(log1yrisk = log(prescr_1yrisk),
           logit1yrisk = log(prescr_1yrisk / (1 - prescr_1yrisk)))


# lcp_pre_emph <- lcp_pre %>%
#     left_join(t0_nlst_emp, by = "pid") %>%
#     identity()
# full_join(t1_nlst_emp, by = "pid") %>%
# full_join(t2_nlst_emp, by = "pid")

lcp_pre %>% write_rds(here("data/nlst_abn_lcp_pre.rds"))
lcp_pre %>% write_csv(here("data/nlst_abn_lcp_pre.csv"))
