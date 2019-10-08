library(here)
library(dplyr)
library(readr)

nlst <- read_rds(here("data/nlst.rds"))
abn <- read_rds(here("data/abn_lrads_merged.rds"))
lcp <- read_csv(here("data/nov17_max_lcp_for_nih.csv"))

lcp <- lcp %>%
    drop_na(max_lcp_score)

nlst_abn <- nlst %>%
    filter(screen_group == "CT") %>%
    merge(abn, by = "pid") %>%
    merge(lcp, by = "pid")

write_csv(nlst_abn, "data/nlst_abn_lcp.csv")
