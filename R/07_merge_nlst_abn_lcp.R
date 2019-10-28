library(here)
library(dplyr)
library(readr)

# di created by script 6
di <- read_rds(here("data/data_interval_abn.rds"))
lcp <- read_csv(here("data/nov17_max_lcp_for_nih.csv"))

lcp <- lcp %>%
    drop_na(max_lcp_score) %>% 
    mutate(interval = yr+1) %>%
    inner_join(di, by = c("pid", "interval"))

write_csv(lcp, "data/nlst_abn_lcp.csv")
write_rds(lcp, "data/nlst_abn_lcp.rds")
