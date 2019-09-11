rds <- readRDS("data/abn_lrads_lag_prescr.rds")

write.csv(rds, file = "data/abn_lrads_lag_prescr.csv")
