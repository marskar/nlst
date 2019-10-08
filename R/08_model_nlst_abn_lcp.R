library(here)
library(h2o)

h2o.init()

nlst_abn_lcp <- h2o.importFile(here("data/nlst_abn_lcp.csv"))
nlst_abn_lcp$case <- as.factor(nlst_abn_lcp$case)
aml <- h2o.automl(y = "case", training_frame = nlst_abn_lcp)

lb <- aml@leaderboard