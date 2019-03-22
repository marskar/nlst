# Convert Hilary.Rdata to RDS
# install.packages("here")
library(here)

load(here("data/hilary.RData"))  # Load NLST and PLCO data
saveRDS(nlst, file = "data/nlst.rds")
saveRDS(plco, file = "data/plco.rds")
