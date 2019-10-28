# Convert Hilary.Rdata to RDS
# install.packages("here")
library(here)
library(readr)

# Load NLST and PLCO data
# Provided by Li Cheung
load(here("data/hilary.RData"))
write_rds(nlst, here("data/nlst.rds"))
write_rds(plco, here("data/plco.rds"))
