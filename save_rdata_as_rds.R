# Convert Hilary.Rdata to RDS

load("hilary.RData")  # Load NLST and PLCO data
saveRDS(nlst, file = "nlst.rds")
saveRDS(plco, file = "plco.rds")
