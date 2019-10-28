#### Load packages ####
library(dplyr)
library(haven)
library(here)
library(purrr)
library(readr)

# Define function for reshaping lungrads data
lrads_extract <- function(df, x) {
    df %>%
        select(c(pid, LRcat = paste0('slungrad', x))) %>%
        mutate("interval" = rep(x, nrow(.)))
}

# Load lung-RADS dataset https://biometry.nci.nih.gov/plcosub/lung/2009-00516-201204-0017-lung-mortality-risk/nlst/lungrads1.sas7bdat/@@download/file/lungrads1.sas7bdat from https://biometry.nci.nih.gov/plcosub/lung/2009-00516-201204-0017-lung-mortality-risk/nlst/lungrads1.sas7bdat/view
# Use above function to reshape lungrads data
# Save RDS file
read_sas(here("data/lungrads1.sas7bdat")) %>%
    map_dfr(0:2, lrads_extract, df = .) %>%
    write_rds(here('data/lungrads.rds'))
