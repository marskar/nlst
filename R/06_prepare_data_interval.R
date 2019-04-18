#### Setup ####

# List packages to be loaded (and installed if needed)
packages <-
    c(
        "lmtest",
        "here",
        "dplyr",
        "ggplot2",
        "survival",
        "gmodels",
        "devtools",
        "geepack",
        "MESS",
        "psych",
        "Hmisc",
        "glmnet",
        "boot"
    )

# List packages that are not installed
not_installed <-
    packages[!(packages %in% installed.packages()[, "Package"])]

# Install packages that are not installed
if (length(not_installed))
    install.packages(not_installed)

# Load all packages
lapply(packages, require, character.only = TRUE)

# Read in the Kovalchik prediction function
source(here("R/kovalchik.R"))

devtools::install_github('marskar/coxph_risk')
devtools::install_github('marskar/lcmodels')

# Function to create "not in" operator
'%!in%' <- function(x, y)
    ! ('%in%'(x, y))


# Load NLST data
nlst_ct_t1_neg = readRDS('data/nlst_ct_t1_neg.rds')
nlst_ct_t2_neg = readRDS('data/nlst_ct_t2_neg.rds')
nlst_ct_t3_neg = readRDS('data/nlst_ct_t3_neg.rds')

### Construct dataset to model risk of ALL interval cancers across all 3 screens
  # Interval variable in data.interval datasets is 1 for T0-T1, 2 for T1-T2, and 3 for post-T2 intervals
data.frame(pid=c(nlst_ct_t1_neg$pid, nlst_ct_t2_neg$pid, nlst_ct_t3_neg$pid),
                            case=c(nlst_ct_t1_neg$case_T1_interval, nlst_ct_t2_neg$case_T2_interval, nlst_ct_t3_neg$case_T3_interval),
                            interval=c(rep(1,times=nrow(nlst_ct_t1_neg)), rep(2, times=nrow(nlst_ct_t2_neg)), rep(3, times=nrow(nlst_ct_t3_neg)))) %>% 
# Merge this back with covariates from NLST. Add screening history variable.
merge(nlst.sub, by="pid", all.x=T)
merge(select(nlst.CT, pid, hist.T0.T1, hist.T1.T2), by="pid", all.x=T, all.y=F) %>% 
data.interval$screen.hist <- ifelse(data.interval$interval==2, data.interval$hist.T0.T1, ifelse(data.interval$interval==3, data.interval$hist.T1.T2, NA))
data.interval$screen.hist <- factor(data.interval$screen.hist, levels=c(1,2,3,4), labels=c("Neg-Neg","Neg-Pos","Pos-Neg","Pos-Pos"))  
data.interval <- select(data.interval, -c(hist.T0.T1,hist.T1.T2))                        ## delete??
data.interval <- arrange(data.interval, pid, interval)
# Update age, quit-years, and smoke-years by adding a year for T1 and T2
data.interval <- mutate(data.interval, age=ifelse(interval==2, age+1, ifelse(interval==3, age+2, age)),
                smkyears=ifelse(interval==2 & qtyears==0, smkyears+1, ifelse(interval==3 & qtyears==0, smkyears+2, smkyears)),
                qtyears=ifelse(interval==2 & qtyears>0, qtyears+1, ifelse(interval==3 & qtyears>0, qtyears+2, qtyears)))
data.interval <- mutate(data.interval, pkyears.cont=cpd*smkyears/20)  # using new smoke-years, update pack-years, then re-categorize
data.interval <- mutate(data.interval, pkyears.cat=as.factor(ifelse(pkyears.cont>=30 & pkyears.cont<40, "[30,40)",
                  ifelse(pkyears.cont>=40 & pkyears.cont<50, "[40,50)", ifelse(pkyears.cont>=50 & pkyears.cont<999,"[50,Inf]",NA)))))
