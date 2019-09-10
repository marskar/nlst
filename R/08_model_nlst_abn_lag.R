# Read in the Kovalchik prediction function
library(survival)
library(glmnet)
library(dplyr)
library(here)
library(tidyr)
source(here("R/kovalchik.R"))
plco <- readRDS("data/plco.rds")
lag = readRDS(file = "data/nlst_abn_lag.rds")

plco_control <-
    subset(plco, control.group == 1) # control arm of PLCO who had no chest xray

LCRAT <- coxph(
    Surv(incidence.years, case) ~
        female + race + edu6 + fam.lung.trend + emp + I(bmi <=
                                                            18.5) + I(cpd > 20) + as.factor(pkyears.cat) +
        I(log(age)) + I(log(bmi)) + I(log(qtyears + 1)) + smkyears,
    data = plco_control
)
cox_death <- coxph(
    Surv(years.followed, other.cause.death) ~
        female + race + edu6 + emp + I(bmi <= 18.5) + I(cpd >
                                                            20) + as.factor(pkyears.cat) + I((age) ^ 2) + I((bmi - 25) ^ 2) +
        I(log(qtyears + 1)) + smkyears,
    data = plco_control
)

colname_vector <- c(names(lag)[36:162])

lag$prescr.1yrisk <- risk.kovalchik(0, 1, lag, LCRAT, cox_death) 
lag <- lag %>% 
    mutate(
        log1yrisk = log(prescr.1yrisk),
        logit1yrisk = log(prescr.1yrisk / (1 - prescr.1yrisk))
    ) %>% 
    drop_na(colname_vector) %>% 
    filter_all(all_vars(is.finite(.)))

nrow(lag)

int_full <- glm(paste("case ~ logit1yrisk -1 +",paste0(colname_vector, collapse="+")), data=lag, family=binomial(link='logit'), na.action = "na.exclude")
bsw.int <- step(int_full, direction="backward", scope = list(lower = case ~ logit1yrisk -1, upper = int_full))

  # Look at a model including these 4 effects
summary(glm(case ~ log1yrisk + log1yrisk:other.above + log1yrisk:benign.nodule + log1yrisk:consolidation + log1yrisk:adenopathy -1, data=data_interval_abn, family=binomial(link='log')))

x = model.matrix(case ~ logit1yrisk + logit1yrisk:. - 1, data = lag[, c("case", "logit1yrisk", colname_vector)])
nrow(x)
nrow(lag)
cv.lasso = cv.glmnet(x,
                     lag,
                     alpha = 1,
                     family = "binomial")
out <-
    glmnet(x,
           data.screen.abn.pos$case,
           alpha = 1,
           family = "binomial")
predict(out,
        type = "coefficients",
        s = (cv.lasso$lambda.min + cv.lasso$lambda.1se) / 2)
# Fit this model - log scale
mod.vars.lasso <-
    glm(
        case ~ log1yrisk:any.upper + log1yrisk:any.spiculation + log1yrisk:any.smooth + log1yrisk:any.mixed +
            log1yrisk:any.other.att + log1yrisk:any.growth + log1yrisk:diam.cat -
            1,
        data = data.screen.abn.pos,
        family = binomial(link = 'log')
    )
summary(mod.vars.lasso)


View(lag)
