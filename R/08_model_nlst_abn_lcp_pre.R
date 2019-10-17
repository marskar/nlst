# Read in the Kovalchik prediction function
library(survival)
library(glmnet)
library(dplyr)
library(stats)
library(here)
library(readr)
library(tidyr)
source(here("R/kovalchik.R"))
plco <- readRDS("data/plco.rds")
lcp = read_rds("data/nlst_abn_lcp.rds")

# Control arm of PLCO who had no chest xray
plco_control <- subset(plco, control.group == 1)

LCRAT <- coxph(
    Surv(incidence.years, case)
    ~ female
    + race
    + edu6
    + fam.lung.trend
    + emp
    + I(bmi <= 18.5)
    + I(cpd > 20)
    + as.factor(pkyears.cat)
    + I(log(age))
    + I(log(bmi))
    + I(log(qtyears + 1))
    + smkyears,
    data = plco_control
)

cox_death <- coxph(
    Surv(years.followed, other.cause.death)
    ~ female
    + race
    + edu6
    + emp
    + I(bmi <= 18.5)
    + I(cpd > 20)
    + as.factor(pkyears.cat)
    + I((age) ^ 2)
    + I((bmi - 25) ^ 2)
    + I(log(qtyears + 1))
    + smkyears,
    data = plco_control
)

names(lcp)

lcp$prescr_1yrisk <- risk.kovalchik(0, 1, lcp, LCRAT, cox_death)
lcp %>% write_rds("data/nlst_abn_lcp_pre.rds")
lcp %>% write_csv("data/nlst_abn_lcp_pre.csv")

library("h2o") 
h2o.init()

data_file <- "data/nlst_abn_lcp_pre.csv"
data <- h2o.importFile(data_file)
dim(data)
names(data)

data$case <- as.factor(data$case)  #encode the binary repsonse as a factor
h2o.levels(data$case)  # show the factor levels

splits <- h2o.splitFrame(data = data, 
                         ratios = c(0.7),  # partition data into 70%, 15% chunks
                         destination_frames = c("train", "test"), # frame ID (not required)
                         seed = 1)  # setting a seed will guarantee reproducibility
train <- splits[[1]]
test <- splits[[2]]

y <- "case"
# Remove the columns that are irrelevant or correlated with the outcome
x <- setdiff(names(data), c(
    y,
    "candx_days",
    "cancyr",
    "canc_rpt_link",
    "incidence.years",
    "years.followed",
    "lung.cancer.death",
    "pid",
    "screen_group",
    "rndgroup",
    "truefalse_scrnres_ly0",
    "truefalse_scrnres_ly1",
    "truefalse_scrnres_ly2",
    "center",
    "lss",
    "scr_res0",
    "scr_res1",
    "scr_res2",
    "scr_iso0",
    "scr_iso1",
    "scr_iso2",
    "other.cause.death"
))
print(x)
length(x)
dim(data)

prerisks <- rep("prescr_1yrisk", length(x)-1)
interact_pairs <- mapply(c, prerisks, x[-length(x)], SIMPLIFY=FALSE)
# TODO
# 0. LCRAT only
glm_fit_lcrat_only <- h2o.glm(
    x = "prescr_1yrisk",
    y = y,
    training_frame = train,
    family = "binomial",
    lambda_search = TRUE
)

summary(glm(case ~ log1yrisk:consolidation + log1yrisk:adenopathy -1, data=data, family=binomial(link='log')))

# 1. LCRAT + CT
glm_fit1 <- h2o.glm(
    x = c(
        "prescr_1yrisk",
        "any_growth",
        "emphysema",
        "consolidation",
        "adenopathy"
        ),
    y = y,
    interaction_pairs = list(
        c("prescr_1yrisk", "any_growth"),
        c("prescr_1yrisk", "emphysema"),
        c("prescr_1yrisk", "consolidation"),
        c("prescr_1yrisk", "adenopathy")
    ),
    training_frame = train,
    family = "binomial",
    lambda_search = TRUE
)

glm_fit1 <- h2o.glm(
    x = c(
        "prescr_1yrisk",
        "any_growth",
        "emphysema",
        "consolidation",
        "adenopathy",
        "max_lcp_score"
        ),
    y = y,
    interaction_pairs = list(
        c("prescr_1yrisk", "any_growth"),
        c("prescr_1yrisk", "emphysema"),
        c("prescr_1yrisk", "consolidation"),
        c("prescr_1yrisk", "adenopathy"),
        c("prescr_1yrisk", "max_lcp_score")
    ),
    training_frame = train,
    family = "binomial",
    lambda_search = TRUE
)

glm_fit_all_pairs <- h2o.glm(
    x = x,
    y = y,
    interaction_pairs = interact_pairs,
    training_frame = train,
    family = "binomial",
    lambda_search = TRUE
)

glm_fit_all <- h2o.glm(
    x = x,
    y = y,
    interactions = x,
    training_frame = train,
    family = "binomial",
    lambda_search = TRUE
)
# 2. LCRAT + CT + LCP_SCORE
# 3. LCRAT + CT + EMPH_SCORE
# 4. LCRAT + CT + LCP_SCORE + EMPH_SCORE
# Describe these models
# Methods: LCRAT + CT (Hilary) + LCP_SCORE (Optellum)
# First, prescreen_risk interactions with all binary features
# Second, prescreen_risk interactions with all discrete (categorical) features
glm_fit1 <- h2o.glm(
    x = x,
    y = y,
    training_frame = train,
    family = "binomial",
    lambda_search = TRUE
)  #Like glm() and glmnet(), h2o.glm() has the family argument

glm_fit1@model$model_summary
coefs = h2o.coef(glm_fit_lcrat_only)
coefs = h2o.coef(glm_fit1)
h2o.coef(glm_fit_all)
coefs[abs(coefs) > 0]
h2o.varimp_plot(glm_fit1)

