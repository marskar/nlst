library(readr)
library(parsnip)
library(h2o) 
library(here) 
h2o.init()

rds <- "data/nlst_abn_lcp_pre.rds"
data <- read_rds(here(rds))

csv <- "data/nlst_abn_lcp_pre.csv"
h2o_data <- h2o.importFile(here(csv))
dim(h2o_data)
names(h2o_data)

h2o_data$case <- as.factor(h2o_data$case)  #encode the binary response as a factor
h2o.levels(h2o_data$case)  # show the factor levels

splits <- h2o.splitFrame(data = h2o_data, 
                         ratios = c(0.7),  # partition data into 70%, 15% chunks
                         destination_frames = c("train", "test"), # frame ID (not required)
                         seed = 1)  # setting a seed will guarantee reproducibility
train <- splits[[1]]
test <- splits[[2]]

y <- "case"
# Remove the columns that are irrelevant or correlated with the outcome
# x <- setdiff(names(h2o_data), c(
#     y,
#     "candx_days",
#     "cancyr",
#     "canc_rpt_link",
#     "incidence.years",
#     "years.followed",
#     "lung.cancer.death",
#     "pid",
#     "screen_group",
#     "rndgroup",
#     "truefalse_scrnres_ly0",
#     "truefalse_scrnres_ly1",
#     "truefalse_scrnres_ly2",
#     "center",
#     "lss",
#     "scr_res0",
#     "scr_res1",
#     "scr_res2",
#     "scr_iso0",
#     "scr_iso1",
#     "scr_iso2",
#     "other.cause.death"
# ))
x <- c(
    "log1yrisk",
    "emphysema",
    "adenopathy",
    "consolidation"
)

print(x)
length(x)
dim(h2o_data)

# prerisks <- rep("log1yrisk", length(x)-1)
# interact_pairs <- mapply(c, prerisks, x[-length(x)], SIMPLIFY=FALSE)
# TODO
# 0. LCRAT only
glm_lcrat_only <-
    glm(case ~ log1yrisk - 1,
        data = data,
        family = binomial(link = 'log'),
        na.action = na.exclude)
summary(glm_lcrat_only)
AIC(glm_lcrat_only)
pROC::auc(data$case, predict(glm_lcrat_only, type = "response"))

h2o_glm_lcrat_only <- h2o.glm(
    x = "log1yrisk",
    y = y,
    training_frame = train,
    family = "binomial",
    lambda_search = TRUE
)


# 1. LCRAT + CT
glm_fit_lcrat_emph <-
    glm(
        case
        ~ log1yrisk
        + log1yrisk:emphysema
        - 1,
        data = data,
        family = binomial(link = 'log'),
        na.action = na.exclude
    )
summary(glm_lcrat_emph)
AIC(glm_lcrat_emph)
pROC::auc(data$case, predict(glm_lcrat_emph, type = "response"))

glm_lcrat_emph_cons <-
    glm(
        case
        ~ log1yrisk
        + log1yrisk:emphysema
        + log1yrisk:consolidation
        - 1,
        data = data,
        family = binomial(link = 'log'),
        na.action = na.exclude
    )
summary(glm_lcrat_emph_cons)
AIC(glm_lcrat_emph_cons)
pROC::auc(data$case, predict(glm_lcrat_emph_cons, type = "response"))

glm_lcrat_emph_cons_aden <-
    glm(
        case
        ~ log1yrisk
        + log1yrisk:emphysema
        + log1yrisk:consolidation
        + log1yrisk:adenopathy
        - 1,
        data = data,
        family = binomial(link = 'log'),
        na.action = na.exclude
    )
summary(glm_lcrat_emph_cons_aden)
AIC(glm_lcrat_emph_cons_aden)
pROC::auc(data$case, predict(glm_lcrat_emph_cons_aden, type = "response"))

# negative only
# only interactions
# check the n (33k negative screens, same dataset?)
# use log 1yr risk instead of prescreen risk
# use log link in glm function
glm_fit_lcrat_lcp <-
    glm(
        case
        ~ log1yrisk:max_lcp_score
        - 1,
        data = data,
        family = binomial(link = 'log'),
        na.action = na.exclude
    )

glm_fit_lcrat_lcp_emph <-
    glm(
        case
        ~ log1yrisk:max_lcp_score
        + log1yrisk:emphysema
        - 1,
        data = data,
        family = binomial(link = 'log'),
        na.action = na.exclude
    )

glm_fit_lcrat_lcp_emph_pemph <-
    glm(
        case
        ~ log1yrisk:max_lcp_score
        + log1yrisk:emphysema
        + log1yrisk:I(log(p_emph))
        - 1,
        data = df,
        family = binomial(link = 'log'),
        na.action = na.exclude
    )


h2o_glm_lcrat_ct <- h2o.glm(
    x = c(
        "log1yrisk",
        "any_growth",
        "emphysema",
        "consolidation",
        "adenopathy"
        ),
    y = y,
    interaction_pairs = list(
        c("log1yrisk", "any_growth"),
        c("log1yrisk", "emphysema"),
        c("log1yrisk", "consolidation"),
        c("log1yrisk", "adenopathy")
    ),
    training_frame = train,
    family = "binomial",
    lambda_search = TRUE
)
h2o.varimp_plot(h2o_glm_lcrat_ct)

x = c(
    "log1yrisk",
    "emphysema",
    "consolidation",
    "adenopathy"
    )

h2o_glm_lcrat_ct_lcp <- h2o.glm(
    x = x,
    y = y,
    interaction_pairs = list(
        c("log1yrisk", "any_growth"),
        c("log1yrisk", "emphysema"),
        c("log1yrisk", "consolidation"),
        c("log1yrisk", "adenopathy"),
        c("log1yrisk", "max_lcp_score")
    ),
    training_frame = train,
    family = "binomial",
    lambda_search = TRUE
)
h2o.varimp_plot(h2o_glm_lcrat_ct_lcp)


x = c(
    "log1yrisk",
    "emphysema",
    "consolidation",
    "adenopathy",
    "max_lcp_score"
    )

h2o_glm_lcrat_ct_all <- h2o.glm(
    x = x,
    y = y,
    interactions = x,
    training_frame = train,
    family = "binomial",
    lambda_search = TRUE
)
h2o.varimp_plot(h2o_glm_lcrat_ct_all)

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

# TODO 
# for nested models, use lrt (get p-value)
# for non-nested models, look for decrease in deviance, AIC, BIC 
# First step: LCDRAT + lcp_score (composite of nodule features) and emphysema (case ~ prerisk + prerisk:lcp + prerisk:emph + prerisk:adenopathy + prerisk:consolidation) Negative screens only with lcp score and for positive screens nodule features (any spiculation, longest diameter)
# Second step: build top "only interaction terms" model
# Second step: build top "no interaction terms" model
# Second step: build "a mix interaction terms"
# both main effect and interaction
