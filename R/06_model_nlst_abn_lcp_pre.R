library(readr)
library(broom)
library(dplyr)
library(tidyr)
library(ggplot2)
library(here) 
library(mlr)
library(glmnet)
library(coefplot)
library(leaps)
library(tidyposterior)
library(rsample)

data <-
    read_rds(here("data/nlst_abn_lcp_pre.rds")) %>%
    select(
        case,
        logit1yrisk,
        diam_cat,
        any_growth,
        emphysema,
        consolidation,
        adenopathy,
        any_upper,
        any_right_mid,
        any_lingula,
        any_mixed,
        any_spiculation,
        any_poor_def,
        max_lcp_score,
        any_margin_unab
           ) %>% 
    mutate(diam_cat = as.factor(diam_cat)) %>% 
    drop_na()
    # replace_na(list(0))

# Best subset selection using leaps package
reg <-
    leaps::regsubsets(
        x = case ~ logit1yrisk * .,
        data = data,
        nvmax = 50,
        method = "exhaustive"
    )
regsum <- summary(reg)

## Best model according to BIC
qplot(seq(length(regsum$bic)), regsum$bic)
best_bic <- which(regsum$bic == min(regsum$bic))
names(regsum$which[best_bic, regsum$which[best_bic,]])
## Best model according to Cp
qplot(seq(length(regsum$cp)), regsum$cp)
best_cp <- which(regsum$cp == min(regsum$cp))
names(regsum$which[best_cp, regsum$which[best_cp,]])

# Other metrics
qplot(seq(length(regsum$rss)), regsum$rss)
qplot(seq(length(regsum$rsq)), regsum$rsq)
qplot(seq(length(regsum$adjr2)), regsum$adjr2)

# Best model using Lasso
x = model.matrix(case ~ logit1yrisk + logit1yrisk*. -1, data = data)
cv.lasso = cv.glmnet(x, data$case, alpha = 1, family = "binomial")
qplot(x = cv.lasso$lambda, y = cv.lasso$nzero)
plot(cv.lasso)
coefpath(cv.lasso)
coefplot(cv.lasso)

## The lambda resulting in minimal error
## does not reduce any coefficients to zero
coef(cv.lasso, s = "lambda.min")
## most regularized model which has an error
## within one standard error of the minimum
coef(cv.lasso, s = "lambda.1se")

# Tidy models
folds <- vfold_cv(data, v = 10)
folds$splits


# MLR - Wrapper
lrnr = makeLearner("classif.logreg")
task = makeClassifTask(id = "nlst", data, "case")
ctrl = makeFeatSelControlExhaustive()
rdesc = makeResampleDesc("CV", iters = 10)
sfeats = selectFeatures(lrnr, task, resampling = rdesc, control = ctrl, show.info = FALSE)
sfeats$x
sfeats$y
analyzeFeatSelResult(sfeats)

# MLR - Filter
fv = generateFilterValuesData(task, method = c("FSelectorRcpp_information.gain"))

plotFilterValues(fv) + coord_flip()

plotpval <- function(model) {
    print(paste("AIC =", AIC(model)))
    ggplot(tidy(model), aes(reorder(term, -log10(p.value)), -log10(p.value))) +
        geom_col(aes(fill = term)) +
        coord_flip() +
        theme(legend.position = "none") +
        ylab("Negative log10 p-value") +
        xlab("")
}



# 0. LCRAT only
# glm_lcrat_only_log <-
#     glm(case ~ log1yrisk - 1,
#         data = data,
#         family = binomial(link = 'log'),
#         na.action = na.exclude)
# summary(glm_lcrat_only_log)
# 
# 
# # 1. LCRAT + CT interacted log ----
# glm_interacted_log <-
#     glm(
#         case
#         ~ log1yrisk
#         + log1yrisk:diam_cat
#         + log1yrisk:any_growth
#         + log1yrisk:emphysema
#         + log1yrisk:consolidation
#         + log1yrisk:adenopathy
#         + log1yrisk:any_upper
#         + log1yrisk:any_right_mid
#         + log1yrisk:any_lingula
#         + log1yrisk:any_mixed
#         + log1yrisk:any_spiculation
#         + log1yrisk:any_poor_def
#         + log1yrisk:any_margin_unab
#         # + log1yrisk:max_lcp_score
#         - 1,
#         data = data,
#         family = binomial(link = 'log'),
#         na.action = na.exclude,
#         start = rep(0.001, 13)
#     )
# summary(glm_interacted_log)
# 
# 
# glm_interacted_lcp_log  <-
#     glm(
#         case
#         ~ log1yrisk
#         + log1yrisk:diam_cat
#         + log1yrisk:any_growth
#         + log1yrisk:emphysema
#         + log1yrisk:consolidation
#         + log1yrisk:adenopathy
#         + log1yrisk:any_upper
#         + log1yrisk:any_right_mid
#         + log1yrisk:any_lingula
#         + log1yrisk:any_mixed
#         + log1yrisk:any_spiculation
#         + log1yrisk:any_poor_def
#         + log1yrisk:any_margin_unab
#         + log1yrisk:max_lcp_score
#         - 1,
#         data = data,
#         family = binomial(link = 'log'),
#         na.action = na.exclude,
#         start = rep(0.001, 14)
#     )
# summary(glm_interacted_lcp_log)
# plotpval(glm_interacted_lcp_log)
# lmtest::lrtest(glm_interacted_log, glm_interacted_lcp_log)
# 
# 
# # 2. LCRAT + CT uninteracted log ----
# glm_uninteracted_log <-
#     glm(
#         case
#         ~ log1yrisk
#         + diam_cat
#         + any_growth
#         + emphysema
#         + consolidation
#         + adenopathy
#         + any_upper
#         + any_right_mid
#         + any_lingula
#         + any_mixed
#         + any_spiculation
#         + any_poor_def
#         + any_margin_unab
#         # + log1yrisk:max_lcp_score
#         - 1,
#         data = data,
#         family = binomial(link = 'log'),
#         na.action = na.exclude,
#         start = rep(0.001, 13)
#     )
# summary(glm_uninteracted_log)
# 
# glm_uninteracted_lcp_log  <-
#     glm(
#         case
#         ~ log1yrisk
#         + diam_cat
#         + any_growth
#         + emphysema
#         + consolidation
#         + adenopathy
#         + any_upper
#         + any_right_mid
#         + any_lingula
#         + any_mixed
#         + any_spiculation
#         + any_poor_def
#         + any_margin_unab
#         + log1yrisk:max_lcp_score
#         - 1,
#         data = data,
#         family = binomial(link = 'log'),
#         na.action = na.exclude,
#         start = rep(0.001, 14)
#     )
# summary(glm_uninteracted_lcp_log)
# plotpval(glm_uninteracted_lcp_log)
# lmtest::lrtest(glm_uninteracted_log, glm_uninteracted_lcp_log)
# 
# 
# # 3. LCRAT + CT hybrid log ----
# glm_hybrid_log <-
#     glm(
#         case
#         ~ log1yrisk
#         + log1yrisk*diam_cat
#         + log1yrisk*any_growth
#         + log1yrisk*emphysema
#         + log1yrisk*consolidation
#         + log1yrisk*adenopathy
#         + log1yrisk*any_upper
#         + log1yrisk*any_right_mid
#         + log1yrisk*any_lingula
#         + log1yrisk*any_mixed
#         + log1yrisk*any_spiculation
#         + log1yrisk*any_poor_def
#         + log1yrisk*any_margin_unab
#         # + log1yrisk*max_lcp_score
#         - 1,
#         data = data,
#         family = binomial(link = 'log'),
#         na.action = na.exclude,
#         start = rep(0.001, 25)
#     )
# summary(glm_hybrid_log)
# 
# glm_hybrid_lcp_log <-
#     glm(
#         case
#         ~ log1yrisk
#         + log1yrisk*diam_cat
#         + log1yrisk*any_growth
#         + log1yrisk*emphysema
#         + log1yrisk*consolidation
#         + log1yrisk*adenopathy
#         + log1yrisk*any_upper
#         + log1yrisk*any_right_mid
#         + log1yrisk*any_lingula
#         + log1yrisk*any_mixed
#         + log1yrisk*any_spiculation
#         + log1yrisk*any_poor_def
#         + log1yrisk*any_margin_unab
#         + log1yrisk*max_lcp_score
#         - 1,
#         data = data,
#         family = binomial(link = 'log'),
#         na.action = na.exclude,
#         start = rep(0.001, 27)
#     )
# summary(glm_hybrid_lcp_log)
# plotpval(glm_hybrid_lcp_log)
# lmtest::lrtest(glm_hybrid_log, glm_hybrid_lcp_log)


# 0. LCRAT only logit
glm_lcrat_only_logit <-
    glm(case ~ logit1yrisk - 1,
        data = data,
        family = binomial(link = 'logit'),
        na.action = na.exclude)
summary(glm_lcrat_only_logit)


# 1. LCRAT + CT interacted logit ----
glm_interacted_logit <-
    glm(
        case
        ~ logit1yrisk
        + logit1yrisk*diam_cat
        + logit1yrisk*any_growth
        + logit1yrisk*emphysema
        + logit1yrisk*consolidation
        + logit1yrisk*adenopathy
        + logit1yrisk*any_upper
        + logit1yrisk*any_right_mid
        + logit1yrisk*any_lingula
        + logit1yrisk*any_mixed
        + logit1yrisk*any_spiculation
        + logit1yrisk*any_poor_def
        + logit1yrisk*any_margin_unab
        # + logit1yrisk:max_lcp_score
        - 1,
        data = data,
        family = binomial(link = 'logit'),
        na.action = na.exclude,
    )
summary(glm_interacted_logit)
levels(data$diam_cat)



glm_interacted_lcp_logit <-
    glm(
        case
        ~ logit1yrisk
        + logit1yrisk:diam_cat
        + logit1yrisk:any_growth
        + logit1yrisk:emphysema
        + logit1yrisk:consolidation
        + logit1yrisk:adenopathy
        + logit1yrisk:any_upper
        + logit1yrisk:any_right_mid
        + logit1yrisk:any_lingula
        + logit1yrisk:any_mixed
        + logit1yrisk:any_spiculation
        + logit1yrisk:any_poor_def
        + logit1yrisk:any_margin_unab
        + max_lcp_score
        + max_lcp_score:diam_cat
        + max_lcp_score:any_growth
        + max_lcp_score:emphysema
        + max_lcp_score:consolidation
        + max_lcp_score:adenopathy
        + max_lcp_score:any_upper
        + max_lcp_score:any_right_mid
        + max_lcp_score:any_lingula
        + max_lcp_score:any_mixed
        + max_lcp_score:any_spiculation
        + max_lcp_score:any_poor_def
        + max_lcp_score:any_margin_unab
        - 1,
        data = data,
        family = binomial(link = 'logit'),
        na.action = na.exclude,
    )
summary(glm_interacted_lcp_logit)
plotpval(glm_interacted_lcp_logit)
lmtest::lrtest(glm_interacted_logit, glm_interacted_lcp_logit)


# 2. LCRAT + CT uninteracted logit ----
glm_uninteracted_logit <-
    glm(
        case
        ~ logit1yrisk
        + diam_cat
        + any_growth
        + emphysema
        + consolidation
        + adenopathy
        + any_upper
        + any_right_mid
        + any_lingula
        + any_mixed
        + any_spiculation
        + any_poor_def
        + any_margin_unab
        # + log1yrisk:max_lcp_score
        + 1,
        data = data,
        family = binomial(link = 'logit'),
        na.action = na.exclude,
    )
summary(glm_uninteracted_logit)

glm_uninteracted_lcp_logit <-
    glm(
        case
        ~ logit1yrisk
        + diam_cat
        + any_growth
        + emphysema
        + consolidation
        + adenopathy
        + any_upper
        + any_right_mid
        + any_lingula
        + any_mixed
        + any_spiculation
        + any_poor_def
        + any_margin_unab
        + max_lcp_score
        - 1,
        data = data,
        family = binomial(link = 'logit'),
        na.action = na.exclude,
    )
summary(glm_uninteracted_lcp_logit)
plotpval(glm_uninteracted_lcp_logit)
lmtest::lrtest(glm_uninteracted_logit, glm_uninteracted_lcp_logit)

# TODO 
# Find out what are the variables in abnlist.neg.int in  the negatives script
# Define final without max_lcp_score (AIC and LRT)
# And then redo model selection with the score (AIC and LRT)
# Set up new call with Optellum

# 3. LCRAT + CT hybrid logit ----
glm_hybrid_logit <-
    glm(
        case
        ~ logit1yrisk
        + logit1yrisk*diam_cat
        + logit1yrisk*any_growth
        + logit1yrisk*emphysema
        + logit1yrisk*consolidation
        + logit1yrisk*adenopathy
        + logit1yrisk*any_upper
        + logit1yrisk*any_right_mid
        + logit1yrisk*any_lingula
        + logit1yrisk*any_mixed
        + logit1yrisk*any_spiculation
        + logit1yrisk*any_poor_def
        + logit1yrisk*any_margin_unab
        # + log1yrisk*max_lcp_score
        - 1,
        data = data,
        family = binomial(link = 'logit'),
        na.action = na.exclude,
    )
summary(glm_hybrid_logit)

glm_hybrid_lcp_logit <-
    glm(
        case
        ~ logit1yrisk
        + logit1yrisk*diam_cat
        + logit1yrisk*any_growth
        + logit1yrisk*emphysema
        + logit1yrisk*consolidation
        + logit1yrisk*adenopathy
        + logit1yrisk*any_upper
        + logit1yrisk*any_right_mid
        + logit1yrisk*any_lingula
        + logit1yrisk*any_mixed
        + logit1yrisk*any_spiculation
        + logit1yrisk*any_poor_def
        + logit1yrisk*any_margin_unab
        + logit1yrisk*max_lcp_score
        - 1,
        data = data,
        family = binomial(link = 'logit'),
        na.action = na.exclude,
    )
summary(glm_hybrid_lcp_logit)
plotpval(glm_hybrid_lcp_logit)
lmtest::lrtest(glm_hybrid_logit, glm_hybrid_lcp_logit)


# glm_lcrat_emph_cons <-
#     glm(
#         case
#         ~ log1yrisk
#         + log1yrisk:emphysema
#         + log1yrisk:consolidation
#         - 1,
#         data = data,
#         family = binomial(link = 'log'),
#         na.action = na.exclude
#     )
# summary(glm_lcrat_emph_cons)
# plotpval(glm_lcrat_emph_cons)
#        
# glm_lcrat_emph_cons_aden <-
#     glm(
#         case
#         ~ log1yrisk
#         + log1yrisk:emphysema
#         + log1yrisk:consolidation
#         + log1yrisk:adenopathy
#         - 1,
#         data = data,
#         family = binomial(link = 'log'),
#         na.action = na.exclude,
#         start = rep(0.001, 13)
#     )
# summary(glm_lcrat_emph_cons_aden)
# plotpval(glm_lcrat_emph_cons_aden)
# 
# # negative only
# # only interactions
# # check the n (33k negative screens, same dataset?)
# # use log 1yr risk instead of prescreen risk
# # use log link in glm function
# glm_lcrat_pemph <-
#     glm(
#         case
#         ~ log1yrisk
#         # + log1yrisk:emphysema
#         + log1yrisk:I(log(p_emph))
#         - 1,
#         data = data,
#         family = binomial(link = 'log'),
#         na.action = na.exclude
#     )
# summary(glm_lcrat_pemph)
# plotpval(glm_lcrat_pemph)
# 
# glm_lcrat_emph_pemph <-
#     glm(
#         case
#         ~ log1yrisk
#         + log1yrisk:emphysema
#         + log1yrisk:p_emph
#         - 1,
#         data = data,
#         family = binomial(link = 'log'),
#         na.action = na.exclude
#     )
# plotpval(glm_lcrat_emph_pemph)
# 
# glm_lcrat_lcp <-
#     glm(
#         case
#         ~ log1yrisk
#         + log1yrisk:max_lcp_score
#         - 1,
#         data = data,
#         family = binomial(link = 'log'),
#         na.action = na.exclude
#     )
# plotpval(glm_lcrat_lcp)
# 
# glm_lcrat_lcp <-
#     glm(
#         case
#         ~ log1yrisk:max_lcp_score
#         - 1,
#         data = data,
#         family = binomial(link = 'log'),
#         na.action = na.exclude
#     )
# plotpval(glm_lcrat_lcp)
# names(data)
# glm_lcrat_lcp_emph <-
#     glm(
#         case
#         ~ logit1yrisk
#         + logit1yrisk:max_lcp_score
#         + logit1yrisk:emphysema
#         + logit1yrisk:consolidation
#         + logit1yrisk:adenopathy
#         - 1,
#         data = data,
#         family = binomial(link = 'logit'),
#         na.action = na.exclude
#     )
# summary(glm_lcrat_lcp_emph)
# plotpval(glm_lcrat_lcp_emph)
# names(data)
# 
#         
# glm_logit <-
#     glm(
#         case
#         ~ logit1yrisk:diam_cat
#         + logit1yrisk:any_growth
#         + logit1yrisk:emphysema
#         + logit1yrisk:consolidation
#         + logit1yrisk:adenopathy
#         + logit1yrisk:any_upper
#         + logit1yrisk:I(any_right_mid == 1 | any_lingula == 1)
#         + logit1yrisk:any_mixed
#         + logit1yrisk:any_spiculation
#         + logit1yrisk:I(any_poor_def == 1 | any.margin.unab == 1)
#         # + logit1yrisk:max_lcp_score
#         - 1,
#         data = data,
#         family = binomial(link = 'logit')
#     )
# 
# summary(glm_logit)
# tidy(glm_logit)
# tail(data$max_lcp_score)
# str(data)
# glm_fit_lcrat_lcp_emph_pemph <-
#     glm(
#         case
#         ~ log1yrisk
#         + log1yrisk:max_lcp_score
#         + log1yrisk:emphysema
#         + log1yrisk:I(log(p_emph))
#         - 1,
#         data = data,
#         family = binomial(link = 'log'),
#         na.action = na.exclude
#     )
# # plotpval(glm_lcrat_lcp_emph)
# 
# x = c(
#     "log1yrisk",
#     "emphysema",
#     "consolidation",
#     "adenopathy"
#     )


# TODO 
# for nested models, use lrt (get p-value)
# for non-nested models, look for decrease in deviance, AIC, BIC 
# First step: LCDRAT + lcp_score (composite of nodule features) and emphysema (case ~ prerisk + prerisk:lcp + prerisk:emph + prerisk:adenopathy + prerisk:consolidation) Negative screens only with lcp score and for positive screens nodule features (any spiculation, longest diameter)
# Second step: build top "only interaction terms" model
# Second step: build top "no interaction terms" model
# Second step: build "a mix interaction terms"
# both main effect and interaction
