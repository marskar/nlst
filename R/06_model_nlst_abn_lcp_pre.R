library(readr)
library(here)
library(dplyr)
library(rsample)
library(recipes)
library(parsnip)
# remotes::install_github("tidymodels/tune")
library(tune)
library(yardstick)
library(ggplot2)
data("two_class_example")
two_class_example
data <-
    read_rds(here("data/nlst_abn_lcp_pre.rds")) %>%
    select(
        # pid,
        # diam_cat,
        case_at_next_screen,
        # screen_result,
        # interval,
        log1yrisk,
        logit1yrisk,
        # longest_diam,
        any_growth,
        emphysema,
        consolidation,
        adenopathy,
        any_upper,
        any_right_mid,
        any_lingula,
        any_mixed,
        # any_spiculation,
        # any_poor_def,
        max_lcp_score,
        # any_margin_unab
           ) %>%
    mutate(case_at_next_screen = as.factor(case_at_next_screen)) %>%
    identity()
    # filter(is.finite(longest_diam))

set.seed(93599150)
logit_mod <- logistic_reg(mode = "classification") %>% set_engine("glm")
ctrl <- control_resamples(save_pred = TRUE)
cv_train <- vfold_cv(data = data, v = 5, strata = "case_at_next_screen")

resampled <- fit_resamples(case_at_next_screen ~ ., logit_mod, resamples = cv_train, control = ctrl)
preds <- collect_predictions(resampled)
preds %>% filter(case_at_next_screen == 1) %>% arrange(desc(.pred_1))
preds %>% arrange(desc(.pred_1))
preds %>% conf_mat(truth=case_at_next_screen, estimate=.pred_class)
preds
df_roc <- preds %>% roc_curve(truth=case_at_next_screen, estimate=.pred_1) %>% arrange(.threshold)
df_pr <- preds %>% pr_curve(truth=case_at_next_screen, estimate=.pred_1) %>% arrange(.threshold)
df_pr %>% tail
df_roc %>% head
df_roc <- df_roc %>% rename(recall = sensitivity)
lorenz <- df_roc %>% filter(is.finite(.threshold)) %>% mutate(proportion_considered_high_risk = .threshold / max(.threshold))
lorenz %>% 
ggplot() +
    aes(x=proportion_considered_high_risk, y=sensitivity) +
    geom_line() +
    geom_abline(intercept=0)

mean(na.omit(preds$.pred_1) < 0.001)
preds$.pred_1 %in% df_roc$.threshold
dim(na.omit(df_roc))
dim(na.omit(df_pr))
dim(na.omit(preds))
sum(preds$.pred_1 < 0.000318)
df_roc %>% na.omit %>% mutate(risk = na.omit(preds$.pred_1))
full_join(df_pr, df_roc, by="recall")
preds %>% roc_curve(truth=case_at_next_screen, estimate=.pred_1) %>% autoplot()
preds %>% pr_curve(truth=case_at_next_screen, estimate=.pred_1) %>% autoplot()

preds %>% lift_curve(truth=case_at_next_screen, estimate=.pred_1) %>% autoplot()
preds %>% gain_curve(truth=case_at_next_screen, estimate=.pred_1) %>% autoplot()

logit_mod <- 
    logistic_reg(mode = "classification") %>%
    set_engine(engine = "glm")

## compute mod on kept part
cv_fit <- function(splits, mod, ...) {
    res_mod <-
        fit(mod, case_at_next_screen ~ ., data = analysis(splits), family = binomial)
    return(res_mod)
}

## get predictions on holdout sets
cv_pred <- function(splits, mod){
    # Save the 10%
    holdout <- assessment(splits)
    pred_assess <- bind_cols(truth = holdout$case_at_next_screen, predict(mod, new_data = holdout))
    return(pred_assess)
}

## get probs on holdout sets
cv_prob <- function(splits, mod){
    holdout <- assessment(splits)
    prob_assess <- bind_cols(truth = as.factor(holdout$case_at_next_screen), 
                             predict(mod, new_data = holdout, type = "prob"))
    return(prob_assess)
}
cv_train %>% mutate(res_mod = map(splits, .f = cv_fit, logit_mod)) ## fit model
res_cv_train <- 
    cv_train %>% 
    mutate(res_mod = map(splits, .f = cv_fit, logit_mod), ## fit model
           res_pred = map2(splits, res_mod, .f = cv_pred), ## predictions
           res_prob = map2(splits, res_mod, .f = cv_prob)) ## probabilities

simple_recipe <- function(dataset) {
    recipe(price ~ ., data = dataset) %>%
        step_center(all_numeric()) %>%
        step_scale(all_numeric()) %>%
        step_dummy(all_nominal())
}

library(stringr)
library(leaps)

library(broom)
library(tidyr)
library(ggplot2)
library(mlr)
library(glmnet)
library(coefplot)
library(tidyposterior)
library(rsample)
library(lme4)

data <-
    read_rds(here("data/nlst_abn_lcp_pre.rds")) %>%
    select(
        pid,
        # diam_cat,
        case_at_next_screen,
        # screen_result,
        # interval,
        log1yrisk,
        logit1yrisk,
        longest_diam,
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
    filter(is.finite(longest_diam))
    # mutate(interval = as.factor(interval))

names(data)
# Model selection decides whether variable is a main effect or interaction

get_best_bic_formula <- function(form, data, yname) {
    regsumm <-
        summary(leaps::regsubsets(
            x = form,
            data = data,
            nvmax = 50,
            method = "exhaustive"
        ))
    coefs <- regsumm$which[which.min(regsumm$bic), -1]
    predictors <- names(which(coefs == TRUE)) %>%
        str_remove("\\d$") %>%
        unique() %>%
        paste(collapse = "+")
    as.formula(paste(yname, "~", predictors))
}
form <- get_best_bic_formula(case_at_next_screen ~ logit1yrisk * . + 1, data, "case_at_next_screen")
form

glm_best_bic <-
    glm(
        form,
        data = data,
        family = binomial(link = 'logit'),
        na.action = na.exclude,
    )
summary(glm_best_bic)

tidy(coef(glm_best_bic))
coefplot(glm_best_bic)

# Test homogeneity
# Interact interval
# Use case_at_next_screen as the dependent variable (y)
glm_interacted_logit <-
    glm(
        case_at_next_screen
        ~ interval*logit1yrisk
        # + interval*longest_diam
        # + interval*any_growth
        # + interval*emphysema
        # + interval*consolidation
        # + interval*adenopathy
        # + interval*any_upper
        # + interval*any_right_mid
        # + interval*any_lingula
        # + interval*any_mixed
        # + interval*any_spiculation
        # + interval*any_poor_def
        # + interval*any_margin_unab
        # + interval*max_lcp_score
        # + longest_diam
        # + any_growth
        # + emphysema
        # + consolidation
        # + adenopathy
        # + any_upper
        # + any_right_mid
        # + any_lingula
        # + any_mixed
        # + any_spiculation
        # + any_poor_def
        # + any_margin_unab
        + interval*max_lcp_score
        - 1,
        data = data,
        family = binomial(link = 'logit'),
        na.action = na.exclude,
    )
summary(glm_interacted_logit)
sum(data$case_at_next_screen)
mixed <- glmer(
    case
    ~ logit1yrisk
    + interval
    + max_lcp_score
    + (1
       + logit1yrisk
       + max_lcp_score
       | interval),
        data = data,
        family = binomial(link = 'logit'),
        na.action = na.exclude
       ) 
summary(mixed)
# Best subset selection using leaps package
reg <-
    leaps::regsubsets(
        x = case ~ logit1yrisk * . + 1,
        data = data,
        nvmax = 50,
        method = "exhaustive"
    )
regsumm <- summary(reg)
bic <- regsumm$bic
cp <- regsumm$cp
bool <- regsumm$which

## Best model according to BIC
qplot(seq(length(bic)), bic)
best_bic <- which.min(bic)
varlist <- names(bool[best_bic, bool[best_bic,]])
tidy(coef(reg, best_bic))
varlist[-1]
best_bic
coefs <- bool[best_bic,-1]
coefs
# Get outcome variable
#form <- as.formula(object$call[[2]])
#outcome <- all.vars(form)[1]
# Get model predictors
predictors <- names(which(coefs == TRUE))
predictors <- paste(predictors, collapse = "+")
form <- as.formula(paste0("case_at_next_screen", "~", predictors))
form
## Best model according to Cp
qplot(seq(length(cp)), cp)
best_cp <- which.min(cp)
best_cp
varlist <- names(bool[best_cp, bool[best_cp,]])
coef(reg, which.min(cp))
varlist[-1]

# Other metrics
qplot(seq(length(regsumm$rss)), regsumm$rss)
qplot(seq(length(regsumm$rsq)), regsumm$rsq)
qplot(seq(length(regsumm$adjr2)), regsumm$adjr2)

# Best model using Lasso
x = model.matrix(case ~ logit1yrisk * . + 1,
                 data = data)
cv.lasso = cv.glmnet(x, data$case, alpha = 1, family = "binomial")
qplot(x = cv.lasso$lambda, y = cv.lasso$nzero)
plot(cv.lasso)
coefpath(cv.lasso)
coefplot(cv.lasso)
cv.lasso$lambda.min
summary(cv.lasso$glmnet.fit)
## The lambda resulting in minimal error
## does not reduce any coefficients to zero
tidy(coef(cv.lasso, s = "lambda.min"))
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
glm_best_bic_lcp <-
    glm(
        case
        ~ logit1yrisk
        + logit1yrisk:longest_diam
        + logit1yrisk:any_poor_def
        # + logit1yrisk:max_lcp_score
        + adenopathy
        + any_upper
        + any_mixed
        # + max_lcp_score
        + 1,
        data = data,
        family = binomial(link = 'logit'),
        na.action = na.exclude,
    )
broom::tidy(glm_best_bic_lcp)
summary(glm_best_bic_lcp)
# AUC, AIC
glm_best_cp <-
    glm(
        case
        ~ logit1yrisk
        + logit1yrisk:diam_cat
        + logit1yrisk:adenopathy
        + logit1yrisk:any_growth
        + logit1yrisk:any_mixed
        + logit1yrisk:any_right_mid
        + logit1yrisk:any_spiculation
        + logit1yrisk:consolidation
        + logit1yrisk:any_margin_unab
        + diam_cat
        + consolidation
        + adenopathy
        + any_right_mid
        + any_lingula
        + any_spiculation
        + any_poor_def
        + 1,
        data = data,
        family = binomial(link = 'logit'),
        na.action = na.exclude,
    )
summary(glm_best_cp)


glm_best_bic <-
    glm(
        case
        ~ logit1yrisk
        + logit1yrisk:diam_cat
        + logit1yrisk:any_poor_def
        + logit1yrisk:max_lcp_score
        + diam_cat
        + adenopathy
        + any_upper
        + any_mixed
        + max_lcp_score
        + 1,
        data = data,
        family = binomial(link = 'logit'),
        na.action = na.exclude,
    )
tidy(coef(glm_best_bic))
coefplot(glm_best_bic)
glm_best_bic_lcp <-
    glm(
        case
        ~ logit1yrisk
        + logit1yrisk:longest_diam
        + logit1yrisk:any_poor_def
        # + logit1yrisk:max_lcp_score
        + adenopathy
        + any_upper
        + any_mixed
        # + max_lcp_score
        + 1,
        data = data,
        family = binomial(link = 'logit'),
        na.action = na.exclude,
    )
broom::tidy(glm_best_bic_lcp)
summary(glm_best_bic_lcp)
# AUC, AIC
glm_best_cp <-
    glm(
        case
        ~ logit1yrisk
        + logit1yrisk:diam_cat
        + logit1yrisk:adenopathy
        + logit1yrisk:any_growth
        + logit1yrisk:any_mixed
        + logit1yrisk:any_right_mid
        + logit1yrisk:any_spiculation
        + logit1yrisk:consolidation
        + logit1yrisk:any_margin_unab
        + diam_cat
        + consolidation
        + adenopathy
        + any_right_mid
        + any_lingula
        + any_spiculation
        + any_poor_def
        + 1,
        data = data,
        family = binomial(link = 'logit'),
        na.action = na.exclude,
    )
summary(glm_best_cp)


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
