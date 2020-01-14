# Libraries ----
library(readr)
library(here)
library(dplyr)
library(leaps)
library(stringr)
library(rsample)
library(broom)
library(recipes)
library(parsnip)
# remotes::install_github("tidymodels/tune")
library(tune)
library(yardstick)
library(ggplot2)

# Data ----
data <-
    read_rds(here("data/nlst_abn_lcp_pre.rds")) %>%
    mutate(case_at_next_screen = as.factor(case_at_next_screen))

# Formulas ----
lcp_form <- as.formula(case_at_next_screen ~ max_lcp_score)

lcp_lcrat_form <- as.formula(case_at_next_screen ~ max_lcp_score + logit1yrisk)

bss_data <-
    read_rds(here("data/nlst_abn_lcp_pre.rds")) %>%
    select(
        pid,
        case_at_next_screen,
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
        any_margin_unab
           )

names(bss_data)
# Model selection decides whether variable is a main effect or interaction

get_best_bic_formula <- function(form, data, yname) {
    regsumm <-
        summary(leaps::regsubsets(
            x = form,
            data = bss_data,
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

lcrat_ct_form <- get_best_bic_formula(
    case_at_next_screen
    ~ logit1yrisk * . + 1,
    bss_data,
    "case_at_next_screen")
lcrat_ct_form

lcp_lcrat_ct_form <- update(lcrat_ct_form, ~ . + max_lcp_score)

# Model setup ----
logit_mod <- logistic_reg(mode = "classification") %>% set_engine("glm")

# Refit to all data to get one model & get summary for each model
# TODO 1. Send model summary tables for all three models (Table 1 for paper)

# Train models ----
    # 1. LCP
lcp <- logit_mod %>% fit(lcp_form, data = data)
lcp
tidy(lcp)
    # 2. LCRAT+CT (including nodule features) main effects & interactions
lcrat_ct <- logit_mod %>% fit(lcrat_ct_form, data = data)
lcrat_ct
tidy(lcrat_ct)
    # 3. LCP+LCRAT
lcp_lcrat <- logit_mod %>% fit(lcp_lcrat_form, data = data)
lcp_lcrat
tidy(lcp_lcrat)
    # 4. LCP+LCRAT+CT (including nodule features) main effects & interactions
lcp_lcrat_ct <- logit_mod %>% fit(lcp_lcrat_ct_form, data = data)
lcp_lcrat_ct
tidy(lcp_lcrat_ct)

# Predict ----
    # 1. LCP
lcp_pred <- lcp %>% predict(new_data = data,
                            type = "prob")
    # 2. LCRAT+CT (including nodule features) main effects & interactions
lcrat_ct_pred <- lcrat_ct %>%
    predict(new_data = data, type = "prob")
lcrat_ct_pred <- lcrat_ct %>%
    predict(new_data = data, type = "prob")
    # 3. LCP+LCRAT
lcp_lcrat_pred <- lcp_lcrat %>%
    predict(new_data = data, type = "prob")
    # 4. LCP+LCRAT+CT (including nodule features) main effects & interactions
lcp_lcrat_ct_pred <- lcp_lcrat_ct %>%
    predict(new_data = data, type = "prob")

# df = pred + true ----

true_pred <- bind_cols(truth = data$case_at_next_screen,
                       lcp = lcp_pred$.pred_1,
                       lcrat_ct = lcrat_ct_pred$.pred_1,
                       lcp_lcrat = lcp_lcrat_pred$.pred_1,
                       lcp_lcrat_ct = lcp_lcrat_ct_pred$.pred_1)

true_pred %>% glimpse()

# Plot ----

options(yardstick.event_first = FALSE)
gain_curve(true_pred, truth = truth, lcp) %>% autoplot()
gain_curve(true_pred, truth = truth, lcrat_ct) %>% autoplot()
gain_curve(true_pred, truth = truth, lcp_lcrat) %>% autoplot()
gain_curve(true_pred, truth = truth, lcp_lcrat_ct) %>% autoplot()

# Cross-validation ----
ctrl <- control_resamples(save_pred = TRUE)
cv_train <- vfold_cv(data = data, v = 5, strata = "case_at_next_screen")

# Fit and predict ----

# Calculate sensitivity and get predictions 

# TODO 2. Send table of sensitivities and percent exceeding threshold for all three models

# 1. LCP

lcp_mod <-
    fit_resamples(
        case_at_next_screen ~ max_lcp_score,
        logit_mod,
        resamples = cv_train,
        control = ctrl
    ) %>%
    identity()
summary(lcp_mod)
    collect_predictions()%>%
    identity()

    roc_curve(truth=case_at_next_screen, estimate=.pred_1) %>% 
    filter(is.finite(.threshold)) %>% 
    mutate(max_pred = max(.threshold, na.rm=TRUE)) %>% 
    mutate(min_pred = min(.threshold, na.rm=TRUE)) %>% 
    mutate(min_max_diff = max_pred - min_pred) %>% 
    mutate(min_max_pred = (.threshold - min_pred) / min_max_diff) %>% 
    mutate(model = "lcp")

# 2. LCRAT+CT (including nodule features) main effects & interactions
lcrat_pred <-
    fit_resamples(
        case_at_next_screen
        ~ logit1yrisk
        + emphysema
        + consolidation,
# Put in features based on leaps best subset selection without lcp
        logit_mod,
        resamples = cv_train,
        control = ctrl
    ) %>%
    collect_predictions() %>%
    roc_curve(truth=case_at_next_screen, estimate=.pred_1) %>% 
    filter(is.finite(.threshold)) %>% 
    mutate(max_pred = max(.threshold, na.rm=TRUE)) %>% 
    mutate(min_pred = min(.threshold, na.rm=TRUE)) %>% 
    mutate(min_max_diff = max_pred - min_pred) %>% 
    mutate(min_max_pred = (.threshold - min_pred) / min_max_diff) %>% 
    mutate(model = "lcrat")

# 3. LCP+LCRAT

# 4. LCP+LCRAT+CT (including nodule features) main effects & interactions
# Send model summaries 
lcp_lcrat_pred <-
    fit_resamples(
        case_at_next_screen
        ~ max_lcp_score
        + logit1yrisk
        + emphysema
        + consolidation,
        logit_mod,
        resamples = cv_train,
        control = ctrl
    ) %>%
    collect_predictions() %>%
    roc_curve(truth=case_at_next_screen, estimate=.pred_1) %>% 
    filter(is.finite(.threshold)) %>% 
    mutate(max_pred = max(.threshold, na.rm=TRUE)) %>% 
    mutate(min_pred = min(.threshold, na.rm=TRUE)) %>% 
    mutate(min_max_diff = max_pred - min_pred) %>% 
    mutate(min_max_pred = (.threshold - min_pred) / min_max_diff) %>% 
    mutate(model = "lcp_lcrat")

lcp_pred %>% autoplot()
lcrat_pred %>% autoplot()
lcp_lcrat_pred %>% autoplot()
lcp_pred %>% glimpse()
lcrat_pred %>% glimpse()

all_pred <- bind_rows(lcp_pred, lcrat_pred, lcp_lcrat_pred)
all_pred %>% glimpse()

all_pred %>% 
    ggplot() +
    aes(x=min_max_pred, y = sensitivity, color = model) +
    geom_line(size = 1.2) +
    geom_segment(x = 0, xend = 1,
                 y = 0, yend =1,
                 linetype="dashed",
                 color="black") +
    xlab("Proportion of participants in annual versus biennial screening") +
    ylab("Sensitivity (Recall)")

ggsave("lorenz.png")

# ROC curve ----

all_pred %>% 
    ggplot() +
    aes(x=1 - specificity, y = sensitivity, color = model) +
    geom_line() +
    geom_segment(x = 0, xend = 1, y = 0, yend =1, linetype="dashed", color="black") +
    xlab("1 - specificity")

ggsave("roc.png")

# Gain curve ----

options(yardstick.event_first = FALSE)

lcp_gain <- fit_resamples(
    case_at_next_screen
    ~ max_lcp_score,
    logit_mod,
    resamples = cv_train,
    control = ctrl
) %>%
    collect_predictions() %>%
    gain_curve(truth=case_at_next_screen, estimate=.pred_1) %>%
    mutate(model = "lcp")

lcrat_gain <- fit_resamples(
    case_at_next_screen
    ~ logit1yrisk
    + emphysema
    + consolidation,
    logit_mod,
    resamples = cv_train,
    control = ctrl
) %>%
    collect_predictions() %>%
    gain_curve(truth=case_at_next_screen, estimate=.pred_1) %>% 
    mutate(model = "lcrat")

lcp_lcrat_gain <- fit_resamples(
    case_at_next_screen
    ~ max_lcp_score
    + logit1yrisk
    + emphysema
    + consolidation,
    logit_mod,
    resamples = cv_train,
    control = ctrl
) %>%
    collect_predictions() %>%
    gain_curve(truth=case_at_next_screen, estimate=.pred_1) %>% 
    mutate(model = "lcp_lcrat")

lcp_gain %>% autoplot()
lcrat_gain %>% autoplot()
lcp_lcrat_gain %>% autoplot()

all_gain <- bind_rows(lcp_gain, lcrat_gain, lcp_lcrat_gain)

all_gain %>% 
    ggplot() +
    aes(x=.percent_tested, y = .percent_found, color = model) +
    geom_line(size = 1.2) +
    geom_segment(x = 0, xend = 100,
                 y = 0, yend = 100,
                 linetype="dashed",
                 color="black") +
    xlab("Percent tested") +
    ylab("Percent found")

ggsave("gain.png")
