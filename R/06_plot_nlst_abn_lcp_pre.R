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

data <-
    read_rds(here("data/nlst_abn_lcp_pre.rds")) %>%
    mutate(case_at_next_screen = as.factor(case_at_next_screen))


logit_mod <- logistic_reg(mode = "classification") %>% set_engine("glm")
ctrl <- control_resamples(save_pred = TRUE)
cv_train <- vfold_cv(data = data, v = 5, strata = "case_at_next_screen")

# Fit and predict
lcp_pred <-
    fit_resamples(
        case_at_next_screen ~ max_lcp_score,
        logit_mod,
        resamples = cv_train,
        control = ctrl
    ) %>%
    collect_predictions() %>%
    roc_curve(truth=case_at_next_screen, estimate=.pred_1) %>% 
    filter(is.finite(.threshold)) %>% 
    mutate(max_pred = max(.threshold, na.rm=TRUE)) %>% 
    mutate(prop_max_pred = .threshold / max_pred) %>% 
    mutate(model = "lcp")

lcrat_pred <-
    fit_resamples(
        case_at_next_screen
        ~ logit1yrisk
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
    mutate(prop_max_pred = .threshold / max_pred) %>% 
    mutate(model = "lcrat")

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
    mutate(prop_max_pred = .threshold / max_pred) %>% 
    mutate(model = "lcp_lcrat")

lcp_pred %>% autoplot()
lcrat_pred %>% autoplot()
lcp_lcrat_pred %>% autoplot()
lcp_roc %>% glimpse()

all_pred <- bind_rows(lcp_pred, lcrat_pred, lcp_lcrat_pred)
all_pred %>% glimpse()

all_pred %>% 
    ggplot() +
    aes(x=prop_max_pred, y = sensitivity, color = model) +
    geom_line() +
    xlab("Proportion of participants in biennial versus annual screening")

all_pred %>% 
    ggplot() +
    aes(x=1 - specificity, y = sensitivity, color = model) +
    geom_line() +
    geom_abline(linetype="dashed") +
    xlab("1 - specificity")
