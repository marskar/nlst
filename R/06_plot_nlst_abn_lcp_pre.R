# Libraries ----
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

# Data ----
data <-
    read_rds(here("data/nlst_abn_lcp_pre.rds")) %>%
    mutate(case_at_next_screen = as.factor(case_at_next_screen))


# Model settings ----
logit_mod <- logistic_reg(mode = "classification") %>% set_engine("glm")
ctrl <- control_resamples(save_pred = TRUE)
cv_train <- vfold_cv(data = data, v = 5, strata = "case_at_next_screen")

# Fit and predict ----
# Refit to all data to get one model & get summary for each model
# TODO 1. Send model summary tables for all three models (Table 1 for paper)
    # 1. LCP
    # 2. LCRAT+CT (including nodule features) main effects & interactions
    # 3. LCP+LCRAT
    # 4. LCP+LCRAT+CT (including nodule features) main effects & interactions

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
