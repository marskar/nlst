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
library(gt)
#remotes::install_github("rstudio/gt")
# remotes::install_github("tidymodels/tune")
library(tune)
library(yardstick)
library(ggplot2)

# Data ----
data <-
    read_rds(here("data/nlst_abn_lcp_pre.rds")) %>%
    mutate(case_at_next_screen = as.factor(case_at_next_screen)) %>% 
    filter(is.finite(longest_diam))

# Formulas ----

# bss_data <-
#     read_rds(here("data/nlst_abn_lcp_pre.rds")) %>%
#     select(
#         pid,
#         case_at_next_screen,
#         logit1yrisk,
#         diam_cat,
#         any_growth,
#         emphysema,
#         consolidation,
#         adenopathy,
#         any_upper,
#         any_right_mid,
#         any_lingula,
#         any_mixed,
#         any_spiculation,
#         any_poor_def,
#         any_margin_unab
#            )

# names(bss_data)
# Model selection decides whether variable is a main effect or interaction

# get_best_bic_formula <- function(form, data, yname) {
#     regsumm <-
#         summary(leaps::regsubsets(
#             x = form,
#             data = bss_data,
#             nvmax = 50,
#             method = "exhaustive"
#         ))
#     coefs <- regsumm$which[which.min(regsumm$bic), -1]
#     predictors <- names(which(coefs == TRUE)) %>%
#         str_remove("\\d$") %>%
#         unique() %>%
#         paste(collapse = "+")
#     as.formula(paste(yname, "~", predictors))
# }

# lcrat_ct_form <- get_best_bic_formula(
#     case_at_next_screen
#     ~ logit1yrisk * . + 1,
#     bss_data,
#     "case_at_next_screen")
# lcrat_ct_form


# lcrat_diam_form <- as.formula(
#             case_at_next_screen
#             ~ logit1yrisk*diam_cat
#             + any_growth
#             + emphysema
#             + consolidation
#             + adenopathy
#             + any_upper
#             + any_right_mid
#             + any_lingula
#             + any_mixed
#             + any_spiculation
#             + any_poor_def
#             + any_margin_unab
# )
# 
# lcrat_grow_form <- as.formula(
#             case_at_next_screen
#             ~ logit1yrisk*any_growth
#             + diam_cat
#             + emphysema
#             + consolidation
#             + adenopathy
#             + any_upper
#             + any_right_mid
#             + any_lingula
#             + any_mixed
#             + any_spiculation
#             + any_poor_def
#             + any_margin_unab
# )
# 
# lcrat_emph_form <- as.formula(
#             case_at_next_screen
#             ~ logit1yrisk*emphysema
#             + diam_cat
#             + any_growth
#             + consolidation
#             + adenopathy
#             + any_upper
#             + any_right_mid
#             + any_lingula
#             + any_mixed
#             + any_spiculation
#             + any_poor_def
#             + any_margin_unab
# )
# 
# lcrat_cons_form <- as.formula(
#             case_at_next_screen
#             ~ logit1yrisk*consolidation
#             + diam_cat
#             + any_growth
#             + emphysema
#             + adenopathy
#             + any_upper
#             + any_right_mid
#             + any_lingula
#             + any_mixed
#             + any_spiculation
#             + any_poor_def
#             + any_margin_unab
# )
# 
# lcrat_aden_form <- as.formula(
#             case_at_next_screen
#             ~ logit1yrisk*adenopathy
#             + diam_cat
#             + any_growth
#             + emphysema
#             + consolidation
#             + any_upper
#             + any_right_mid
#             + any_lingula
#             + any_mixed
#             + any_spiculation
#             + any_poor_def
#             + any_margin_unab
# )
# 
# lcrat_uppr_form <- as.formula(
#             case_at_next_screen
#             ~ logit1yrisk*any_upper
#             + diam_cat
#             + any_growth
#             + emphysema
#             + consolidation
#             + adenopathy
#             + any_right_mid
#             + any_lingula
#             + any_mixed
#             + any_spiculation
#             + any_poor_def
#             + any_margin_unab
# )
# 
# lcrat_rmid_form <- as.formula(
#             case_at_next_screen
#             ~ logit1yrisk*any_right_mid
#             + diam_cat
#             + any_growth
#             + emphysema
#             + consolidation
#             + adenopathy
#             + any_upper
#             + any_lingula
#             + any_mixed
#             + any_spiculation
#             + any_poor_def
#             + any_margin_unab
# )
# 
# lcrat_ling_form <- as.formula(
#             case_at_next_screen
#             ~ logit1yrisk*any_lingula
#             + diam_cat
#             + any_growth
#             + emphysema
#             + consolidation
#             + adenopathy
#             + any_upper
#             + any_right_mid
#             + any_mixed
#             + any_spiculation
#             + any_poor_def
#             + any_margin_unab
# )
# 
# lcrat_mixd_form <- as.formula(
#             case_at_next_screen
#             ~ logit1yrisk*any_mixed
#             + diam_cat
#             + any_growth
#             + emphysema
#             + consolidation
#             + adenopathy
#             + any_upper
#             + any_right_mid
#             + any_lingula
#             + any_spiculation
#             + any_poor_def
#             + any_margin_unab
# )
# 
# lcrat_spic_form <- as.formula(
#             case_at_next_screen
#             ~ logit1yrisk*any_spiculation
#             + diam_cat
#             + any_growth
#             + emphysema
#             + consolidation
#             + adenopathy
#             + any_upper
#             + any_right_mid
#             + any_lingula
#             + any_mixed
#             + any_poor_def
#             + any_margin_unab
# )
# 
# lcrat_poor_form <- as.formula(
#             case_at_next_screen
#             ~ logit1yrisk*any_poor_def
#             + diam_cat
#             + any_growth
#             + emphysema
#             + consolidation
#             + adenopathy
#             + any_upper
#             + any_right_mid
#             + any_mixed
#             + any_lingula
#             + any_spiculation
#             + any_margin_unab
# )
# 
# lcrat_marg_form <- as.formula(
#             case_at_next_screen
#             ~ logit1yrisk*any_margin_unab
#             + diam_cat
#             + any_growth
#             + emphysema
#             + consolidation
#             + adenopathy
#             + any_upper
#             + any_right_mid
#             + any_mixed
#             + any_lingula
#             + any_spiculation
#             + any_poor_def
# )

lcp_form <- as.formula(case_at_next_screen ~ max_lcp_score)

lcp_lcrat_form <- as.formula(case_at_next_screen ~ max_lcp_score + logit1yrisk)

lcrat_ct_form <- as.formula(
            case_at_next_screen
            ~ logit1yrisk
            + longest_diam
            # + emphysema
            # + consolidation
            # + adenopathy
            + any_growth # TODO any_growth may be misleading in T0: turn into factor
            + any_upper
            + any_right_mid
            + any_lingula
            + any_mixed
            + any_spiculation
            + any_poor_def
            + any_margin_unab
)

lcp_lcrat_ct_form <- update(lcrat_ct_form, ~ . + max_lcp_score)

# Model setup ----
logit_mod <-
    logistic_reg(mode = "classification") %>%
    set_engine("glm")

# Refit to all data to get one model & get summary for each model
# TODO 1. Send model summary tables for all three models (Table 1 for paper)

# Model output functions ----
my_glance <- function(x, ...) {
    tibble(
        # model_name = as_label(enquo(x)),
        null_deviance = x$fit$null.deviance,
        null = x$fit$df.null,
        aic = x$fit$aic,
        deviance = x$fit$deviance,
        residual = x$fit$df.residual
    )
}

# TODO put in splines, try transformation of max_lcp_score
my_augment <- function(x, ...) {
    tibble(
    # model_name = as_label(enquo(x)),
    actual = x$fit$y,
    predicted = x$fit$fitted.values,
    residuals = x$fit$residuals,
    )
} 

# Train models ----
    # 1. LCP
lcp <- logit_mod %>% fit(lcp_form, data = data)
lcp
tidy_lcp <- tidy(lcp) %>% select(variable=term, coefficient=estimate, p_value=p.value)


    # 2. LCRAT+CT 
lcrat_ct <- logit_mod %>% fit(lcrat_ct_form, data = data)
# lcrat_grow <- logit_mod %>% fit(lcrat_grow_form, data = data)
# lcrat_diam <- logit_mod %>% fit(lcrat_diam_form, data = data)
# lcrat_emph <- logit_mod %>% fit(lcrat_emph_form, data = data)
# lcrat_cons <- logit_mod %>% fit(lcrat_cons_form, data = data)
# lcrat_aden <- logit_mod %>% fit(lcrat_aden_form, data = data)
# lcrat_uppr <- logit_mod %>% fit(lcrat_uppr_form, data = data)
# lcrat_rmid <- logit_mod %>% fit(lcrat_rmid_form, data = data)
# lcrat_ling <- logit_mod %>% fit(lcrat_ling_form, data = data)
# lcrat_spic <- logit_mod %>% fit(lcrat_spic_form, data = data)
# lcrat_mixd <- logit_mod %>% fit(lcrat_mixd_form, data = data)
# lcrat_poor <- logit_mod %>% fit(lcrat_poor_form, data = data)
# lcrat_marg <- logit_mod %>% fit(lcrat_marg_form, data = data)

tidy_lcrat_ct <- tidy(lcrat_ct)
# tidy_lcrat_grow <- tidy(lcrat_grow) %>% mutate(model = "grow")
# tidy_lcrat_diam <- tidy(lcrat_diam) %>% mutate(model = "diam")
# tidy_lcrat_emph <- tidy(lcrat_emph) %>% mutate(model = "emph")
# tidy_lcrat_cons <- tidy(lcrat_cons) %>% mutate(model = "cons")
# tidy_lcrat_aden <- tidy(lcrat_aden) %>% mutate(model = "aden")
# tidy_lcrat_uppr <- tidy(lcrat_uppr) %>% mutate(model = "uppr")
# tidy_lcrat_rmid <- tidy(lcrat_rmid) %>% mutate(model = "rmid")
# tidy_lcrat_ling <- tidy(lcrat_ling) %>% mutate(model = "ling")
# tidy_lcrat_spic <- tidy(lcrat_spic) %>% mutate(model = "spic")
# tidy_lcrat_mixd <- tidy(lcrat_mixd) %>% mutate(model = "mixd")
# tidy_lcrat_poor <- tidy(lcrat_poor) %>% mutate(model = "poor")
# tidy_lcrat_marg <- tidy(lcrat_marg) %>% mutate(model = "marg")

# tidy_all <- rbind(
#     tidy_lcrat_grow,
#     tidy_lcrat_diam,
#     tidy_lcrat_emph,
#     tidy_lcrat_cons,
#     tidy_lcrat_aden,
#     tidy_lcrat_uppr,
#     tidy_lcrat_rmid,
#     tidy_lcrat_ling,
#     tidy_lcrat_spic,
#     tidy_lcrat_mixd,
#     tidy_lcrat_poor,
#     tidy_lcrat_marg
# )

    # 2. LCRAT+CT diam
lcrat_diam <- logit_mod %>% fit(lcrat_diam_form, data = data)
lcrat_diam
tidy_lcrat_diam <- tidy(lcrat_diam)
    # 2. LCRAT+CT grow
lcrat_grow <- logit_mod %>% fit(lcrat_grow_form, data = data)
lcrat_grow
tidy_lcrat_grow <- tidy(lcrat_grow)
    # 3. LCP+LCRAT
lcp_lcrat <- logit_mod %>% fit(lcp_lcrat_form, data = data)
lcp_lcrat
tidy_lcp_lcrat <- tidy(lcp_lcrat)
    # 4. LCP+LCRAT+CT (including nodule features) main effects & interactions
lcp_lcrat_ct <- logit_mod %>% fit(lcp_lcrat_ct_form, data = data)
lcp_lcrat_ct
tidy_lcp_lcrat_ct <- tidy(lcp_lcrat_ct)
lcp$fit

pvals <- anova(lcp$fit, lcp_lcrat$fit, lcp_lcrat_ct$fit, test = "LRT")[5][[1]][2:3]
pvals <- c(NA, NA, pvals)

glance(lcp$fit)
# Glance ----
mod_table <- 
    bind_rows(
    glance(lcrat_ct$fit) %>% mutate(model = "LCRAT+CT"),
    glance(lcp$fit) %>% mutate(model = "LCP"),
    glance(lcp_lcrat$fit) %>% mutate(model = "LCP+LCRAT"),
    glance(lcp_lcrat_ct$fit) %>% mutate(model = "LCP+LCRAT+CT")
) %>% 
    select(model, AIC, deviance) %>%
    mutate(lrt_p_value = pvals) %>% 
    gt() %>%
    tab_header(
        title = "Nested Models"
    ) %>%
    fmt_number(
        columns = vars(AIC, deviance),
        decimals = 0
    ) %>% 
    fmt_scientific(
        columns = vars(lrt_p_value),
    ) %>% 
    fmt_missing(
        columns = vars(lrt_p_value),
        missing_text = ""
    ) %>%
    gt::cols_label(
        model="Model",
        deviance="Deviance",
        lrt_p_value="LRT p-value"
    )
mod_table

obs_table <- augment(lcp$fit) 
tidy(lcp)
tidy(lcp$fit)
var_table <- 
    bind_rows(
    tidy(lcp) %>% mutate(model = "LCP"),
    tidy(lcrat_ct) %>% mutate(model = "LCRAT+CT"),
    tidy(lcp_lcrat) %>% mutate(model = "LCP+LCRAT"),
    tidy(lcp_lcrat_ct) %>% mutate(model = "LCP+LCRAT+CT")
) %>% 
    select(model, variable = term, coefficient = estimate, standard_error = std.error, p_value = p.value) %>% 
    gt() %>%
    tab_header(
        title = "Variables",
    ) %>%
    fmt_scientific(
        columns = vars(p_value),
    ) %>% 
    fmt_number(
        columns = vars(coefficient, standard_error),
        decimals = 3
    ) %>%
    gt::cols_label(
        model="Model",
        variable="Variable",
        coefficient="Coefficient",
        standard_error="Standard Error",
        p_value="p-value"
    )

var_table
remotes::install_github("rstudio/gt")

augment(lcp$fit)
my_augment(lcp)


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
lcp_gain <- 
    gain_curve(true_pred, truth = truth, lcp) %>%
        mutate(model = "lcp_lcrat")
lcrat_gain <- 
    gain_curve(true_pred, truth = truth, lcrat_ct) %>%
        mutate(model = "lcp_lcrat")
lcrat_gain <- 
    gain_curve(true_pred, truth = truth, lcp_lcrat) %>%
        mutate(model = "lcp_lcrat")
gain_curve(true_pred, truth = truth, lcp_lcrat_ct) %>%
    mutate(model = "lcp_lcrat")

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
