# Libraries ----
library(readr)
# install.packages('here')
library(here)
library(dplyr)
# library(leaps)
library(stringr)
library(rsample)
library(broom)
library(recipes)
library(parsnip)
library(purrr)
# install.packages('remotes')
# remotes::install_github("rstudio/gt")
# remotes::install_github("tidymodels/tune")
# library(gt)
# library(tune)
library(yardstick)
library(ggplot2)
# install.packages("plotly")
library(plotly)

# Data ----
data <-
    read_rds(here("data/nlst_abn_lcp_pre.rds")) %>%
    mutate(case_at_next_screen = as.factor(case_at_next_screen),
           any_growth = as.factor(any_growth)) %>%
    filter(is.finite(longest_diam))

# Formulas ----

lcp_form <- as.formula(case_at_next_screen ~ max_lcp_score)

lcp_lcrat_form <- as.formula(case_at_next_screen ~ max_lcp_score + logit1yrisk)

lcrat_ct_form <- as.formula(
            case_at_next_screen
            ~ logit1yrisk
            + longest_diam
            # + emphysema
            # + consolidation
            # + adenopathy
            + any_growth # any_growth is a factor to not mislead in T0
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
    model_name = as_label(enquo(x)),
    actual = x$fit$y,
    predicted = x$fit$fitted.values,
    )
}

# Train models ----
# 1. LCP
lcp <- logit_mod %>% fit(lcp_form, data = data)
lcp
tidy_lcp <- tidy(lcp) %>% select(variable=term, coefficient=estimate, p_value=p.value)
# 2. LCRAT+CT
lcrat_ct <- logit_mod %>% fit(lcrat_ct_form, data = data)
tidy_lcrat_ct <- tidy(lcrat_ct)
# 3. LCP+LCRAT
lcp_lcrat <- logit_mod %>% fit(lcp_lcrat_form, data = data)
lcp_lcrat
tidy_lcp_lcrat <- tidy(lcp_lcrat)
# 4. LCP+LCRAT+CT (including nodule features) main effects & interactions
lcp_lcrat_ct <- logit_mod %>% fit(lcp_lcrat_ct_form, data = data)
lcp_lcrat_ct
tidy_lcp_lcrat_ct <- tidy(lcp_lcrat_ct)

pvals <- anova(lcp$fit, lcp_lcrat$fit, lcp_lcrat_ct$fit, test = "LRT")[5][[1]][2:3]
pvals <- c(NA, NA, pvals)

# Model table ----
# TODO add AUC to model table
get_model_table <- function (x) {
    model_name <- quo_name(enquo(x))
    glance(x$fit) %>%
        mutate(model = model_name,
               auc = roc_auc_vec(as.factor(x$fit$y),
                                 x$fit$fitted.values)
               )
}

test <- function (x) {
    quo_name(enquo(x))
}

test(lcp)

models <- c(lcrat_ct, lcp, lcp_lcrat, lcp_lcrat_ct)
fits <- c(lcrat_ct$fit, lcp$fit, lcp_lcrat$fit, lcp_lcrat_ct$fit)

mod_table <-
    bind_rows(
    get_model_table(lcrat_ct),
    get_model_table(lcp),
    get_model_table(lcp_lcrat),
    get_model_table(lcp_lcrat_ct)
    )

mod_table

head(data)

my_augment(lcp_lcrat_ct$fit)

        glance(lcrat_ct$fit) %>%
            mutate(model = "LCRAT+CT",
                   auc = roc_auc(my_augment(lcp_lcrat_ct),
                                 truth = as.factor(actual),
                                 predicted)[[3]]),
        glance(lcp$fit) %>%
            mutate(model = "LCP",
                   auc = roc_auc(my_augment(lcp_lcrat_ct),
                                 truth = as.factor(actual),
                                 predicted)[[3]]),
        glance(lcp_lcrat$fit) %>%
            mutate(model = "LCP+LCRAT",
                   auc = roc_auc(my_augment(lcp_lcrat_ct),
                                 truth = as.factor(actual),
                                 predicted)[[3]]),
        glance(lcp_lcrat_ct$fit) %>%
            mutate(model = "LCP+LCRAT+CT",
                   auc = roc_auc(my_augment(lcp_lcrat_ct),
                                 truth = as.factor(actual),
                                 predicted)[[3]])
) %>%
select(model, AIC, auc, deviance) %>%
    mutate(lrt_p_value = pvals) %>%
    identity()

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

# Variable table ----
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

# TODO make risk table as in https://academic.oup.com/jnci/article/111/9/996/5445482 first 3 columns only:
# columns 2 and 3 should absolute numbers and precentages

tail(lor_lcp)
lor_lcp <-
    my_augment(lcp) %>%
    arrange(desc(predicted)) %>%
    # recall (sensitivity)
    mutate(true_positives = cumsum(actual)) %>%
    mutate(positives = sum(actual)) %>%
    mutate(recall = true_positives / positives) %>%
    # precision (positive predictive value)
    mutate(predicted_positives = row_number()) %>%
    mutate(precision = true_positives / predicted_positives) %>%
    # fallout (false positive rate)
    mutate(false_positives = predicted_positives - true_positives) %>%
    mutate(negatives = n() - positives) %>%
    mutate(fallout = false_positives / negatives) %>%
    mutate(model = "lcp")

sum(data$case_at_next_screen)
lor_lcrat_ct <-
    my_augment(lcrat_ct) %>% 
    arrange(desc(predicted)) %>% 
    mutate(proportion_of_n = row_number() / n()) %>%
    mutate(proportion_of_y = cumsum(actual) / sum(actual)) %>%
    mutate(model = "lcrat_ct")
    
lor_lcp_lcrat <- 
    my_augment(lcp_lcrat) %>% 
    arrange(desc(predicted)) %>% 
    mutate(proportion_of_n = row_number() / n()) %>%
    mutate(proportion_of_y = cumsum(actual) / sum(actual)) %>%
    mutate(model = "lcp_lcrat")

lor_lcp_lcrat_ct <- 
    my_augment(lcp_lcrat_ct) %>% 
    arrange(desc(predicted)) %>% 
    mutate(proportion_of_n = row_number() / n()) %>%
    mutate(proportion_of_y = cumsum(actual) / sum(actual)) %>%
    mutate(model = "lcp_lcrat_ct")

all_lor <- bind_rows(lor_lcrat_ct, lor_lcp, lor_lcp_lcrat, lor_lcp_lcrat_ct)
all_lor %>% 
    ggplot() +
    aes(x=proportion_of_n, y = proportion_of_y, color = model) +
    geom_line(size = 1.2) +
    geom_segment(x = 0, xend = 1,
                 y = 0, yend =1,
                 linetype="dashed",
                 color="black") 
    # xlab("Proportion of participants in annual versus biennial screening") +
    # ylab("Sensitivity (Recall)") %>% 
    ggplotly(ggplot2::last_plot())

ggsave("lorenz.png")

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

merge
