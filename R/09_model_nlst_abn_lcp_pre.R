library(readr)
library(broom)
library(dplyr)
library(ggplot2)
library(here) 

data <- read_rds(here("data/nlst_abn_lcp_pre.rds"))

summarize_model <- function(model) {
    print(paste("AIC =", AIC(model)))
    print(paste("AUC =",
                pROC::auc(
                    data$case,
                    predict(model, type = "response")
                )))
    ggplot(tidy(model), aes(term, estimate)) +
        geom_col(aes(fill = term)) +
        coord_flip() +
        theme(legend.position = "none")
}


# 0. LCRAT only
glm_lcrat_only <-
    glm(case ~ log1yrisk - 1,
        data = data,
        family = binomial(link = 'log'),
        na.action = na.exclude)
summarize_model(glm_lcrat_only)


# 1. LCRAT + CT interacted
glm_interacted <-
    glm(
        case
        ~ log1yrisk
        + log1yrisk:diam.cat
        + log1yrisk:any.growth
        + log1yrisk:emphysema
        + log1yrisk:consolidation
        + log1yrisk:adenopathy
        + log1yrisk:any.upper
        + log1yrisk:I(any.right.mid == 1 | any.lingula == 1)
        + log1yrisk:any.mixed
        + log1yrisk:any.spiculation
        + log1yrisk:I(any.poor.def == 1 | any.margin.unab == 1)
        # + log1yrisk:max_lcp_score
        - 1,
        data = data,
        family = binomial(link = 'log'),
        na.action = na.exclude
    )
summary(glm_interacted)

glm_interacted_lcp <-
    glm(
        case
        ~ log1yrisk
        + log1yrisk:diam.cat
        + log1yrisk:any.growth
        + log1yrisk:emphysema
        + log1yrisk:consolidation
        + log1yrisk:adenopathy
        + log1yrisk:any.upper
        + log1yrisk:I(any.right.mid == 1 | any.lingula == 1)
        + log1yrisk:any.mixed
        + log1yrisk:any.spiculation
        + log1yrisk:I(any.poor.def == 1 | any.margin.unab == 1)
        + log1yrisk:max_lcp_score
        - 1,
        data = data,
        family = binomial(link = 'log'),
        na.action = na.exclude
    )

summary(glm_interacted_lcp)
lmtest::lrtest(glm_interacted, glm_interacted_lcp)

# 2. LCRAT + CT uninteracted
glm_uninteracted <-
    glm(
        case
        ~ log1yrisk
        + diam.cat
        + any.growth
        + emphysema
        + consolidation
        + adenopathy
        + any.upper
        + I(any.right.mid == 1 | any.lingula == 1)
        + any.mixed
        + any.spiculation
        + log1yrisk:I(any.poor.def == 1 | any.margin.unab == 1)
        # + log1yrisk:max_lcp_score
        - 1,
        data = data,
        family = binomial(link = 'log'),
        na.action = na.exclude
    )
summary(glm_uninteracted)

glm_uninteracted_lcp <-
    glm(
        case
        ~ log1yrisk
        + diam.cat
        + any.growth
        + emphysema
        + consolidation
        + adenopathy
        + any.upper
        + I(any.right.mid == 1 | any.lingula == 1)
        + any.mixed
        + any.spiculation
        + log1yrisk:I(any.poor.def == 1 | any.margin.unab == 1)
        + log1yrisk:max_lcp_score
        - 1,
        data = data,
        family = binomial(link = 'log'),
        na.action = na.exclude
    )

summary(glm_uninteracted_lcp)
lmtest::lrtest(glm_uninteracted, glm_uninteracted_lcp)

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
summarize_model(glm_lcrat_emph_cons)
       
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
summarize_model(glm_lcrat_emph_cons_aden)

# negative only
# only interactions
# check the n (33k negative screens, same dataset?)
# use log 1yr risk instead of prescreen risk
# use log link in glm function
glm_lcrat_pemph <-
    glm(
        case
        ~ log1yrisk
        # + log1yrisk:emphysema
        + log1yrisk:I(log(p_emph))
        - 1,
        data = data,
        family = binomial(link = 'log'),
        na.action = na.exclude
    )
summary(glm_lcrat_pemph)
summarize_model(glm_lcrat_pemph)

glm_lcrat_emph_pemph <-
    glm(
        case
        ~ log1yrisk
        + log1yrisk:emphysema
        + log1yrisk:p_emph
        - 1,
        data = data,
        family = binomial(link = 'log'),
        na.action = na.exclude
    )
summarize_model(glm_lcrat_emph_pemph)

glm_lcrat_lcp <-
    glm(
        case
        ~ log1yrisk
        + log1yrisk:max_lcp_score
        - 1,
        data = data,
        family = binomial(link = 'log'),
        na.action = na.exclude
    )
summarize_model(glm_lcrat_lcp)

glm_lcrat_lcp <-
    glm(
        case
        ~ log1yrisk:max_lcp_score
        - 1,
        data = data,
        family = binomial(link = 'log'),
        na.action = na.exclude
    )
summarize_model(glm_lcrat_lcp)
names(data)
glm_lcrat_lcp_emph <-
    glm(
        case
        ~ logit1yrisk
        + logit1yrisk:max_lcp_score
        + logit1yrisk:emphysema
        + logit1yrisk:consolidation
        + logit1yrisk:adenopathy
        - 1,
        data = data,
        family = binomial(link = 'logit'),
        na.action = na.exclude
    )
summary(glm_lcrat_lcp_emph)
summarize_model(glm_lcrat_lcp_emph)
names(data)

data <- data %>% 
    mutate(
        diam_cat = case_when(
            longest_diam == 0 ~ 1,
            longest_diam > 0 & longest_diam <= 5 ~ 2,
            longest_diam > 5 & longest_diam <= 7 ~ 3,
            longest_diam > 7 & longest_diam <= 10 ~ 4,
            longest_diam > 10 & longest_diam <= 13 ~ 5,
            longest_diam > 13 & longest_diam < 100 ~ 6
        ))
        
glm_base <-
    glm(
        case
        ~ logit1yrisk:diam_cat
        + logit1yrisk:any.growth
        + logit1yrisk:emphysema
        + logit1yrisk:consolidation
        + logit1yrisk:adenopathy
        + logit1yrisk:any.upper
        + logit1yrisk:I(any.right.mid == 1 | any.lingula == 1)
        + logit1yrisk:any.mixed
        + logit1yrisk:any.spiculation
        + logit1yrisk:I(any.poor.def == 1 | any.margin.unab == 1)
        # + logit1yrisk:max_lcp_score
        - 1,
        data = data,
        family = binomial(link = 'logit')
    )

summary(glm_base)
tidy(glm_base)
tail(data$max_lcp_score)
sum(is.na(data$max_lcp_score))
glm_fit_lcrat_lcp_emph_pemph <-
    glm(
        case
        ~ log1yrisk
        + log1yrisk:max_lcp_score
        + log1yrisk:emphysema
        + log1yrisk:I(log(p_emph))
        - 1,
        data = data,
        family = binomial(link = 'log'),
        na.action = na.exclude
    )
# summarize_model(glm_lcrat_lcp_emph)

x = c(
    "log1yrisk",
    "emphysema",
    "consolidation",
    "adenopathy"
    )

h2o_glm_lcrat_ct <- h2o.glm(
    x = x,
    y = y,
    interaction_pairs = list(
        c("log1yrisk", "emphysema"),
        c("log1yrisk", "consolidation"),
        c("log1yrisk", "adenopathy")
    ),
    training_frame = train,
    family = "binomial",
    lambda_search = TRUE
)
h2o_summarize_model(h2o_glm_lcrat_ct)

x = c(
    "log1yrisk",
    "emphysema",
    "consolidation",
    "adenopathy",
    "max_lcp_score"
    )

h2o_glm_lcrat_ct_lcp <- h2o.glm(
    x = x,
    y = y,
    interaction_pairs = list(
        c("log1yrisk", "emphysema"),
        c("log1yrisk", "consolidation"),
        c("log1yrisk", "adenopathy"),
        c("log1yrisk", "max_lcp_score")
    ),
    training_frame = train,
    family = "binomial",
    lambda_search = TRUE
)
h2o_summarize_model(h2o_glm_lcrat_ct_lcp)


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
h2o_summarize_model(h2o_glm_lcrat_ct_all)

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
h2o_summarize_model(h2o_glm_lcrat_ct_all)

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


# TODO 
# for nested models, use lrt (get p-value)
# for non-nested models, look for decrease in deviance, AIC, BIC 
# First step: LCDRAT + lcp_score (composite of nodule features) and emphysema (case ~ prerisk + prerisk:lcp + prerisk:emph + prerisk:adenopathy + prerisk:consolidation) Negative screens only with lcp score and for positive screens nodule features (any spiculation, longest diameter)
# Second step: build top "only interaction terms" model
# Second step: build top "no interaction terms" model
# Second step: build "a mix interaction terms"
# both main effect and interaction
