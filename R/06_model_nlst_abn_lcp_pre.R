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

# 3. LCRAT + CT hybrid
glm_hybrid <-
    glm(
        case
        ~ log1yrisk
        + log1yrisk*diam.cat
        + log1yrisk*any.growth
        + log1yrisk*emphysema
        + log1yrisk*consolidation
        + log1yrisk*adenopathy
        + log1yrisk*any.upper
        + log1yrisk*I(any.right.mid == 1 | any.lingula == 1)
        + log1yrisk*any.mixed
        + log1yrisk*any.spiculation
        + log1yrisk*I(any.poor.def == 1 | any.margin.unab == 1)
        # + log1yrisk*max_lcp_score
        - 1,
        data = data,
        family = binomial(link = 'log'),
        na.action = na.exclude
    )
glimpse(tidy(glm_interacted))

glm_hybrid_lcp <-
    glm(
        case
        ~ log1yrisk
        + log1yrisk*diam.cat
        + log1yrisk*any.growth
        + log1yrisk*emphysema
        + log1yrisk*consolidation
        + log1yrisk*adenopathy
        + log1yrisk*any.upper
        + log1yrisk*I(any.right.mid == 1 | any.lingula == 1)
        + log1yrisk*any.mixed
        + log1yrisk*any.spiculation
        + log1yrisk*I(any.poor.def == 1 | any.margin.unab == 1)
        + log1yrisk*max_lcp_score
        - 1,
        data = data,
        family = binomial(link = 'log'),
        na.action = na.exclude
    )
summary(glm_hybrid_lcp)
glimpse(tidy(glm_hybrid_lcp))
lmtest::lrtest(glm_hybrid, glm_hybrid_lcp)

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

        
glm_logit <-
    glm(
        case
        ~ logit1yrisk:diam.cat
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

summary(glm_logit)
tidy(glm_logit)
tail(data$max_lcp_score)
str(data)
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


# TODO 
# for nested models, use lrt (get p-value)
# for non-nested models, look for decrease in deviance, AIC, BIC 
# First step: LCDRAT + lcp_score (composite of nodule features) and emphysema (case ~ prerisk + prerisk:lcp + prerisk:emph + prerisk:adenopathy + prerisk:consolidation) Negative screens only with lcp score and for positive screens nodule features (any spiculation, longest diameter)
# Second step: build top "only interaction terms" model
# Second step: build top "no interaction terms" model
# Second step: build "a mix interaction terms"
# both main effect and interaction
