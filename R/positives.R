#### Setup ####

# List packages to be loaded (and installed if needed)
packages <-
    c(
        "readr",
        "lmtest",
        "here",
        "dplyr",
        "ggplot2",
        "survival",
        "gmodels",
        "devtools",
        "geepack",
        "MESS",
        "psych",
        "Hmisc",
        "glmnet",
        "boot",
        "pROC"
    )

# List packages that are not installed
not_installed <-
    packages[!(packages %in% installed.packages()[, "Package"])]

# Install packages that are not installed
if (length(not_installed))
    install.packages(not_installed)

# Load all packages
lapply(packages, require, character.only = TRUE)


data_screen_abn_pos <-
    readRDS(here('data/data_screen_abn_pos.rds')) %>%
    # Make a categorical variable for diameter
    mutate(
        diam_cat = case_when(
            longest.diam == 0 ~ 1,
            longest.diam > 0 & longest.diam <= 5 ~ 2,
            longest.diam > 5 & longest.diam <= 7 ~ 3,
            longest.diam > 7 & longest.diam <= 10 ~ 4,
            longest.diam > 10 & longest.diam <= 13 ~ 5,
            longest.diam > 13 & longest.diam < 100 ~ 6
        )
    ) %>%
    mutate(diam_cat = factor(
        diam_cat,
        levels = c(1:6),
        labels = c("0", "4-5", "6-7", "8-10", "11-13", "14+")
    ))

# Run the main models
# Fit an overall model with one exponent

glm_screen_pos <-
    glm(case ~ log1yrisk - 1,
        data = data_screen_abn_pos,
        family = binomial(link = 'log'))


# Model with abnormalities
glm_screen_pos_abn_log <-
    glm(
        case ~ log1yrisk:diam_cat + log1yrisk:any.growth + log1yrisk:any.upper +
            log1yrisk:I(any.right.mid == 1 |
                            any.lingula == 1) + log1yrisk:any.mixed +
            log1yrisk:any.spiculation + log1yrisk:I(any.poor.def ==
                                                        1 | any.margin.unab == 1) - 1,
        data = data_screen_abn_pos,
        family = binomial(link = 'log')
    )
# data_screen_abn_pos$post_risk_abn_pos <- fitted.values(glm_screen_pos_abn_log)

abn <- readRDS(here('data/abn_lrads_merged.rds'))
nlst_ct <- readRDS(here('data/nlst_ct.rds'))

# Will need this vector for exploratory analysis of abnormalities (24 Jan 2017 - also create interaction vectors)
abnlist <- abnlist.pos <- c(names(abn)[5:32], "diam_cat") # omits longest.diam and any.nodule, but includes diam.cat which subsumes both. 12 Oct 2018 now includes any.new.nodule
abnlist.pos.int <-
    lapply(abnlist.pos, function(x) {
        substitute(logit1yrisk:i, list(i = as.name(x)))
    })

# Make dataset of unique individuals for descriptive table of false-positives
all.subj.pos <-
    dplyr::filter(nlst_ct, pid %in% data_screen_abn_pos$pid)
all.subj.pos <-
    mutate(
        all.subj.pos,
        age.cat = as.factor(ifelse(
            age >= 55 & age < 60,
            "55-59",
            ifelse(
                age >= 60 & age < 65,
                "60-64",
                ifelse(age >= 65 &
                           age < 70, "65-69", ifelse(age >= 70 & age < 75, "70-74", NA))
            )
        )),
        qtyears.cat = as.factor(ifelse(
            qtyears == 0,
            "Current smoker",
            ifelse(
                qtyears > 0 & qtyears <= 5,
                "1-5",
                ifelse(
                    qtyears > 5 &
                        qtyears <= 10,
                    "6-10",
                    ifelse(qtyears > 10 & qtyears < 99, "11 or more", NA)
                )
            )
        )),
        bmi.cat = as.factor(ifelse(
            bmi > 0 &
                bmi < 18.5,
            "Underweight",
            ifelse(
                bmi >= 18.5 & bmi < 25,
                "Normal",
                ifelse(bmi >= 25 &
                           bmi < 30, "Overweight", ifelse(bmi >= 30, "Obese", NA))
            )
        )),
        cpd.cat = as.factor(ifelse(
            cpd > 0 & cpd < 20, "<20", ifelse(
                cpd >= 20 & cpd < 30,
                "20-29",
                ifelse(cpd >= 30 &
                           cpd < 40, "30-39", ifelse(cpd >= 40 & cpd < 99, "40+", NA))
            )
        )),
        smkyears.cat = as.factor(ifelse(
            smkyears > 0 &
                smkyears < 30,
            "<30",
            ifelse(
                smkyears >= 30 & smkyears < 40,
                "30-39",
                ifelse(
                    smkyears >= 40 &
                        smkyears < 50,
                    "40-49",
                    ifelse(smkyears >= 50 & smkyears < 99, "50+", NA)
                )
            )
        ))
    )


# --------------------------- run code to this line for data setup ----------------------------- #


#### Descriptive stats ####

# Descriptive characteristics - table 1
nrow(all.subj.pos)  # number of unique individuals
quantile(all.subj.pos$prescr.1yrisk.T0, probs = c(0.25, 0.5, 0.75))  # median IQR of pre-screening risk
CrossTable((
    data_screen_abn_pos %>% group_by(pid) %>% summarise(times.in.analysis =
                                                            n())
)$times.in.analysis) # times included in next-screen analysis
CrossTable(all.subj.pos$female, missing.include = T)
CrossTable(all.subj.pos$age.cat, missing.include = T)
CrossTable(all.subj.pos$race, missing.include = T)   # 0 white, 1 black, 2 hispanic, 3 other
CrossTable(all.subj.pos$edu6, missing.include = T)   # see codebook
CrossTable(all.subj.pos$bmi.cat, missing.include = T)
CrossTable(all.subj.pos$fam.lung.trend, missing.include = T)  # none, 1, 2+
CrossTable(all.subj.pos$qtyears.cat, missing.include = T)
CrossTable(all.subj.pos$pkyears.cat, missing.include = T)
CrossTable(all.subj.pos$smkyears.cat, missing.include = T)
CrossTable(all.subj.pos$cpd.cat, missing.include = T)
CrossTable(all.subj.pos$emp, missing.include = T)

# Prevalence of abnormalities - table 2
lapply(data_screen_abn_pos[data_screen_abn_pos$interval == 2, abnlist.pos], function(x)
    prop.table(table(x))) # at T1 abn screen

# Prevalence and risks - table 3
# Get contrasts for diameter and confidence intervals.
summary(glm_screen_pos_abn_log)
lapply(coefficients(glm_screen_pos_abn_log)[1:6], function(x)
    x - coefficients(glm_screen_pos_abn_log)[2])
confint(glm_screen_pos_abn_log)
# Risks
tab3risk <-
    function(var, cat) {
        # this function prints the prevalence as a percentage, then median/IQR risk
        df <-
            subset(
                dplyr::filter(data_screen_abn_pos, interval == 2),
                select = c("post.risk.abn.pos", var)
            )
        df <- df[df[, 2] %in% cat, ]
        print(round(100 * (nrow(df) / nrow(
            dplyr::filter(data_screen_abn_pos, interval == 2)
        )), 3))
        100 * round(quantile(df[, 1], probs = c(0.25, 0.5, 0.75)), 5)
    }
tab3risk("diam.cat", "0")
tab3risk("diam.cat", "4-5")
tab3risk("diam.cat", "6-7")
tab3risk("diam.cat", "8-10")
tab3risk("diam.cat", "11-13")
tab3risk("diam.cat", "14+")
tab3risk("any.upper", 1)
tab3risk("any.growth", 1)
tab3risk("any.mixed", 1)
tab3risk("any.spiculation", 1)
# Prevalence/probs for any right mid or lingula
100 * nrow(dplyr::filter(
    data_screen_abn_pos,
    interval == 2 &
        (any.right.mid == 1 |
             any.lingula == 1)
)) / nrow(filter(data_screen_abn_pos, interval == 2))
100 * with(
    dplyr::filter(
        data_screen_abn_pos,
        interval == 2 &
            (any.right.mid == 1 |
                 any.lingula == 1)
    ),
    quantile(post.risk.abn.pos, probs = c(0.25, 0.5, 0.75))
)
# Prevalence/probs for any poorly defined margins or any unable to characterize margins
100 * nrow(dplyr::filter(
    data_screen_abn_pos,
    interval == 2 &
        (any.poor.def == 1 |
             any.margin.unab == 1)
)) / nrow(filter(data_screen_abn_pos, interval == 2))
100 * with(
    dplyr::filter(
        data_screen_abn_pos,
        interval == 2 &
            (any.poor.def == 1 |
                 any.margin.unab == 1)
    ),
    quantile(post.risk.abn.pos, probs = c(0.25, 0.5, 0.75))
)

#### Properties of NLST screening, abnormal screens ####
# Is pre-screening risk important?
glm_screen_pos_nopsr <-
    glm(case ~ 1, data = data_screen_abn_pos, family = binomial(link = 'log'))
summary(glm_screen_pos_nopsr)
summary(glm_screen_pos)
lrtest(glm_screen_pos_nopsr, glm_screen_pos)
# Does the interval matter? no (p=0.13)
glm_screen_pos_by_int <-
    glm(
        case ~ log1yrisk:as.factor(interval) - 1,
        data = data_screen_abn_pos,
        family = binomial(link = 'log')
    )
summary(glm_screen_pos_by_int)
1 - pchisq(
    glm_screen_pos$deviance - glm_screen_pos_by_int$deviance,
    length(glm_screen_pos_by_int$coefficients) - length(glm_screen_pos$coefficients)
)
with(data_screen_abn_pos, table(interval))
# Do previous screens matter? no (p=0.62)
glm_screen_pos_2levels <-
    glm(
        case ~ log1yrisk:screen_comb - 1,
        data = dplyr::filter(data_screen_abn_pos, interval == 2 &
                                 !is.na(screen_comb)),
        family = binomial(link = 'log')
    )
summary(glm_screen_pos_2levels)
confint(glm_screen_pos_2levels)
glm_screen_pos_1level <-
    glm(
        case ~ log1yrisk - 1,
        data = dplyr::filter(data_screen_abn_pos, interval == 2 &
                                 !is.na(screen_comb)),
        family = binomial(link = 'log')
    )
summary(glm_screen_pos_1level)
1 - pchisq(
    glm_screen_pos_1level$deviance - glm_screen_pos_2levels$deviance,
    length(glm_screen_pos_2levels$coefficients) - length(glm_screen_pos_1level$coefficients)
)
with(filter(data_screen_abn_pos, interval == 2), table(screen_comb))
# Residual effects of LCRAT variables examined below

#### Effects of abnormalities and nodules ####
# Backwards stepwise: selects any.lingula, any.mixed, any.other.att, any.right.mid, any.poor.def, any.upper, any.margin.unab, any.spiculation. diam.cat, any.growth
scr_pos_full <-
    glm(
        paste(
            "case ~ logit1yrisk -1 +",
            paste(abnlist.pos.int, collapse = "+"),
            sep = ""
        ),
        data = dplyr::filter(data_screen_abn_pos),
        family = binomial(link = 'logit')
    )
bsw_scr_pos <-
    step(
        scr_pos_full,
        direction = "backward",
        scope = list(lower = case ~ logit1yrisk - 1, upper = scr_pos_full)
    )
# Fit this model
mod_vars_bsw <-
    glm(
        case ~ logit1yrisk:any.lingula + logit1yrisk:any.mixed + logit1yrisk:any.right.mid + logit1yrisk:any.poor.def + logit1yrisk:any.upper +
            logit1yrisk:any.other.att + logit1yrisk:any.margin.unab + logit1yrisk:any.spiculation + logit1yrisk:diam_cat + logit1yrisk:any.growth -
            1,
        data = data_screen_abn_pos,
        family = binomial(link = 'logit')
    )
summary(mod_vars_bsw)
  # get contrast amounts for diameter variable (for stepwise model in table)
lapply(coefficients(mod_vars_bsw)[9:14], function(x)
    x - coefficients(mod_vars_bsw)[10])
    # Get LRT p-value for diam.cat: p<0.000001
mod_vars_bsw_nodiamcat <-
    glm(
        case ~ logit1yrisk + logit1yrisk:any.lingula + logit1yrisk:any.mixed + logit1yrisk:any.right.mid + logit1yrisk:any.poor.def + logit1yrisk:any.upper +
            logit1yrisk:any.other.att + logit1yrisk:any.margin.unab + logit1yrisk:any.spiculation + logit1yrisk:any.growth -
            1,
        data = data_screen_abn_pos,
        family = binomial(link = 'logit')
    )
lrtest(mod_vars_bsw_nodiamcat, mod_vars_bsw)
# Lasso with intermediate lambda - selects any.upper, any.spiculation, any.smooth, any.mixed, any.other.att, any.growth, 4 of 6 diam.cat categories 
# why isn't this omitting the intercept? 
set.seed(61116)
x = model.matrix(case ~ logit1yrisk + logit1yrisk:. - 1, data = data_screen_abn_pos[, c("case", "logit1yrisk", abnlist.pos)])
cv.lasso = cv.glmnet(x,
                     data_screen_abn_pos$case,
                     alpha = 1,
                     family = "binomial")
out <-
    glmnet(x,
           data_screen_abn_pos$case,
           alpha = 1,
           family = "binomial")
predict(out,
        type = "coefficients",
        s = (cv.lasso$lambda.min + cv.lasso$lambda.1se) / 2)
    # Fit this model - log scale
mod.vars.lasso <- glm(case ~ log1yrisk:any.upper + log1yrisk:any.spiculation + log1yrisk:any.smooth + log1yrisk:any.mixed +
                        log1yrisk:any.other.att + log1yrisk:any.growth + log1yrisk:diam.cat -1, data=data_screen_abn_pos, family=binomial(link='log'))
summary(mod.vars.lasso)
    # get contrast amounts for diameter variable (for lasso model in table)
lapply(coefficients(mod.vars.lasso)[7:12], function(x) x-coefficients(mod.vars.lasso)[8])
    # Get LRT p-value for diam.cat: p<0.00001
mod.vars.lasso.nodiamcat <- glm(case ~ log1yrisk + log1yrisk:any.upper + log1yrisk:any.spiculation + log1yrisk:any.smooth + log1yrisk:any.mixed +
                        log1yrisk:any.other.att + log1yrisk:any.growth -1, data=data_screen_abn_pos, family=binomial(link='log'))
lrtest(mod.vars.lasso.nodiamcat, mod.vars.lasso)
# Try fitting a final model:
  # 1) combine right-mid + lingula location as one effect. 2) combine poor-def and cannot-determine margins as one effect.
  # 3) fit a model with all variables selected by either approach. 4) remove any.smooth (soft tissue attenuation) where beta=0.01 and p=0.73
  # 5) remove any.other.att - it's hard to interpret; we don't know what it means
  # 6) try replacing diam.cat with diam.cat.GG to look at whether GG <20/>=20 mm is important (AIC is much worse)
  # 7) try adding number of nodules to confirm it doesn't contribute
  # 8) try adding any.new.nodule to confirm it doesn't contribute
  # 9) try adding an effect for new nodules 4-7mm in diameter to allow the effect of diameter to vary. doesn't contribute
  # 10) add an effect for black race to see how much it adds and whether it changes other parameters
glm.screen.pos.abn.dev <- glm(case ~ log1yrisk:diam.cat + log1yrisk:any.growth +     # replace diam.cat with diam.cat.GG to test AIC
                      log1yrisk:any.upper + log1yrisk:I(any.right.mid==1|any.lingula==1) +
                      log1yrisk:any.mixed + 
                     # log1yrisk:any.other.att +        # REMOVE - too difficult to interpret
                     # log1yrisk:any.smooth +       # REMOVE - beta=0.01, p=0.73
                     # log1yrisk:nodule.count +     # to confirm it does not contribute: p=0.29
                     # log1yrisk:any.new.nodule +   # to confirm it does not contribute: p=0.77
                     # log1yrisk:any.new.nodule.4.7.mm +   # to confirm it does not contribute: p=0.86
                     # log1yrisk:I(race==1) +  # p=0.006, change in exponent -0.13. other coefficients don't change.
                       log1yrisk:any.spiculation + log1yrisk:I(any.poor.def==1|any.margin.unab==1) -1, 
                      data=data_screen_abn_pos, family=binomial(link='log'))
summary(glm.screen.pos.abn.dev)
# Get LRT p-value for diam.cat
glm.screen.pos.abn.dev.nodiamcat <- glm(case ~ log1yrisk + log1yrisk:any.growth +
                                log1yrisk:any.upper + log1yrisk:I(any.right.mid==1|any.lingula==1) +
                                log1yrisk:I(any.new.nodule==T & diam.cat %in% c("4-5","6-7")) +
                                log1yrisk:any.mixed + log1yrisk:any.spiculation + log1yrisk:I(any.poor.def==1|any.margin.unab==1) -1, 
                              data=data_screen_abn_pos, family=binomial(link='log'))
lrtest(glm.screen.pos.abn.dev.nodiamcat, glm.screen.pos.abn.dev)
# Compare this model to one that ignores pre-screening risk
glm.screen.pos.abn.nopsr <-  glm(case ~ diam.cat + any.growth + any.upper + 
                      I(any.right.mid==1|any.lingula==1) + any.mixed + any.spiculation + 
                      I(any.new.nodule==T & diam.cat %in% c("4-5","6-7")) +
                      I(any.poor.def==1|any.margin.unab==1) -1, data=data_screen_abn_pos, family=binomial(link='log'))
summary(glm.screen.pos.abn.nopsr)
# Check for residual effects of LCRAT variables.
titles <- c("var","null model # param", "extended model # param", "LRT p-value", "check same # obs")
mat.out.ns <- matrix(rep(NA),nrow=length(varlist),ncol=length(titles))
mod.without <- glm(case ~ log1yrisk:diam.cat + log1yrisk:any.growth + log1yrisk:any.upper + log1yrisk:I(any.right.mid==1|any.lingula==1) + log1yrisk:any.mixed + log1yrisk:any.spiculation + log1yrisk:I(any.poor.def==1|any.margin.unab==1) + log1yrisk:I(any.new.nodule==T & diam.cat %in% c("4-5","6-7")) -1, data=data_screen_abn_pos, family=binomial(link='log'))
for (x in seq_along(varlist)) {
  mod.with <- glm(substitute(case ~ log1yrisk:diam.cat + log1yrisk:any.growth + log1yrisk:any.upper + log1yrisk:I(any.right.mid==1|any.lingula==1) + log1yrisk:any.mixed + log1yrisk:any.spiculation + log1yrisk:I(any.poor.def==1|any.margin.unab==1) + log1yrisk:i -1, list(i=as.name(varlist[x]))), data=data_screen_abn_pos, family=binomial(link='log'))
  print(summary(mod.with))
  mat.out.ns[x,] <- c(varlist[x], length(mod.without$coefficients), sum(!is.na(mod.with$coefficients)), 1-pchisq(mod.without$deviance-mod.with$deviance, df=sum(!is.na(mod.with$coefficients))-length(mod.without$coefficients)), I(length(mod.without$residuals)==length(mod.with$residuals)))
}
rbind(titles, mat.out.ns)  # These are the LRT p-values, but in the manuscript I used the Wald p-values from the models.
  # The issue of how to model diameter is explored in more detail in analysis_nlst_v10.R. 
      # I tried various approaches and transformation and found that a categorical approach is best (lowest AIC).
  # I tried replacing any.growth with growth.3l, which separates those at T0-T1 where growth can't really be assessed.
      # But the coefficient was essentially exactly the same for the no group and the not-evaluable group, so I switched back to any.growth.


#### Additional analyses ####

### Comparison with GEE - this impacts the SEs negligibly ###
# Screen-false-positives model
summary(geeglm(case ~ log1yrisk:diam.cat + log1yrisk:any.growth + log1yrisk:any.upper + 
                log1yrisk:I(any.right.mid==1|any.lingula==1) + log1yrisk:any.mixed + 
                log1yrisk:any.spiculation + log1yrisk:I(any.poor.def==1|any.margin.unab==1) -1, 
                id=pid, data=data_screen_abn_pos, family=binomial(link='log'), corstr="exchangeable", waves=interval))
summary(glm.screen.pos.abn.log)

### Calibration and validation analyses ###
# 10-fold cross-validated calibration
set.seed(61116)
data_screen_abn_pos$randgrp <- base::sample(1:10, nrow(data_screen_abn_pos), replace=T)
data_screen_abn_pos$cvpred <- NA
for (i in 1:10) {
  fit <- glm(formula = case ~ log1yrisk:diam.cat + log1yrisk:any.growth + log1yrisk:any.upper + 
               log1yrisk:I(any.right.mid==1|any.lingula==1) + log1yrisk:any.mixed + 
               log1yrisk:any.spiculation + log1yrisk:I(any.poor.def==1|any.margin.unab==1) -1, 
             family = binomial(link = "log"), data = dplyr::filter(data_screen_abn_pos, randgrp!=i)) 
  data_screen_abn_pos[data_screen_abn_pos$randgrp==i,]$cvpred <- predict(fit, newdata=data_screen_abn_pos[data_screen_abn_pos$randgrp==i,], type="response")
}
data_screen_abn_pos <- mutate(data_screen_abn_pos, cvpred.ntile = ntile(cvpred, 5))
data_screen_abn_pos %>% group_by(cvpred.ntile) %>% summarise(pred.cases= sum(cvpred), obs.cases = sum(case))
c(sum(data_screen_abn_pos$cvpred), sum(data_screen_abn_pos$case))
poisson.test(round(sum(data_screen_abn_pos$cvpred),0), sum(data_screen_abn_pos$case), alternative="two.sided") # p-value, requires rounding

### Calculate AUCs ###
## Optimism-corrected AUCs - have to use logistic models for this, and have to actually add the interaction terms to the dataset.
  # http://thestatsgeek.com/2014/10/04/adjusting-for-optimismoverfitting-in-measures-of-predictive-ability-using-bootstrapping/
# install.packages('rms')
library(rms)
data_screen_abn_pos <- mutate(data_screen_abn_pos, logit1yriskdiamcat0 = logit1yrisk*I(diam.cat=="0"), 
                              logit1yriskdiamcat4_5 = logit1yrisk*I(diam.cat=="4-5"),
                              logit1yriskdiamcat6_7 = logit1yrisk*I(diam.cat=="6-7"),
                              logit1yriskdiamcat8_10 = logit1yrisk*I(diam.cat=="8-10"),
                              logit1yriskdiamcat11_13 = logit1yrisk*I(diam.cat=="11-13"),
                              logit1yriskdiamcat14 = logit1yrisk*I(diam.cat=="14+"),
                              logit1yriskanygrowth = logit1yrisk*any.growth,
                              logit1yriskanyupper = logit1yrisk*any.upper, 
                              logit1yriskanyrightmidanylingula = logit1yrisk*I(any.right.mid==1 | any.lingula==1),
                              logit1yriskanymixed = logit1yrisk*any.mixed, logit1yriskanyspiculation = logit1yrisk*any.spiculation,
                              logit1yriskanypoordefanymarginunab = logit1yrisk*I(any.poor.def==1 | any.margin.unab==1))
mod.ns.pos <- lrm(case ~ logit1yriskdiamcat0 + logit1yriskdiamcat4_5 + logit1yriskdiamcat6_7 +
                    logit1yriskdiamcat8_10 + logit1yriskdiamcat11_13 + logit1yriskdiamcat14 + logit1yriskanygrowth +
                    logit1yriskanyupper + logit1yriskanyrightmidanylingula + logit1yriskanymixed + logit1yriskanyspiculation +
                    logit1yriskanypoordefanymarginunab -1, x=T, y=T, data=data_screen_abn_pos)
set.seed(61116)
validate(mod.ns.pos, B=1000)
c(0.5*(0.5938+1), 0.5*(0.5767+1)) # AUC = 0.5(Dxy+1). 
  # the above gives naive (index.orig) and optimism-corrected (index.corrected) AUCs. 0.79 is OC-AUC



#### Figure 1 - risk density #### 

prq <- with(dplyr::filter(data_screen_abn_pos, interval==2), quantile(post.risk.abn.pos, probs=c(0.5, 0.75, 0.95)))
prq_ar_end <- c(prq[1]-0.003, prq[2]-0.003, prq[3]+0.003)
prq

png(file=here("abnormals_density.png"),width=1200,height=850)
ggplot(data=dplyr::filter(data_screen_abn_pos, interval==2)) + theme_bw() +
  theme(axis.title = element_text(size=30), axis.text = element_text(size=28), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_y_continuous(name="Density", breaks=NULL) +
  scale_x_continuous(name="\nPredicted lung cancer risk", limits=c(0, 0.08), 
                     breaks=c(0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08), 
                     labels=c("0%","1%","2%","3%","4%","5%","6%","7%","8%")) +
  geom_density(aes(x=prescr.1yrisk, y=..density..), fill="gray60", colour="gray60", alpha=0.3) +
  geom_density(aes(x=post.risk.abn.pos, y=..density..), fill="midnightblue", colour="midnightblue", alpha=0.5) + 
  geom_segment(aes(x=prq[1], xend=prq[1], y=0, yend=135), color="midnightblue") +
  geom_segment(aes(x=prq[1], xend=prq_ar_end[1], y=135, yend=135), arrow=arrow(), colour="midnightblue") +
  geom_segment(aes(x=prq[2], xend=prq[2], y=0, yend=95), color="midnightblue") +
  geom_segment(aes(x=prq[2], xend=prq_ar_end[2], y=95, yend=95), arrow=arrow(), colour="midnightblue") +
  geom_segment(aes(x=prq[3], y=0, xend=prq[3], yend=55), color="midnightblue") +
  geom_segment(aes(x=prq[3], xend=prq_ar_end[3], y=55, yend=55), arrow=arrow(), colour="midnightblue") +
  annotate(geom="text", x=0.033, y=178, size=13, colour="gray60", label="1-year pre-screening risk before abnormal CT") +
  annotate(geom="text", x=0.037, y=30, size=13, colour="midnightblue", label="Next-screen risk\nafter abnormal CT") +
  annotate(geom="text", x=0.018, y=150, size=9, label="50% of abnormal CTs have \n next-screen risk below 0.90%", color="midnightblue") +
  annotate(geom="text", x=0.030, y=110, size=9, label="75% of abnormal CTs have \n next-screen risk below 2.1%", color="midnightblue") +
  annotate(geom="text", x=0.069, y=71, size=9, label="5% of abnormal CTs have \n next-screen risk above 7.1%", color="midnightblue") +
  NULL
dev.off()
# Look at percentiles, look at those in a reasonable risk range for continued annual screening
with(dplyr::filter(data_screen_abn_pos, interval==2), quantile(post.risk.abn.pos, probs=seq(0,1,0.1)))
with(dplyr::filter(data_screen_abn_pos, interval==2), CrossTable(I(post.risk.abn.pos>0.002), I(post.risk.abn.pos<0.03)))

# Look at the 5% with highest risk (>7.1%)
View(dplyr::filter(data_screen_abn_pos, interval==2 & post.risk.abn.pos>0.071))
with(dplyr::filter(data_screen_abn_pos, interval==2 & post.risk.abn.pos>0.071), CrossTable(diam.cat))
with(dplyr::filter(data_screen_abn_pos, interval==2 & post.risk.abn.pos>0.071), CrossTable(any.upper))
with(dplyr::filter(data_screen_abn_pos, interval==2 & post.risk.abn.pos>0.071), CrossTable(I(any.right.mid == 1 | any.lingula == 1)))
with(dplyr::filter(data_screen_abn_pos, interval==2 & post.risk.abn.pos>0.071), CrossTable(any.mixed))
with(dplyr::filter(data_screen_abn_pos, interval==2 & post.risk.abn.pos>0.071), CrossTable(any.spiculation))
with(dplyr::filter(data_screen_abn_pos, interval==2 & post.risk.abn.pos>0.071), CrossTable(I(any.poor.def == 1 | any.margin.unab == 1)))
with(dplyr::filter(data_screen_abn_pos, interval==2 & post.risk.abn.pos>0.071), CrossTable(any.growth))

# Look at the 8% with a risk <0.20% (e.g. potential extended interval)
View(dplyr::filter(data_screen_abn_pos, interval==2 & post.risk.abn.pos<0.002))
      # Among these 8% look at 5-year pre-screening risk (1.9% threshold is proposed by LCRAT)
with(dplyr::filter(data_screen_abn_pos, interval==2 & post.risk.abn.pos<0.002), quantile(I(5*prescr.1yrisk), probs=seq(0,1,.1)))
with(dplyr::filter(data_screen_abn_pos, interval==2), quantile(I(5*prescr.1yrisk), probs=seq(0,1,.1)))
with(dplyr::filter(data_screen_abn_pos, interval==2 & post.risk.abn.pos<0.002), CrossTable(diam.cat))

