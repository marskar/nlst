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

# Read in the Kovalchik prediction function
source(here("R/kovalchik.R"))

devtools::install_github('marskar/coxph_risk')
devtools::install_github('marskar/lcmodels')

# Function to create "not in" operator
'%!in%' <- function(x, y) {
    !('%in%'(x, y))
}


# Merge in Wes' T0 data
nlst_emp <- read_csv(here('data/T0_data.csv'))
data_screen_abn_neg <- readRDS(here('data/data_screen_abn_neg.rds'))
data_screen_abn_neg_emp <-
    merge(data_screen_abn_neg, nlst_emp, by = "pid")
data_screen_abn_negt0_emp <-
    filter(data_screen_abn_neg_emp, interval == 1)
data_screen_abn_negt1_emp <-
    filter(data_screen_abn_neg_emp, interval == 2)
data_screen_abn_negt0_emp$weight <-
    ifelse(data_screen_abn_negt0_emp$case == 1, 1, 20)

abn <- readRDS(here('data/abn_lrads_merged.rds'))
data_interval_abn <- readRDS(here('data/data_interval_abn.rds'))
data_interval_abn_emp <-
    merge(data_interval_abn, nlst_emp, by = "pid")
data_interval_abn_t0_emp <-
    filter(data_interval_abn_emp, interval == 1)
data_interval_abn_t0_emp$weight <-
    ifelse(data_interval_abn_t0_emp$case == 1, 1, 20)


# Run the main models

# # With specific CT findings
# glm_int_abn <-
#     glm(
#         case ~ log1yrisk + log1yrisk:adenop.consol - 1,
#         data = data_interval_abn,
#         family = binomial(link = 'log')
#     )
# data_interval_abn$post_risk_abn <- fitted.values(glm_int_abn)
#
# # data_screen_abn_emph_neg_drop_na <-
# #     data_screen_abn_emph_neg %>% drop_na(p_emph)
# # data_interval_abn_emp <-
# #     data_interval_abn_emp %>% tidyr::drop_na(p_emph)
#
#
# glm_screen_neg_abn <-
#     glm(
#         case ~ log1yrisk + log1yrisk:consolidation + log1yrisk:emphysema - 1,
#         data = data_screen_abn_emph_neg_drop_na,
#         family = binomial(link = 'log'),
#         na.action = na.exclude
#     )
#
# data_screen_abn_emph_neg$post_risk_abn <-
#     fitted.values(glm_screen_neg_abn)
#
# glm_screen_neg_abn <-
#     glm(
#         case ~ log1yrisk - 1,
#         data = data_screen_abn_emph_neg_drop_na,
#         family = binomial(link = 'log'),
#         na.action = na.exclude
#     )
#
# predictions <- predict(glm_screen_neg_abn, type = "response")
# pROC::auc(data_screen_abn_emph_neg_drop_na$case, predictions)
#
# glm_screen_neg_abn_emph_logit <-
#     glm(
#         case ~ log1yrisk + I(logit(p_emph)) - 1,
#         data = data_screen_abn_emph_neg_drop_na,
#         family = binomial(link = 'log'),
#         na.action = na.exclude
#     )
#
# glm_screen_neg_abn_emph <-
#     glm(
#         case ~ log1yrisk + p_emph - 1,
#         data = data_screen_abn_emph_neg_drop_na,
#         family = binomial(link = 'log'),
#         na.action = na.exclude
#     )
# predictions_emph <-
#     predict(glm_screen_neg_abn_emph, type = "response")
# predictions_emph_logit <-
#     predict(glm_screen_neg_abn_emph_logit, type = "response")
# predictions_physician_emph_logit <-
#     predict(glm_screen_neg_abn, type = "response")
# pROC::auc(data_screen_abn_emph_neg_drop_na$case, predictions_emph)
# pROC::auc(data_screen_abn_emph_neg_drop_na$case,
#           predictions_emph_logit)
# pROC::auc(data_screen_abn_emph_neg_drop_na$case,
#           predictions_physician_emph_logit)
#
# summary(glm_screen_neg_abn)
# summary(glm_screen_neg_abn_emph)
# summary(glm_screen_neg_abn_emph_logit)
# summary(glm_screen_neg_abn)


#all data
# glm_screen_neg_abn <-
#     glm(
#         case ~ log1yrisk + log1yrisk:consolidation + log1yrisk:emphysema - 1,
#         data = data_screen_abn_neg,
#         family = binomial(link = 'log'),
#         na.action = na.exclude
#     )
# 
# glm_screen_neg_abn_sum <- summary(glm_screen_neg_abn)
# coef(glm_screen_neg_abn_sum)
# 
# data_screen_abn_neg <- data_screen_abn_neg %>%
#     mutate(post_risk_abn = fitted.values(glm_screen_neg_abn))
# summary(glm_screen_neg_abn)
# 
# 
# glm_interval <-
#     glm(case ~ log1yrisk - 1,
#         data = data_interval_abn,
#         family = binomial(link = 'log'))

glm_screen_neg <-
    glm(case ~ log1yrisk - 1,
        data = data_screen_abn_neg_emp,
        family = binomial(link = 'log'))
summary(glm_screen_neg)

glm_screen_neg_pemph <-
    glm(case ~ log1yrisk + I(logit(p_emph)) - 1,
        data = data_screen_abn_neg_emp,
        family = binomial(link = 'log'),
        na.action = na.exclude)
summary(glm_screen_neg_pemph)

glm_screen_neg_pemph_cons <-
    glm(
        case ~ log1yrisk + I(logit(p_emph)) + log1yrisk:consolidation - 1,
        data = data_screen_abn_neg_emp,
        family = binomial(link = 'log'),
        na.action = na.exclude
    )
summary(glm_screen_neg_pemph_cons)

glm_screen_neg_pemph_cons_emp <-
    glm(
        case ~ log1yrisk + I(logit(p_emph)) +log1yrisk:consolidation + log1yrisk:emphysema - 1,
        data = data_screen_abn_neg_emp,
        family = binomial(link = 'log'),
        na.action = na.exclude
    )
summary(glm_screen_neg_pemph_cons_emp)

glm_screen_neg_abn_emp <-
    glm(
        case ~ log1yrisk + log1yrisk:consolidation + log1yrisk:emphysema - 1,
        data = data_screen_abn_neg_emp,
        family = binomial(link = 'log'),
        na.action = na.exclude
    )
summary(glm_screen_neg_abn_emp)
# AIC: 862.8, residual deviance =856.8

# glm_screen_neg_abn_emp1 <- glm(case ~ log1yrisk + log1yrisk:consolidation + log1yrisk:score -1, data=data_screen_abn_neg_emp, family=binomial(link='log'), na.action=na.exclude)
# data_screen_abn_neg_emp$post_risk_abn_emp1 <- fitted.values(glm_screen_neg_abn_emp1)
#AIC: 862.7, residual deviance =856.7

# glm_screen_neg_abn_emp2 <- glm(case ~ log1yrisk + log1yrisk:consolidation + log1yrisk:score2 -1, data=data_screen_abn_neg_emp, family=binomial(link='log'), na.action=na.exclude)
# data_screen_abn_neg_emp$post_risk_abn_emp2 <- fitted.values(glm_screen_neg_abn_emp2)
#AIC: 862.8, residual deviance =856.8

### restricted to screen detected cancers at T1:
glm_screen_negt0_abn_emp <-
    glm(
        case ~ log1yrisk + log1yrisk:consolidation + log1yrisk:emphysema - 1,
        data = data_screen_abn_negt0_emp,
        family = binomial(link = 'log'),
        na.action = na.exclude
    )
data_screen_abn_negt0_emp$post_risk_abn <-
    fitted.values(glm_screen_negt0_abn_emp)
glm_screen_negt0_abn_emp_sum <- summary(glm_screen_negt0_abn_emp)
coef(glm_screen_negt0_abn_emp_sum)
summary(glm_screen_negt0_abn_emp)
# AIC: 761.9, residual deviance =755.9

# Interval cancer model
# T0-T1 - change interval for T1-T2, post-T2
with(filter(data_screen_abn_negt0_emp, interval == 1),
     roc(case, post_risk_abn, ci = T))


# glm.screen.negT0.abn.emp1 <- glm(case ~ log1yrisk + log1yrisk:consolidation + log1yrisk:score -1, data=data_screen_abn_negt0_emp, family=binomial(link='log'), na.action=na.exclude)
# data_screen_abn_negt0_emp$post_risk_abn_emp1 <- fitted.values(glm.screen.negT0.abn.emp1)
# glm.screen.negT0.abn.emp1.sum<- summary(glm.screen.negT0.abn.emp1)
# coef(glm.screen.negT0.abn.emp1.sum)
# T0-T1 - change interval for T1-T2, post-T2
# with(filter(data_screen_abn_negt0_emp, interval==1), roc(case, post_risk_abn.emp1, ci=T, plot=T))

# glm.screen.negT0.abn.emp2 <- glm(case ~ log1yrisk + log1yrisk:consolidation + log1yrisk:score2 -1, data=data_screen_abn_negt0_emp, family=binomial(link='log'), na.action=na.exclude)
# data_screen_abn_negt0_emp$post_risk_abn.emp2 <- fitted.values(glm.screen.negT0.abn.emp2)
# glm.screen.negT0.abn.emp2.sum<- summary(glm.screen.negT0.abn.emp2)
# coef(glm.screen.negT0.abn.emp2.sum)
# T0-T1 - change interval for T1-T2, post-T2
# with(filter(data_screen_abn_negt0_emp, interval==1), roc(case, post_risk_abn.emp2, ci=T, plot=T))


### restrcited to interval cancers T0-T1
glm_int_abn_emp <- glm(case ~ log1yrisk + log1yrisk:adenop.consol -1, data=data_interval_abn_t0_emp, family=binomial(link='log'))
data_interval_abn_t0_emp$post_risk_abn <- fitted.values(glm_int_abn_emp)
glm_int_abn_emp_sum<- summary(glm_int_abn_emp)
coef(glm_int_abn_emp_sum)

with(filter(data_interval_abn_t0_emp, interval==1), roc(case, post_risk_abn, ci=T, plot=T))  # T0-T1 - change interval for T1-T2, post-T2

glm.int.abn.emp0 <- glm(case ~ log1yrisk + log1yrisk:adenop.consol + log1yrisk:emphysema -1, data=data_interval_abn_t0_emp, family=binomial(link='log'))
data_interval_abn_t0_emp$post_risk_abn0 <- fitted.values(glm.int.abn.emp0)
glm.int.abn.emp0.sum<- summary(glm.int.abn.emp0)
coef(glm.int.abn.emp0.sum)

with(filter(data_interval_abn_t0_emp, interval==1), roc(case, post_risk_abn0, ci=T, plot=T))  # T0-T1 - change interval for T1-T2, post-T2

# glm.int.abn.emp1 <- glm(case ~ log1yrisk + log1yrisk:adenop.consol + log1yrisk:score -1, data=data_interval_abn_t0_emp, family=binomial(link='log'))
data_interval_abn_t0_emp$post_risk_abn1 <- fitted.values(glm.int.abn.emp1)
glm.int.abn.emp1.sum<- summary(glm.int.abn.emp1)
coef(glm.int.abn.emp1.sum)

with(filter(data_interval_abn_t0_emp, interval==1), roc(case, post_risk_abn1, ci=T, plot=T))  # T0-T1 - change interval for T1-T2, post-T2

# glm.int.abn.emp2 <- glm(case ~ log1yrisk + log1yrisk:adenop.consol + log1yrisk:score2 -1, data=data_interval_abn_t0_emp, family=binomial(link='log'))
data_interval_abn_t0_emp$post_risk_abn2 <- fitted.values(glm.int.abn.emp2)
glm.int.abn.emp2.sum<- summary(glm.int.abn.emp2)
coef(glm.int.abn.emp2.sum)

# T0-T1 - change interval for T1-T2, post-T2
with(filter(data_interval_abn_t0_emp, interval==1), roc(case, post_risk_abn2, ci=T, plot=T)) 






par(mfrow=c(2,3))

data_screen_abn_negt0_emp$score_min<-min(data_screen_abn_negt0_emp$score[data_screen_abn_negt0_emp$score>0])
data_screen_abn_negt0_emp$score_max<-max(data_screen_abn_negt0_emp$score[data_screen_abn_negt0_emp$score<1])
data_screen_abn_negt0_emp$score2_min<-min(data_screen_abn_negt0_emp$score2[data_screen_abn_negt0_emp$score2>0])
data_screen_abn_negt0_emp$score2_max<-max(data_screen_abn_negt0_emp$score2[data_screen_abn_negt0_emp$score2<1])



hist(data_screen_abn_negt0_emp$score,  main="original score")
hist(log(data_screen_abn_negt0_emp$score),  main="log(original score)")
hist(logit(data_screen_abn_negt0_emp$score), main="logit(original score)")
hist(data_screen_abn_negt0_emp$score2,  main="no-noise score")
hist(log(data_screen_abn_negt0_emp$score2),  main="log(no-noise score)")
hist(logit(data_screen_abn_negt0_emp$score2),  main="logit(no-noise score)")

data_screen_abn_negt0_emp$scoreb<-ifelse(data_screen_abn_negt0_emp$score==0,data_screen_abn_negt0_emp$score_min,ifelse(data_screen_abn_negt0_emp$score==1,data_screen_abn_negt0_emp$score_max,data_screen_abn_negt0_emp$score))
data_screen_abn_negt0_emp$score2b<-ifelse(data_screen_abn_negt0_emp$score2==0,data_screen_abn_negt0_emp$score2_min,ifelse(data_screen_abn_negt0_emp$score2==1,data_screen_abn_negt0_emp$score2_max,data_screen_abn_negt0_emp$score2))

data_screen_abn_negt0_emp$logscore<-log(data_screen_abn_negt0_emp$scoreb)
data_screen_abn_negt0_emp$logitscore<-logit(data_screen_abn_negt0_emp$scoreb)
data_screen_abn_negt0_emp$logscore2<-log(data_screen_abn_negt0_emp$score2b)
data_screen_abn_negt0_emp$logitscore2<-logit(data_screen_abn_negt0_emp$score2b)


### restricted to screen detected cancers at T1:
glm.screen.negT0.abn.logitemp1 <- glm(case ~ log1yrisk + log1yrisk:consolidation + log1yrisk:logitscore -1, data=data_screen_abn_negt0_emp, family=binomial(link='log'), na.action=na.exclude)
data_screen_abn_negt0_emp$post_risk_abn.logitemp1 <- fitted.values(glm.screen.negT0.abn.logitemp1)
glm.screen.negT0.abn.logitemp1.sum<- summary(glm.screen.negT0.abn.logitemp1)
coef(glm.screen.negT0.abn.logitemp1.sum)

with(filter(data_screen_abn_negt0_emp, interval==1), roc(case, post_risk_abn.logitemp1, ci=T, plot=T))  # T0-T1 - change interval for T1-T2, post-T2


glm.screen.negT0.abn.logitemp2 <- glm(case ~ log1yrisk + log1yrisk:consolidation + log1yrisk:logitscore2 -1, data=data_screen_abn_negt0_emp, family=binomial(link='log'), na.action=na.exclude)
data_screen_abn_negt0_emp$post_risk_abn.logitemp2 <- fitted.values(glm.screen.negT0.abn.logitemp2)
glm.screen.negT0.abn.logitemp2.sum<- summary(glm.screen.negT0.abn.logitemp2)
coef(glm.screen.negT0.abn.logitemp2.sum)

with(filter(data_screen_abn_negt0_emp, interval==1), roc(case, post_risk_abn.logitemp2, ci=T, plot=T))  # T0-T1 - change interval for T1-T2, post-T2

glm.screen.negT0.abn.logemp1 <- glm(case ~ log1yrisk + log1yrisk:consolidation + log1yrisk:logscore -1, data=data_screen_abn_negt0_emp, family=binomial(link='log'), na.action=na.exclude)
data_screen_abn_negt0_emp$post_risk_abn.logemp1 <- fitted.values(glm.screen.negT0.abn.logemp1)
glm.screen.negT0.abn.logemp1.sum<- summary(glm.screen.negT0.abn.logemp1)
coef(glm.screen.negT0.abn.logemp1.sum)

glm.screen.negT0.abn.logemp2 <- glm(case ~ log1yrisk + log1yrisk:consolidation + log1yrisk:logscore2 -1, data=data_screen_abn_negt0_emp, family=binomial(link='log'), na.action=na.exclude)
data_screen_abn_negt0_emp$post_risk_abn.logemp2 <- fitted.values(glm.screen.negT0.abn.logemp2)
glm.screen.negT0.abn.logemp2.sum<- summary(glm.screen.negT0.abn.logemp2)
coef(glm.screen.negT0.abn.logemp2.sum)

### interval cancers
data_interval_abn_t0_emp$score_min<-min(data_interval_abn_t0_emp$score[data_interval_abn_t0_emp$score>0])
data_interval_abn_t0_emp$score_max<-max(data_interval_abn_t0_emp$score[data_interval_abn_t0_emp$score<1])
data_interval_abn_t0_emp$score2_min<-min(data_interval_abn_t0_emp$score2[data_interval_abn_t0_emp$score2>0])
data_interval_abn_t0_emp$score2_max<-max(data_interval_abn_t0_emp$score2[data_interval_abn_t0_emp$score2<1])


data_interval_abn_t0_emp$scoreb<-ifelse(data_interval_abn_t0_emp$score==0,data_interval_abn_t0_emp$score_min,ifelse(data_interval_abn_t0_emp$score==1,data_interval_abn_t0_emp$score_max,data_interval_abn_t0_emp$score))
data_interval_abn_t0_emp$score2b<-ifelse(data_interval_abn_t0_emp$score2==0,data_interval_abn_t0_emp$score2_min,ifelse(data_interval_abn_t0_emp$score2==1,data_interval_abn_t0_emp$score2_max,data_interval_abn_t0_emp$score2))

data_interval_abn_t0_emp$logscore<-log(data_interval_abn_t0_emp$scoreb)
data_interval_abn_t0_emp$logitscore<-logit(data_interval_abn_t0_emp$scoreb)
data_interval_abn_t0_emp$logscore2<-log(data_interval_abn_t0_emp$score2b)
data_interval_abn_t0_emp$logitscore2<-logit(data_interval_abn_t0_emp$score2b)

### restricted to interval cancers T0- T1:
glm.int.abn.logitemp1 <- glm(case ~ log1yrisk + log1yrisk:consolidation + log1yrisk:logitscore -1, data=data_interval_abn_t0_emp, family=binomial(link='log'), na.action=na.exclude)
data_interval_abn_t0_emp$post_risk_abn.logitemp1 <- fitted.values(glm.int.abn.logitemp1)
glm.int.abn.logitemp1.sum<- summary(glm.int.abn.logitemp1)
coef(glm.int.abn.logitemp1.sum)

with(filter(data_interval_abn_t0_emp, interval==1), roc(case, post_risk_abn.logitemp1, ci=T, plot=T))  # T0-T1 - change interval for T1-T2, post-T2

glm.int.abn.logitemp2 <- glm(case ~ log1yrisk + log1yrisk:consolidation + log1yrisk:logitscore2 -1, data=data_interval_abn_t0_emp, family=binomial(link='log'), na.action=na.exclude)
data_interval_abn_t0_emp$post_risk_abn.logitemp2 <- fitted.values(glm.int.abn.logitemp2)
glm.int.abn.logitemp2.sum<- summary(glm.int.abn.logitemp1)
coef(glm.int.abn.logitemp2.sum)

with(filter(data_interval_abn_t0_emp, interval==1), roc(case, post_risk_abn.logitemp2, ci=T, plot=T))  # T0-T1 - change interval for T1-T2, post-T2





# Make dataset of unique individuals for descriptive table of screen-negatives
all.subj.neg <- filter(nlst.CT, pid %in% data_interval_abn$pid | pid %in% data_screen_abn_emph.neg$pid)
all.subj.neg <- mutate(all.subj.neg, age.cat=as.factor(ifelse(age>=55 & age<60, "55-59", ifelse(age>=60 & age<65, "60-64",
                    ifelse(age>=65 & age<70, "65-69", ifelse(age>=70 & age<75, "70-74", NA))))),
                  qtyears.cat=as.factor(ifelse(qtyears==0, "Current smoker", ifelse(qtyears>0 & qtyears<=5, "1-5",
                    ifelse(qtyears>5 & qtyears<=10, "6-10", ifelse(qtyears>10 & qtyears<99, "11 or more", NA))))),
                  bmi.cat=as.factor(ifelse(bmi>0 & bmi<18.5, "Underweight", ifelse(bmi>=18.5 & bmi<25, "Normal",
                    ifelse(bmi>=25 & bmi<30, "Overweight", ifelse(bmi>=30, "Obese", NA))))),
                  cpd.cat=as.factor(ifelse(cpd>0 & cpd<20, "<20", ifelse(cpd>=20 & cpd<30, "20-29",
                    ifelse(cpd>=30 & cpd<40, "30-39", ifelse(cpd>=40 & cpd<99, "40+", NA))))),
                  smkyears.cat=as.factor(ifelse(smkyears>0 & smkyears<30, "<30", ifelse(smkyears>=30 & smkyears<40, "30-39",
                    ifelse(smkyears>=40 & smkyears<50, "40-49", ifelse(smkyears>=50 & smkyears<99, "50+", NA))))))

# --------------------------- run code to this line for data setup ----------------------------- #


#### Descriptive stats ####

# Descriptive table 1 for analysis after NEGATIVE screen
nrow(all.subj.neg)  # number of unique individuals
length(unique(filter(nlst.CT, T0posneg==0 | T1posneg==0 | T2posneg==0)$pid)) # confirm - this is the same # with at least one negative screen
length(unique(data.interval$pid)) # number of unique individuals in interval cancers analysis
CrossTable((data.interval %>% group_by(pid) %>% summarise(times.in.interv.analysis=n()))$times.in.interv.analysis) # times included in interval ca analysis
length(unique(data_screen_abn_neg$pid)) # number of unique individuals in next-screen analysis
CrossTable((data_screen_abn_neg %>% group_by(pid) %>% summarise(times.in.screen.analysis=n()))$times.in.screen.analysis) # times included in next-screen analysis
CrossTable(all.subj.neg$female, missing.include=T)
CrossTable(all.subj.neg$age.cat, missing.include=T)
CrossTable(all.subj.neg$race, missing.include=T)   # 0 white, 1 black, 2 hispanic, 3 other
CrossTable(all.subj.neg$edu6, missing.include=T)   # see codebook
CrossTable(all.subj.neg$bmi.cat, missing.include=T)
CrossTable(all.subj.neg$fam.lung.trend, missing.include=T)  # none, 1, 2+
CrossTable(all.subj.neg$qtyears.cat, missing.include=T)
CrossTable(all.subj.neg$pkyears.cat, missing.include=T)
CrossTable(all.subj.neg$smkyears.cat, missing.include=T)
CrossTable(all.subj.neg$cpd.cat, missing.include=T)
CrossTable(all.subj.neg$emp, missing.include=T)
quantile(all.subj.neg$prescr.1yrisk.T0, probs=c(0.25, 0.5, 0.75))  # median IQR of pre-screening risk

# Other numbers
with(filter(data_screen_abn_neg, interval==2), CrossTable(screen.comb, case))   # numbers considered for Markov assumption test
c(sum(data_interval_abn$case), nrow(data_interval_abn), sum(data_interval_abn$case)/nrow(data_interval_abn)) # interval cancers: cases, # at risk, overall risk
c(sum(data_screen_abn_neg$case), nrow(data_screen_abn_neg), sum(data_screen_abn_neg$case)/nrow(data_screen_abn_neg)) # next-screen cancers after negative: cases, # at risk, overall risk
range_without_outliers(all.subj.neg$prescr.1yrisk)


#### Model development: Interval cancer among NEGATIVES ####

# Interval cancers: overall model (no abnormalities) - properties of screening
# Confirm that pre-screening risk improves the model
int.nopsr <- glm(case ~ 1, data=data_interval_abn, family=binomial(link='log'))
summary(int.nopsr)
int.psr <- glm(case ~ log1yrisk+1, data=data_interval_abn, family=binomial(link='log'))
summary(int.psr)
1-pchisq(int.nopsr$deviance - int.psr$deviance, length(int.psr$coefficients)-length(int.nopsr$coefficients))  # LRT
# Overall model results
  # glm.interval.abn <- glm(case ~ log1yrisk -1, data=data_interval_abn, family=binomial(link='log'))   # run above in setup
  # data_interval_abn$post.risk.interv <- fitted.values(glm.interval)                                   # run above in setup
summary(glm.interval)
confint(glm.interval)
# Does the risk coefficient differ by interval? No (LRT p=0.23). Steps below: fit model, estimate 3 exponents, get p-value, get counts
glm.interval.intervals <- glm(case ~ log1yrisk + log1yrisk:I(as.numeric(interval==2)) + log1yrisk:I(as.numeric(interval==3)) -1, data=data_interval_abn, family=binomial(link='log'))
c(coefficients(glm.interval.intervals)[1], coefficients(glm.interval.intervals)[1]+coefficients(glm.interval.intervals)[2], coefficients(glm.interval.intervals)[1]+coefficients(glm.interval.intervals)[3])
1-pchisq(glm.interval$deviance - glm.interval.intervals$deviance, length(glm.interval.intervals$coefficients)-length(glm.interval$coefficients))
with(data_interval_abn, table(interval))
# Do previous screens matter? No (LRT p=0.99)
glm.int.2levels <- glm(case ~ log1yrisk:as.factor(screen.hist) -1, data=filter(data_interval_abn, interval %in% c(2,3)), family=binomial(link='log'))
summary(glm.int.2levels)
confint(glm.int.2levels)
glm.int.1level <- glm(case ~ log1yrisk -1, data=filter(data_interval_abn, interval %in% c(2,3) & !is.na(screen.hist)), family=binomial(link='log'))
summary(glm.int.1level)
1-pchisq(glm.int.1level$deviance-glm.int.2levels$deviance, df=length(glm.int.2levels$coefficients-length(glm.int.1level$coefficients)))

# Interval cancers: effects of abnormalities
  # Following a negative screen, the relevant CT features are in abnlist.neg. The relevant p-value is for the interaction (i.e. risk differs between the 0 and 1 levels)
# Backwards stepwise selection: selects other.above, benign.nodule, consolidation, adenopathy
int.full <- glm(paste("case ~ logit1yrisk -1 +",paste(abnlist.neg.int, collapse="+"),sep=""), data=data_interval_abn, family=binomial(link='logit'))
bsw.int <- step(int.full, direction="backward", scope = list(lower = case ~ logit1yrisk -1, upper = int.full))
  # Look at a model including these 4 effects
summary(glm(case ~ log1yrisk + log1yrisk:other.above + log1yrisk:benign.nodule + log1yrisk:consolidation + log1yrisk:adenopathy -1, data=data_interval_abn, family=binomial(link='log')))
# Lasso using intermediate lambda: selects adenopathy and consolidation
set.seed(61116)
x  <- model.matrix(case ~ logit1yrisk -1 + logit1yrisk:., data = data_interval_abn[,c("case","logit1yrisk",abnlist.neg)])
cv.lasso <- cv.glmnet(x, data_interval_abn$case, alpha=1, family="binomial")
out <- glmnet(x, data_interval_abn$case, alpha=1, family="binomial")
predict(out, type="coefficients", s=(cv.lasso$lambda.min+cv.lasso$lambda.1se)/2)
  # Look at a model including these two effects
summary(glm(case ~ log1yrisk + log1yrisk:consolidation + log1yrisk:adenopathy -1, data=data_interval_abn, family=binomial(link='log')))
# Based on discussion with Chris: include adenopathy and consolidation. Model as 1 variable (effect size is the same)
# Switch back to log scale for final model for interpretability (these models are run above in data setup section)
    # glm.int.abn.log <- glm(case ~ log1yrisk + log1yrisk:adenop.consol -1, data=data_interval_abn, family=binomial(link='log'))
    # data_interval_abn$post_risk_abn <- fitted.values(glm.int.abn.log)
summary(glm.int.abn)
    # Get estimate and CIs for the exponents
mat <- c(1,0)  # Use this matrix for "Neither noted"
mat <- c(1,1)  # Use this matrix for adenopathy or consolidation
stder <- sqrt(c(t(mat) %*% vcov(glm.int.abn) %*% mat))
c(coefficients(glm.int.abn) %*% mat, (coefficients(glm.int.abn) %*% mat)-1.96*stder, (coefficients(glm.int.abn) %*% mat)+1.96*stder)
# Check for residual effects of LCRAT variables using likelihood ratio tests - the LRT for emp is 0.02, but the Wald is 0.06. We will say p>0.05.
titles <- c("var","null model # param", "extended model # param", "LRT p-value", "check same # obs")
mat.out.interv <- matrix(rep(NA),nrow=length(varlist),ncol=length(titles))
for (x in seq_along(varlist)) {
  mod.without <- glm(case ~ logit1yrisk + logit1yrisk:adenop.consol -1, data=data_interval_abn, family=binomial(link='logit'))
  mod.with <- glm(substitute(case ~ logit1yrisk + logit1yrisk:adenop.consol + logit1yrisk:i -1, list(i=as.name(varlist[x]))), data=data_interval_abn, family=binomial(link='logit'))
  print(summary(mod.with))
  mat.out.interv[x,] <- c(varlist[x], length(mod.without$coefficients), sum(!is.na(mod.with$coefficients)), 1-pchisq(mod.without$deviance-mod.with$deviance, df=sum(!is.na(mod.with$coefficients))-length(mod.without$coefficients)), I(length(mod.without$residuals)==length(mod.with$residuals)))
}
rbind(titles, mat.out.interv)


#### Model development: Next-screen cancer among NEGATIVES ####

# Overall model for next-screen cancer among negatives (no abnormalities) - properties of screening
# Confirm that pre-screening risk improves the model
ns.nopsr <- glm(case ~ 1, data=data_screen_abn_neg, family=binomial(link='log'))
summary(ns.nopsr)
ns.psr <- glm(case ~ log1yrisk+1, data=data_screen_abn_neg, family=binomial(link='log'))
summary(ns.psr)
1-pchisq(ns.nopsr$deviance - ns.psr$deviance, length(ns.psr$coefficients)-length(ns.nopsr$coefficients))
# Overall model results
  # glm.screen.neg <- glm(case ~ log1yrisk -1, data=data_screen_abn_neg, family=binomial(link='log'))  # this is run above the line
  # data_screen_abn_neg$post.risk.neg.overall <- fitted.values(glm.screen.neg)
summary(glm.screen.neg)
confint(glm.screen.neg)
# Does the interval matter? no (p=0.38)
glm.screen.neg.by.int <-  glm(case ~ log1yrisk:as.factor(interval) -1, data=data_screen_abn_neg, family=binomial(link='log'))
summary(glm.screen.neg.by.int)
1-pchisq(glm.screen.neg$deviance - glm.screen.neg.by.int$deviance, length(glm.screen.neg.by.int$coefficients) - length(glm.screen.neg$coefficients))
with(data_screen_abn_neg, table(interval))
# Do previous screens matter? no (p=0.26)
glm.screen.neg.2levels <- glm(case ~ log1yrisk:screen.comb -1, data=filter(data_screen_abn_neg, interval==2 & !is.na(screen.comb)), family=binomial(link='log'))
summary(glm.screen.neg.2levels)
confint(glm.screen.neg.2levels)
glm.screen.neg.1level <- glm(case ~ log1yrisk -1, data=filter(data_screen_abn_neg, interval==2 & !is.na(screen.comb)), family=binomial(link='log'))
summary(glm.screen.neg.1level)
1-pchisq(glm.screen.neg.1level$deviance - glm.screen.neg.2levels$deviance, length(glm.screen.neg.2levels$coefficients) - length(glm.screen.neg.1level$coefficients))
# Do previous screens matter if we ignore pre-screening risk? p=0.14
glm.screen.neg.2levels.nopsr <- glm(case ~ screen.comb -1, data=filter(data_screen_abn_neg, interval==2 & !is.na(screen.comb)), family=binomial(link='log'))
summary(glm.screen.neg.2levels.nopsr)
glm.screen.neg.1level.nopsr <- glm(case ~ 1, data=filter(data_screen_abn_neg, interval==2 & !is.na(screen.comb)), family=binomial(link='log'))
summary(glm.screen.neg.1level.nopsr)
1-pchisq(glm.screen.neg.1level.nopsr$deviance - glm.screen.neg.2levels.nopsr$deviance, length(glm.screen.neg.2levels.nopsr$coefficients) - length(glm.screen.neg.1level.nopsr$coefficients))
# Analysis to compare with Patz result
# Start with everyone with T0 screen being negative, then compare those who are T1-false-positive or
  # negative for the risk of screen-detected cancer at T2. 
patz.m1 <- glm(case ~ log1yrisk -1, data=filter(data.screen.abn, interval==2 & screen.comb %in% c("Neg-Pos","Neg-Neg")), family=binomial(link='log'))
summary(patz.m1)
patz.m2 <- glm(case ~ log1yrisk:screen.comb -1, data=filter(data.screen.abn, interval==2 & screen.comb %in% c("Neg-Pos","Neg-Neg")), family=binomial(link='log'))
summary(patz.m2)
1-pchisq(patz.m1$deviance - patz.m2$deviance, length(patz.m2$coefficients) - length(patz.m1$coefficients))


# Effects of abnormalities for next-screen among negatives
# Backwards stepwise: selects nod6.not.susp, opac.fibr, consolidation, emphysema
scr.neg.full <- glm(paste("case ~ logit1yrisk -1 +",paste(abnlist.neg.int, collapse="+"),sep=""), data=data_screen_abn_neg, family=binomial(link='logit'))
bsw.scr.neg <- step(scr.neg.full, direction="backward", scope = list(lower = case ~ logit1yrisk -1, upper = scr.neg.full))
  # Look at a model including these 4 effects
summary(glm(case ~ log1yrisk + log1yrisk:opac.fibr + log1yrisk:nod6.not.susp + log1yrisk:consolidation + log1yrisk:emphysema -1, data=data_screen_abn_neg, family=binomial(link='log')))
# Lasso using intermediate lambda: selects ONLY logit1yrisk
set.seed(61116)
x  <-  model.matrix(case ~ logit1yrisk -1 + logit1yrisk:. , data = data_screen_abn_neg[,c("case","logit1yrisk",abnlist.neg)])
cv.lasso <- cv.glmnet(x, data_screen_abn_neg$case, alpha=1, family="binomial")
out <- glmnet(x, data_screen_abn_neg$case, alpha=1, family="binomial")
predict(out, type="coefficients", s=(cv.lasso$lambda.min+cv.lasso$lambda.1se)/2)
# From discussion with Chris: keep consolidation and emphysema.
# Switch back to log scale for final model for interpretability (this model is run above in data setup)
    # glm_screen_neg_abn <- glm(case ~ log1yrisk + log1yrisk:consolidation + log1yrisk:emphysema -1, data=data_screen_abn_neg, family=binomial(link='log'))
    # data_screen_abn_neg$post_risk_abn <- fitted.values(glm_screen_neg_abn)
summary(glm_screen_neg_abn)
      # Get estimates and CIs for the exponents
mat <- c(1,0,0) # Use this matrix for "neither noted"
mat <- c(1,1,0) # Use this matrix for consolidation
mat <- c(1,0,1) # Use this matrix for emphysema
stder <- sqrt(c(t(mat) %*% vcov(glm_screen_neg_abn) %*% mat))
c(coefficients(glm_screen_neg_abn) %*% mat, (coefficients(glm_screen_neg_abn) %*% mat)-1.96*stder, (coefficients(glm_screen_neg_abn) %*% mat)+1.96*stder)
# Check for residual effects of LCRAT variables. All p>0.05
titles <- c("var","null model # param", "extended model # param", "LRT p-value", "check same # obs")
mat.out.ns <- matrix(rep(NA),nrow=length(varlist),ncol=length(titles))
for (x in seq_along(varlist)) {
  mod.without <- glm(case ~ logit1yrisk + logit1yrisk:consolidation + logit1yrisk:emphysema -1, data=data_screen_abn_neg, family=binomial(link='logit'))
  mod.with <- glm(substitute(case ~ logit1yrisk + logit1yrisk:consolidation + logit1yrisk:emphysema + logit1yrisk:i -1, list(i=as.name(varlist[x]))), data=data_screen_abn_neg, family=binomial(link='logit'))
  mat.out.ns[x,] <- c(varlist[x], length(mod.without$coefficients), sum(!is.na(mod.with$coefficients)), 1-pchisq(mod.without$deviance-mod.with$deviance, df=sum(!is.na(mod.with$coefficients))-length(mod.without$coefficients)), I(length(mod.without$residuals)==length(mod.with$residuals)))
}
rbind(titles, mat.out.ns)
# What if we account for screening history along with pre-screening risk, emphysema, and consolidation? p=0.34
m1 <- glm(case ~ log1yrisk + log1yrisk:consolidation + log1yrisk:emphysema -1, data=filter(data_screen_abn_neg, interval==2 & !is.na(screen.comb)), family=binomial(link='log'))
summary(m1)
m2 <- glm(case ~ log1yrisk:screen.comb + log1yrisk:consolidation + log1yrisk:emphysema -1, data=filter(data_screen_abn_neg, interval==2 & !is.na(screen.comb)), family=binomial(link='log'))
summary(m2)
1-pchisq(m1$deviance - m2$deviance, length(m2$coefficients) - length(m1$coefficients))






#### Additional analyses (GEE, AUCs, cross-validation, etc) ####

### Comparison with GEE - this impacts the SEs negligibly ###
# Interval cancer models
summary(geeglm(case ~ log1yrisk -1, id=pid, data=data.interval, family=binomial(link='log'), corstr="exchangeable", waves=interval))
summary(glm.interval)
summary(geeglm(case ~ log1yrisk + log1yrisk:adenop.consol -1, id=pid, data=data_interval_abn, family=binomial(link='log'), corstr="exchangeable", waves=interval))
summary(glm.int.abn.log)
# Next-screen among negatives model
summary(geeglm(case ~ log1yrisk -1, id=pid, data=data_screen_abn_neg, family=binomial(link='log'), corstr="exchangeable", waves=interval))
summary(glm.screen.neg)
summary(geeglm(case ~ log1yrisk + log1yrisk:consolidation + log1yrisk:emphysema -1, id=pid, data=data_screen_abn_neg, family=binomial(link='log'), corstr="exchangeable", waves=interval))
summary(glm_screen_neg_abn)


### Calibration and validation analyses ###
# 10-fold cross-validated calibration
# Interval cancers
set.seed(61116)
data_interval_abn$randgrp <- base::sample(1:10, nrow(data_interval_abn), replace=T)
data_interval_abn$cvpred <- NA
for (i in 1:10) {
  fit <- glm(formula = case ~ log1yrisk + log1yrisk:adenop.consol - 1, 
             family = binomial(link = "log"), data = filter(data_interval_abn, randgrp!=i))
  data_interval_abn[data_interval_abn$randgrp==i,]$cvpred <- predict(fit, newdata=data_interval_abn[data_interval_abn$randgrp==i,], type="response")
}
data_interval_abn <- mutate(data_interval_abn, cvpred.ntile = ntile(cvpred, 5))
data_interval_abn %>% group_by(cvpred.ntile) %>% summarise(pred.cases= sum(cvpred), obs.cases = sum(case))
c(sum(data_interval_abn$cvpred), sum(data_interval_abn$case)) # number obs and expected cases
poisson.test(round(sum(data_interval_abn$cvpred),0), sum(data_interval_abn$case), alternative="two.sided") # p-value, requires rounding

# Next-screen among screen-negatives
set.seed(61116)
data_screen_abn_neg$randgrp <- base::sample(1:10, nrow(data_screen_abn_neg), replace=T)
data_screen_abn_neg$cvpred <- NA
for (i in 1:10) {
  fit <- glm(formula = case ~ log1yrisk + log1yrisk:consolidation + log1yrisk:emphysema - 1, 
             family = binomial(link = "log"), data = filter(data_screen_abn_neg, randgrp!=i))
  data_screen_abn_neg[data_screen_abn_neg$randgrp==i,]$cvpred <- predict(fit, newdata=data_screen_abn_neg[data_screen_abn_neg$randgrp==i,], type="response")
}
data_screen_abn_neg <- mutate(data_screen_abn_neg, cvpred.ntile = ntile(cvpred, 5))
data_screen_abn_neg %>% group_by(cvpred.ntile) %>% summarise(pred.cases= sum(cvpred), obs.cases = sum(case))
c(sum(data_screen_abn_neg$cvpred), sum(data_screen_abn_neg$case)) # number obs and expected cases
poisson.test(round(sum(data_screen_abn_neg$cvpred)), sum(data_screen_abn_neg$case), alternative="two.sided")

# 10-fold cross validation to get CV error. The first delta is standard version; second is bias-corrected.
set.seed(61116)
cv.err.int <- cv.glm(data_interval_abn, glm.int.abn, K=10)
cv.err.int$delta
cv.err.screen.neg <- cv.glm(data_screen_abn_neg, glm_screen_neg_abn, K=10)
cv.err.screen.neg$delta


### Calculate AUCs ###
## Regular AUCs. By default, the 95% CI are computed with 2000 stratified bootstrap replicates.
library(pROC)
# Interval cancer model
with(filter(data_interval_abn, interval==1), roc(case, post_risk_abn, ci=T, plot=T))  # T0-T1 - change interval for T1-T2, post-T2
# Next-screen model among negatives
with(filter(data_screen_abn_neg, interval==1), roc(case, post_risk_abn, ci=T, plot=T))   # T1 - change interval for T2
## Optimism-corrected AUCs - have to use logistic models for this, and have to actually add the interaction terms to the dataset.
library(rms)
data_screen_abn_neg <- mutate(data_screen_abn_neg, logit1yriskconsolidation = logit1yrisk*consolidation, logit1yriskemphysema = logit1yrisk*emphysema)
data_interval_abn <- mutate(data_interval_abn, logit1yriskadenopconsol = logit1yrisk*adenop.consol)
# Interval cancer model
mod.int <- lrm(case ~ logit1yrisk + logit1yriskadenopconsol -1, x=T, y=T, data=data_interval_abn)
set.seed(61116)
validate(mod.int, B=1000)
c(0.5*(0.5072+1), 0.5*(0.5082+1)) # AUC = 0.5(Dxy+1). Naive, optimism-corrected AUCs - 0.75 is OC-AUC
# Next-screen model among negatives
mod.ns <- lrm(case ~ logit1yrisk + logit1yriskconsolidation + logit1yriskemphysema -1, x=T, y=T, data=data_screen_abn_neg)
set.seed(61116)
validate(mod.ns, B=1000)     # AUC = 0.5(Dxy+1)
c(0.5*(0.4760+1), 0.5*(0.4689+1)) # AUC = 0.5(Dxy+1). Naive, optimism-corrected AUCs - 0.73 is OC-AUC

#### Figures #### 

# Effect of screen findings on risk of INTERVAL ca among screen-negatives
med.risk.interv.prescr <- median(filter(data_interval_abn, interval==1)$prescr.1yrisk, na.rm=T)
med.risk.interv.post.noac <- median(filter(data_interval_abn, interval==1 & adenop.consol==0)$post_risk_abn, na.rm=T)
med.risk.interv.post.ac <- median(filter(data_interval_abn, interval==1 & adenop.consol==1)$post_risk_abn, na.rm=T)
png(file="/Users/hrobbins827/Documents/PhD/NCI overflow/NLST/Figures/screen_neg_interval_ad_con.png",width=1200,height=850)
ggplot() +
  theme(panel.background = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.title = element_text(size=28),
        axis.text = element_text(colour="black", size=28)) +
  scale_y_continuous(labels=scales::percent, limits=c(0, 0.045)) +
  scale_x_continuous(breaks=NULL) +
  # scale_x_continuous(breaks=c(0,1), labels = c("Negative screen", "1-year interval")) +
  ylab("1-year lung cancer risk, %\n") + xlab("") +
  geom_boxplot(data=subset(data_interval_abn, interval==1), aes(x=0, y=prescr.1yrisk), lwd=1, width=0.4, outlier.shape=NA) +
  geom_boxplot(data=subset(data_interval_abn, interval==1 & adenop.consol==0), aes(x=0.9, y=post_risk_abn), lwd=1, width=0.98*.8, outlier.shape=NA) +
  geom_boxplot(data=subset(data_interval_abn, interval==1 & adenop.consol==1), aes(x=1.1, y=post_risk_abn), lwd=1, width=0.02*.8, outlier.shape=NA) +
  geom_segment(aes(x=0, y=med.risk.interv.prescr, xend=0.9, yend=med.risk.interv.post.noac), linetype="dashed", size=0.6) +
  geom_segment(aes(x=0, y=med.risk.interv.prescr, xend=1.1, yend=med.risk.interv.post.ac), linetype="dashed", size=0.6) +
  annotate(geom="text", x=0.6, y=0.0062, label = "Adenopathy or consolidation (2%)", angle=4, size=9) +
  annotate(geom="text", x=0.45, y=0.0028, label = "Neither noted (98%)", angle=-4, size=9) +
  annotate(geom="text", x=0, y=0.045, label="Pre-screening risk", size=9) +
  annotate(geom="text", x=1, y=0.045, label="Risk during 1-year interval", size=9)
dev.off()
# Numbers for the text
with(data_interval_abn, CrossTable(interval, adenop.consol)) # 1.8% have adenop or consol at T0
# Some percentiles for below the figure
neg.i.psr.q <- quantile(filter(data_interval_abn, interval==1)$prescr.1yrisk, probs=c(0.25, 0.5, 0.75))
neg.i.no.q <- quantile(filter(data_interval_abn, interval==1 & adenop.consol==0)$post_risk_abn, probs=c(0.25, 0.5, 0.75))
neg.i.adcon.q <- quantile(filter(data_interval_abn, interval==1 & adenop.consol==1)$post_risk_abn, probs=c(0.25, 0.5, 0.75))
rbind(neg.i.psr.q, neg.i.no.q, neg.i.adcon.q)  # print the quantiles for each group
c(neg.i.no.q[2]/neg.i.psr.q[2], neg.i.adcon.q[2]/neg.i.psr.q[2]) # median RRs for no, adenop.consol
c(neg.i.no.q[2]-neg.i.psr.q[2], neg.i.adcon.q[2]-neg.i.psr.q[2]) # median RDs for no, adenop.consol
range_without_outliers(filter(data_interval_abn, interval==1 & adenop.consol==0)$post_risk_abn) # this uses my function defined in hilary_functions.R

# Effect of screen findings on risk of SCREEN-DETECTED ca among screen negatives
    # Note: I am not making a boxplot for the N=36 with emphysema and consolidation.
med.risk.screen.neg.prescr <- median(filter(data_screen_abn_neg, interval==1)$prescr.1yrisk)
med.risk.screen.neg.neither <- median(filter(data_screen_abn_neg, interval==1 & emphysema==0 & consolidation==0)$post_risk_abn)
med.risk.screen.neg.emp <-  median(filter(data_screen_abn_neg, interval==1 & emphysema==1)$post_risk_abn)
med.risk.screen.neg.consol <-  median(filter(data_screen_abn_neg, interval==1 & consolidation==1)$post_risk_abn)
png(file="/Users/hrobbins827/Documents/PhD/NCI overflow/NLST/Figures/screen_neg_emp_consol.png",width=1200,height=850)
ggplot() +
  theme(panel.background = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.title = element_text(size=28),
        axis.text = element_text(colour="black", size=28)) +
  scale_y_continuous(labels=scales::percent, limits=c(0, 0.045)) +
  scale_x_continuous(breaks=NULL) +
  # scale_x_continuous(breaks=c(0,1), labels = c("Negative screen", "Screen")) +   # remove labels for manuscript
  ylab("1-year lung cancer risk, %\n") + xlab("") +
  geom_boxplot(data=subset(data_screen_abn_neg, interval==1), aes(x=0, y=prescr.1yrisk), lwd=1, width=0.4, outlier.shape=NA) +
  geom_boxplot(data=subset(data_screen_abn_neg, interval==1 & emphysema==0 & consolidation==0), aes(x=0.8, y=post_risk_abn), lwd=1, width=0.69*1.2, outlier.shape=NA) +
  geom_boxplot(data=subset(data_screen_abn_neg, interval==1 & emphysema==1), aes(x=1, y=post_risk_abn), lwd=1, width=.3*1.2, outlier.shape=NA) +
  geom_boxplot(data=subset(data_screen_abn_neg, interval==1 & consolidation==1), aes(x=1.2, y=post_risk_abn), lwd=1, width=.01*1.2, outlier.shape=NA) +
  geom_segment(aes(x=0, y=med.risk.screen.neg.prescr, xend=0.8, yend=med.risk.screen.neg.neither), linetype="dashed", size=0.6) +
  geom_segment(aes(x=0, y=med.risk.screen.neg.prescr, xend=1.0, yend=med.risk.screen.neg.emp), linetype="dashed", size=0.6) +
  geom_segment(aes(x=0, y=med.risk.screen.neg.prescr, xend=1.2, yend=med.risk.screen.neg.consol-0.0019), linetype="dashed", size=0.6) +
  annotate(geom="text", x=0.3, y=0.0003, label = "Neither (70%)", angle=-4, size=9) +
  annotate(geom="text", x=0.5, y=0.0052, label = "Emphysema (30%)", angle=2, size=9) +
  annotate(geom="text", x=0.65, y=0.0105, label = "Consolidation (0.6%)", angle=12, size=9) +
  annotate(geom="text", x=0, y=0.045, label="Pre-screening risk", size=9) +
  annotate(geom="text", x=1, y=0.045, label="Risk at next screen", size=9)
dev.off()
# Numbers for the text: Among T0-negatives: prevalence of self-reported emphysema (7%), CT emphysema (30%), and consolidation (0.6%)
with(filter(data_screen_abn_neg, interval==1), CrossTable(emp))
with(filter(data_screen_abn_neg, interval==1), CrossTable(emphysema))
with(filter(data_screen_abn_neg, interval==1), CrossTable(consolidation))
with(filter(data_screen_abn_neg, interval==1), CrossTable(emphysema, consolidation))
with(data_screen_abn_neg, CrossTable(interval, I(emphysema==0 & consolidation==0))) # no emphysema NOR consolidation (70% at T0)
  # Thus .738 *.696 = 51% of participants would be screen-negative without emp or consol.
# Some percentiles for below the figure
quantile(filter(data_screen_abn_neg, interval==1)$post_risk_abn, probs=c(0.25, 0.5, 0.75, 0.8, 0.85, 0.9)) # overall quantiles among all negatives
neg.s.psr.q <- quantile(filter(data_screen_abn_neg, interval==1)$prescr.1yrisk, probs=c(0.25, 0.5, 0.75))
neg.s.no.no.q <- quantile(filter(data_screen_abn_neg, interval==1 & emphysema==0 & consolidation==0)$post_risk_abn, probs=c(0.25, 0.5, 0.75))
neg.s.emp.q <- quantile(filter(data_screen_abn_neg, interval==1 & emphysema==1)$post_risk_abn, probs=c(0.25, 0.5, 0.75))
neg.s.consol.q <- quantile(filter(data_screen_abn_neg, interval==1 & consolidation==1)$post_risk_abn, probs=c(0.25, 0.5, 0.75))
rbind(neg.s.psr.q, neg.s.no.no.q, neg.s.emp.q, neg.s.consol.q)   # print the quantiles for each group
rbind(neg.s.no.no.q[2]/neg.s.psr.q[2], neg.s.emp.q[2]/neg.s.psr.q[2], neg.s.consol.q[2]/neg.s.psr.q[2]) # median RRs for no-no, emp, consol
rbind(neg.s.no.no.q[2]-neg.s.psr.q[2], neg.s.emp.q[2]-neg.s.psr.q[2], neg.s.consol.q[2]-neg.s.psr.q[2]) # median RDs for no-no, emp, consol


# Potential risk THRESHOLDS for longer interval after negative screen. Use T0-negatives and cases at T1
  # Note: 64 cancers at T1; only 30 among the 70% with no emp or consol.
# Make a dataset for impact of different risk thresholds among all screen-negatives
num_T0_neg <- nrow(filter(data_screen_abn_neg, interval==1))        # change to interval==2 to look at T2
num_T0_all <- nrow(filter(data.screen.abn, interval==1))            # change to interval==2 to look at T2
cases_T1 <- sum(filter(data_screen_abn_neg, interval==1)$case)      # change to interval==2 to look at T2
n <- perc_negs <- perc_all <- ca_N <- perc_of_ca <- perc_w_ca <- perc_w_emp <- perc_w_consol <- perc_w_emp_or_consol <- vector()
thresholds <- seq(0,0.1,0.0001)
for (i in seq_along(thresholds)) {
  dat <- filter(data_screen_abn_neg, interval==1 & post_risk_abn<=thresholds[i])  # change to interval==2 to look at T2
  n[i] <- nrow(dat)
  perc_negs[i] <- 100*nrow(dat)/num_T0_neg   # percent of all negatives below threshold
  perc_all[i] <- 100*nrow(dat)/num_T0_all    # percent of all individuals below threshold
  ca_N[i] <- sum(dat$case)                   # number with cancer below threshold
  perc_of_ca[i] <- 100*sum(dat$case)/cases_T1  # percent of cancers falling below threshold
  perc_w_ca[i] <- 100*sum(dat$case)/nrow(dat)  # percent of individuals below threshold who have cancer
  perc_w_emp[i] <- 100*sum(dat$emphysema)/nrow(dat)  # percent of individuals below threshold who have CT-emphysema
  perc_w_consol[i] <- 100*sum(dat$consolidation)/nrow(dat) # percent of individuals below threshold who have CT-consolidation
  perc_w_emp_or_consol[i] <- 100*sum(I(dat$emphysema==1 | dat$consolidation==1))/nrow(dat)  # percent of individuals below threshold who have CT-emphysema or consolidation
}
thres.plot.all <- as.data.frame(cbind(threshold=100*thresholds, n, perc_negs, perc_all, ca_N, perc_of_ca, perc_w_ca, perc_w_emp, perc_w_consol, perc_w_emp_or_consol))
# Plot this
thres.of.interest <- c(0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8)
png(file="/Users/hrobbins827/Documents/PhD/NCI overflow/NLST/Figures/negatives_threshold.png",width=1200,height=850)
ggplot() + geom_line(data=thres.plot.all, aes(x=perc_negs, y=perc_of_ca), size=0.8) + 
  theme(panel.background = element_rect(fill=NA), panel.grid.major=element_line(colour="grey88"), panel.grid.minor=element_line(colour="grey88"),
        axis.line = element_line(colour="black"), axis.title = element_text(size=28), axis.text = element_text(colour="black", size=28)) +
  ylab("% of detectable next-screen cancers with delayed diagnosis") + xlab("\n% of screen-negatives with longer-than-annual interval") +
  geom_point(data=subset(thres.plot.all, threshold %in% thres.of.interest), aes(x=perc_negs, y=perc_of_ca), size=4.5) +
  geom_text(data=subset(thres.plot.all, threshold %in% thres.of.interest[1:5]), aes(x=perc_negs, y=perc_of_ca, label=paste("r \u2264", as.character(threshold), "%", sep="")), size=9, vjust=-1, hjust=0.8) +
  geom_text(data=subset(thres.plot.all, threshold %in% thres.of.interest[6:8]), aes(x=perc_negs, y=perc_of_ca, label=paste("r \u2264", as.character(threshold), "%", sep="")), size=9, hjust=-0.25)
dev.off()
# Print some numbers to highlight in the text
filter(thres.plot.all, threshold %in% thres.of.interest)



