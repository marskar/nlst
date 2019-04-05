#### Setup ####

# Read in library functions including Stephanie's coxph.risk which I installed from the local tar.gz file
packages <- c(
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
    "boot"
    )

not_installed <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(not_installed)) install.packages(not_installed)
all(lapply(packages, require, character.only = TRUE))

# Read in the Kovalchik prediction function
source(here("R/kovalchik.R"))

devtools::install_github('marskar/coxph_risk')
devtools::install_github('marskar/lcmodels')

# Function to create "not in" operator
'%!in%' <- function(x,y)!('%in%'(x,y))

# Load PLCO and NLST data sets. New versions provided by Li Cheung on 11 July 2016, 14 July 2016, 20 July 2016, 8 Aug 2016.
# Load NLST and PLCO data
nlst = readRDS('data/nlst.rds')
plco = readRDS('data/plco.rds')

# In the PLCO dataset, impute missing family history values to 0
plco$fam.lung.trend <- ifelse(is.na(plco$fam.lung.trend), 0, plco$fam.lung.trend)
plco_control <- subset(plco, control.group==1) # control arm of PLCO who had no chest xray

# Remove people with <30 pack-years and age<55 or age>74 from NLST
nlst <- nlst %>%
    filter(pkyears.cat != "[0,30)" & age >= 55 & age <= 74) %>%
    # Make a new pack-years variable to get rid of the [0,30) level
    mutate(pkyears.cat.clone = ifelse(
        pkyears.cat == "[30,40)",
        "[30,40)",
        ifelse(
            pkyears.cat == "[40,50)",
            "[40,50)",
            ifelse(pkyears.cat == "[50,Inf]", "[50,Inf]", NA)
        )
    )) %>%
    mutate(pkyears.cat = as.factor(pkyears.cat.clone)) %>%
    # Make a variable for days to diagnosis
    mutate(days_to_dx = ifelse(case == 1, 365 * incidence.years, NA)) %>%
    identity()

# Make a subset of NLST data with the LCRAT variables that we will need later
varlist <-
    c(
        "female",
        "race",
        "edu6",
        "fam.lung.trend",
        "emp",
        "bmi",
        "cpd",
        "pkyears.cat",
        "age",
        "qtyears",
        "smkyears"
    )
nlst.sub <- as.data.frame(cbind(nlst[,varlist], pid=nlst$pid, lss=as.numeric(nlst$lss)))

# To later calculate pre-screening risk, we must first fit the incidence model and other-cause death models in PLCO.
LCRAT <- coxph(Surv(incidence.years, case) ~ 
                 female+race+edu6+fam.lung.trend+emp+I(bmi<=18.5)+I(cpd>20)+as.factor(pkyears.cat)+
                 I(log(age))+I(log(bmi))+I(log(qtyears+1))+smkyears,data = plco_control)
cox.death <- coxph(Surv(years.followed, other.cause.death) ~ 
                     female+race+edu6+emp+I(bmi <= 18.5)+I(cpd>20)+as.factor(pkyears.cat)+I((age)^2)+I((bmi-25)^2)+
                     I(log(qtyears+1))+smkyears, data = plco_control)

lcratvars <-
    c(
        "age",
        "female",
        "smkyears",
        "qtyears",
        "cpd",
        "race",
        "emp",
        "fam.lung.trend",
        "bmi",
        "edu6",
        "pkyears.cat"
    )

nlst_lcrat <- as.data.frame(cbind(nlst[, lcratvars], pid=nlst$pid, lss=as.numeric(nlst$lss)))
lcmodels::lcmodels(nlst_lcrat)

# TODO Data should include all positives, not just false positives
# positives <-
#     nlst %>% filter_at(
#         vars(
#             truefalse_scrnres_ly0,
#             truefalse_scrnres_ly1,
#             truefalse_scrnres_ly2
#         ),
#         all_vars(. %in% c(1, 2, 3))
#     )
# positives

# Subset to CT arm in NLST and make a pos/neg variable for the first, second, and third screens
nlst <- nlst %>% mutate(T0posneg = case_when(
    truefalse_scrnres_ly0 %in% c(4, 5, 6) ~ 0,
    truefalse_scrnres_ly0 %in% c(1, 2, 3) ~ 1
)) %>% mutate(T1posneg = case_when(
    truefalse_scrnres_ly1 %in% c(4, 5, 6) ~ 0,
    truefalse_scrnres_ly1 %in% c(1, 2, 3) ~ 1
)) %>% mutate(T2posneg = case_when(
    truefalse_scrnres_ly2 %in% c(4, 5, 6) ~ 0,
    truefalse_scrnres_ly2 %in% c(1, 2, 3) ~ 1
))

# Subset to CT arm and create screening history variables
nlst_CT <- nlst %>% 
    subset(screen_group == "CT") %>%
    mutate(hist.T0.T1 =
        case_when(
            T0posneg == 0 & T1posneg == 0 ~ 1,
            T0posneg == 0 & T1posneg == 1 ~ 2,
            T0posneg == 1 & T1posneg == 0 ~ 3,
            T0posneg == 1 & T1posneg == 1 ~ 4
        )) %>% 
    mutate(hist.T1.T2 =
        case_when(
            T1posneg == 0 & T2posneg == 0 ~ 1,
            T1posneg == 0 & T2posneg == 1 ~ 2,
            T1posneg == 1 & T2posneg == 0 ~ 3,
            T1posneg == 1 & T2posneg == 1 ~ 4
        )) %>% 
    mutate_at(vars(hist.T0.T1, hist.T1.T2),
              funs(factor(
                  .,
                  levels = c(1, 2, 3, 4),
                  labels = c("Neg-Neg", "Neg-Pos", "Pos-Neg", "Pos-Pos")
              ))) %>%
    identity()

### Create dataset for risk from T0 to T1.
  # 0 inadeq, 1 true-pos, 2 poss true-pos, 3 false-pos, 4 true-neg, 5 poss false-neg, 6 false-neg
# At risk for screen-detected at T1: either false-positive or true-negative at T0, and did not have any of the following at T1:
  # inadequate image, left study, refused, wrong screen, erroneous report of LC, form not submitted (no missing values of scr_res0).
  # Case status: case=1 AND either of (true-pos at T1 or T1 is coded as "not expected: cancer/death in screening window")
nlst.CT.T1.scrisk <-
    dplyr::filter(nlst_CT,
                  truefalse_scrnres_ly0 %in% c(2, 3, 4, 5) &
                      scr_res1 %!in% c(10, 11, 15, 17, 95, 97)) %>%
    mutate(case_T1_screen = if_else(case == 1 &
                                        (
                                            truefalse_scrnres_ly1 == 1 | scr_res1 %in% c(23, 24)
                                        ),
                                    1, 0))
### Create dataset for risk from T1 to T2.
# At risk for screen-detected at T2: either false-positive or true-negative at T1, and did not have any of the following at T2:
# inadequate image, left study, refused, wrong screen, erroneous report of LC, form not submitted (no missing values of scr_res0).
# Case status: case=1 AND either of (true-pos at T2 or T2 is coded as "not expected: cancer/death in screening window")
nlst.CT.T2.scrisk <-
    dplyr::filter(nlst_CT,
                  truefalse_scrnres_ly1 %in% c(2, 3, 4, 5) &
                      scr_res2 %!in% c(10, 11, 15, 17, 95, 97)) %>%
    mutate(case_T2_screen = if_else(case == 1 &
                                        (
                                            truefalse_scrnres_ly2 == 1 | scr_res2 %in% c(23, 24)
                                        ),
                                    1, 0))

# Construct dataset to model risk of ALL screen-detected cancers (at T1 and T2)
data.screen <-
    data.frame(
        pid = c(nlst.CT.T1.scrisk$pid, nlst.CT.T2.scrisk$pid),
        case = c(
            nlst.CT.T1.scrisk$case_T1_screen,
            nlst.CT.T2.scrisk$case_T2_screen
        ),
        screen.result = c(nlst.CT.T1.scrisk$T0posneg, nlst.CT.T2.scrisk$T1posneg),
        interval = c(rep(1, times = nrow(nlst.CT.T1.scrisk)), rep(2, times =
                                                                      nrow(nlst.CT.T2.scrisk)))
    ) %>%
    # Merge this back with covariates from NLST
    merge(nlst.sub, by = "pid", all.x = TRUE) %>%
    # Add a variable for lagged screen result & a 6-level variable for all combinations
    group_by(pid) %>%
    mutate(lag.screen = dplyr::lag(screen.result, order_by = interval)) %>%
    mutate(
        screen.comb = 1 * (interval == 1 & screen.result == 0) +
            2 * (interval == 1 &
                     screen.result == 1) + 3 * (interval == 2 &
                                                    lag.screen == 0 &
                                                    screen.result == 0) +
            4 * (interval == 2 &
                     lag.screen == 0 &
                     screen.result == 1) + 5 * (interval == 2 &
                                                    lag.screen == 1 &
                                                    screen.result == 0) +
            6 * (interval == 2 &
                     lag.screen == 1 & screen.result == 1)
    )

data.screen$screen.comb <- factor(data.screen$screen.comb, levels = c(1,2,3,4,5,6),
                                  labels = c("Neg","Pos","Neg-Neg","Neg-Pos","Pos-Neg","Pos-Pos"))
# Update age, quit-years, and smoke-years by adding a year for T1
data.screen <- mutate(data.screen, age=as.numeric(age), smkyears=as.numeric(smkyears), qtyears=as.numeric(qtyears))
data.screen <- mutate(data.screen, age=ifelse(interval==2, age+1, age),
                      smkyears=ifelse(interval==2 & qtyears==0, smkyears+1, smkyears),
                      qtyears=ifelse(interval==2 & qtyears>0, qtyears+1, qtyears))
data.screen <- mutate(data.screen, pkyears.cont=cpd*smkyears/20) # using new smoke-years, update pack-years, then re-categorize
data.screen <- mutate(data.screen, pkyears.cat = 1*(pkyears.cont>=30 & pkyears.cont<40) + 
                        2*(pkyears.cont>=40 & pkyears.cont<50) + 3*(pkyears.cont>=50 & pkyears.cont<999))
data.screen$pkyears.cat <- factor(data.screen$pkyears.cat, levels=1:3, labels=c("[30,40)","[40,50)","[50,Inf]"))

# Load abnormalities data (person-screen level) and merge with data.screen
  # This dataset was generated by the program prepare_abn_data_vX.R (replaced with v5 26 Nov 2018)
abn <- readRDS(here("data/abn_lrads_merged.rds"))
data.screen.abn <- merge(data.screen, abn, by=c("pid","interval"), all.x=TRUE, all.y=FALSE)
# Replace NAs with 0 (not present) for appropriate variables
replacevars  <- names(abn.pl.all)[!names(abn.pl.all) %in% c("pid","interval","LRcat","LRcatcol.neg","LRcatcol.pos")]
data.screen.abn[replacevars][is.na(data.screen.abn[replacevars])]  <- 0
data.screen.abn <- mutate(data.screen.abn, longest.diam = ifelse(is.infinite(longest.diam), 0, longest.diam))
# Make a variable for including observations in Lung-RADS analysis
data.screen.abn$LR.include <- (data.screen.abn$LRcat %in% c("1","2","3","4A","4B","4X"))
# Create variable for log(diameter)
data.screen.abn$log.diam <- log(data.screen.abn$longest.diam+1)
# Calculate pre-screening risk inside this dataset
data.screen.abn$prescr.1yrisk <- risk.kovalchik(0, 1, data.screen.abn, LCRAT, cox.death)
data.screen.abn <- mutate(data.screen.abn, log1yrisk=log(prescr.1yrisk), logit1yrisk=log(prescr.1yrisk/(1-prescr.1yrisk)))
# This datasets is needed to separately model screen-detected cancers incorporating abnormalities for false-positives
data.screen.abn.pos <- dplyr::filter(data.screen.abn, screen.result==1)
# Make a categorical variable for diameter
data.screen.abn.pos <- mutate(data.screen.abn.pos, diam.cat = 1*(longest.diam==0)+2*(longest.diam>0 & longest.diam<=5)+
                        3*(longest.diam>5 & longest.diam<=7) + 4*(longest.diam>7 & longest.diam<=10) +
                          5*(longest.diam>10 & longest.diam<=13) + 6*(longest.diam>13 & longest.diam<100))
data.screen.abn.pos$diam.cat <- factor(data.screen.abn.pos$diam.cat, levels=c(1:6),labels=c("0","4-5","6-7","8-10","11-13","14+"))
# Make another categorical variable for diameter that separates GG nodules by < >= 20 mm (I am not including this in abn.list)
data.screen.abn.pos <- mutate(data.screen.abn.pos, diam.cat.GG = ifelse(any.GG==0 | longest.diam==0, diam.cat,
                        ifelse(any.GG==1 & longest.diam<20, 7, ifelse(any.GG==1 & longest.diam>=20 & longest.diam<99, 8, NA))))
data.screen.abn.pos$diam.cat.GG <- factor(data.screen.abn.pos$diam.cat.GG, levels=c(1:8),
                        labels=c("0","noGG, 4-5","noGG, 6-7","noGG, 8-10","noGG, 11-13","noGG, 14+","GG, <20","GG, >=20"))
# Create a variable for any.growth that reflects a group in which growth can't be assessed (i.e. screen=T0)
data.screen.abn.pos <- mutate(data.screen.abn.pos, growth.3l = 
                          1*(interval==1) + 2*(interval==2 & any.growth==0) + 3*(interval==2 & any.growth==1))
data.screen.abn.pos$growth.3l <- factor(data.screen.abn.pos$growth.3l, levels = c(1,2,3), labels=c("NA","No","Yes"))
# Will need this vector for exploratory analysis of abnormalities (24 Jan 2017 - also create interaction vectors)
abnlist <- abnlist.pos <- c(names(abn.pl.all)[5:32], "diam.cat") # omits longest.diam and any.nodule, but includes diam.cat which subsumes both. 12 Oct 2018 now includes any.new.nodule
abnlist.pos.int <- lapply(abnlist.pos, function(x) {substitute(logit1yrisk:i, list(i=as.name(x)))})

# Make dataset of unique individuals for descriptive table of false-positives
all.subj.pos <- dplyr::filter(nlst.CT, pid %in% data.screen.abn.pos$pid)
all.subj.pos <- mutate(all.subj.pos, age.cat=as.factor(ifelse(age>=55 & age<60, "55-59", ifelse(age>=60 & age<65, "60-64",
                          ifelse(age>=65 & age<70, "65-69", ifelse(age>=70 & age<75, "70-74", NA))))),
                       qtyears.cat=as.factor(ifelse(qtyears==0, "Current smoker", ifelse(qtyears>0 & qtyears<=5, "1-5",
                          ifelse(qtyears>5 & qtyears<=10, "6-10", ifelse(qtyears>10 & qtyears<99, "11 or more", NA))))),
                       bmi.cat=as.factor(ifelse(bmi>0 & bmi<18.5, "Underweight", ifelse(bmi>=18.5 & bmi<25, "Normal",
                          ifelse(bmi>=25 & bmi<30, "Overweight", ifelse(bmi>=30, "Obese", NA))))),
                       cpd.cat=as.factor(ifelse(cpd>0 & cpd<20, "<20", ifelse(cpd>=20 & cpd<30, "20-29",
                          ifelse(cpd>=30 & cpd<40, "30-39", ifelse(cpd>=40 & cpd<99, "40+", NA))))),
                       smkyears.cat=as.factor(ifelse(smkyears>0 & smkyears<30, "<30", ifelse(smkyears>=30 & smkyears<40, "30-39",
                          ifelse(smkyears>=40 & smkyears<50, "40-49", ifelse(smkyears>=50 & smkyears<99, "50+", NA))))))

# Run the main models
# Fit an overall model with one exponent
glm.screen.pos <- glm(case ~ log1yrisk -1, data=data.screen.abn.pos, family=binomial(link='log'))
# Model with abnormalities
glm.screen.pos.abn.log <- glm(case ~ log1yrisk:diam.cat + log1yrisk:any.growth + log1yrisk:any.upper + 
                                log1yrisk:I(any.right.mid==1|any.lingula==1) + log1yrisk:any.mixed + 
                                log1yrisk:any.spiculation + log1yrisk:I(any.poor.def==1|any.margin.unab==1) -1, 
                                data=data.screen.abn.pos, family=binomial(link='log'))
data.screen.abn.pos$post.risk.abn.pos <- fitted.values(glm.screen.pos.abn.log)



# --------------------------- run code to this line for data setup ----------------------------- #


#### Descriptive stats ####

# Descriptive characteristics - table 1
nrow(all.subj.pos)  # number of unique individuals
quantile(all.subj.pos$prescr.1yrisk.T0, probs=c(0.25, 0.5, 0.75))  # median IQR of pre-screening risk
CrossTable((data.screen.abn.pos %>% group_by(pid) %>% summarise(times.in.analysis=n()))$times.in.analysis) # times included in next-screen analysis
CrossTable(all.subj.pos$female, missing.include=T)
CrossTable(all.subj.pos$age.cat, missing.include=T)
CrossTable(all.subj.pos$race, missing.include=T)   # 0 white, 1 black, 2 hispanic, 3 other
CrossTable(all.subj.pos$edu6, missing.include=T)   # see codebook
CrossTable(all.subj.pos$bmi.cat, missing.include=T)
CrossTable(all.subj.pos$fam.lung.trend, missing.include=T)  # none, 1, 2+
CrossTable(all.subj.pos$qtyears.cat, missing.include=T)
CrossTable(all.subj.pos$pkyears.cat, missing.include=T)
CrossTable(all.subj.pos$smkyears.cat, missing.include=T)
CrossTable(all.subj.pos$cpd.cat, missing.include=T)
CrossTable(all.subj.pos$emp, missing.include=T)

# Prevalence of abnormalities - table 2
lapply(data.screen.abn.pos[data.screen.abn.pos$interval==2,abnlist.pos], function(x) prop.table(table(x))) # at T1 abn screen

# Prevalence and risks - table 3
# Get contrasts for diameter and confidence intervals.
summary(glm.screen.pos.abn.log)
lapply(coefficients(glm.screen.pos.abn.log)[1:6], function(x) x-coefficients(glm.screen.pos.abn.log)[2])
confint(glm.screen.pos.abn.log)
# Risks
tab3risk <- function(var, cat) {      # this function prints the prevalence as a percentage, then median/IQR risk
  df <- subset(dplyr::filter(data.screen.abn.pos, interval==2), select=c("post.risk.abn.pos", var))
  df <- df[df[,2] %in% cat,]
  print(round(100*(nrow(df)/nrow(dplyr::filter(data.screen.abn.pos, interval==2))), 3))
  100*round(quantile(df[,1], probs=c(0.25, 0.5, 0.75)), 5)}
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
100*nrow(dplyr::filter(data.screen.abn.pos, interval==2 & (any.right.mid==1 | any.lingula==1)))/nrow(filter(data.screen.abn.pos, interval==2))
100*with(dplyr::filter(data.screen.abn.pos, interval==2 & (any.right.mid==1 | any.lingula==1)), quantile(post.risk.abn.pos, probs=c(0.25, 0.5, 0.75)))
# Prevalence/probs for any poorly defined margins or any unable to characterize margins
100*nrow(dplyr::filter(data.screen.abn.pos, interval==2 & (any.poor.def==1 | any.margin.unab==1)))/nrow(filter(data.screen.abn.pos, interval==2))
100*with(dplyr::filter(data.screen.abn.pos, interval==2 & (any.poor.def==1 | any.margin.unab==1)), quantile(post.risk.abn.pos, probs=c(0.25, 0.5, 0.75)))

#### Properties of NLST screening, abnormal screens ####
# Is pre-screening risk important?
glm.screen.pos.nopsr <- glm(case ~ 1, data=data.screen.abn.pos, family=binomial(link='log'))
summary(glm.screen.pos.nopsr)
summary(glm.screen.pos)
lrtest(glm.screen.pos.nopsr, glm.screen.pos)
# Does the interval matter? no (p=0.13)
glm.screen.pos.by.int <-  glm(case ~ log1yrisk:as.factor(interval) -1, data=data.screen.abn.pos, family=binomial(link='log'))
summary(glm.screen.pos.by.int)
1-pchisq(glm.screen.pos$deviance - glm.screen.pos.by.int$deviance, length(glm.screen.pos.by.int$coefficients) - length(glm.screen.pos$coefficients))
with(data.screen.abn.pos, table(interval))
# Do previous screens matter? no (p=0.62)
glm.screen.pos.2levels <- glm(case ~ log1yrisk:screen.comb -1, data=dplyr::filter(data.screen.abn.pos, interval==2 & !is.na(screen.comb)), family=binomial(link='log'))
summary(glm.screen.pos.2levels)
confint(glm.screen.pos.2levels)
glm.screen.pos.1level <- glm(case ~ log1yrisk -1, data=dplyr::filter(data.screen.abn.pos, interval==2 & !is.na(screen.comb)), family=binomial(link='log'))
summary(glm.screen.pos.1level)
1-pchisq(glm.screen.pos.1level$deviance - glm.screen.pos.2levels$deviance, length(glm.screen.pos.2levels$coefficients) - length(glm.screen.pos.1level$coefficients))
with(filter(data.screen.abn.pos, interval==2), table(screen.comb))
# Residual effects of LCRAT variables examined below

#### Effects of abnormalities and nodules ####
# Backwards stepwise: selects any.lingula, any.mixed, any.other.att, any.right.mid, any.poor.def, any.upper, any.margin.unab, any.spiculation. diam.cat, any.growth
scr.pos.full <- glm(paste("case ~ logit1yrisk -1 +",paste(abnlist.pos.int, collapse="+"),sep=""), data=dplyr::filter(data.screen.abn.pos), family=binomial(link='logit'))
bsw.scr.pos <- step(scr.pos.full, direction="backward", scope = list(lower = case ~ logit1yrisk -1, upper = scr.pos.full))
    # Fit this model
mod.vars.bsw <- glm(case ~ logit1yrisk:any.lingula + logit1yrisk:any.mixed + logit1yrisk:any.right.mid + logit1yrisk:any.poor.def + logit1yrisk:any.upper +
            logit1yrisk:any.other.att + logit1yrisk:any.margin.unab + logit1yrisk:any.spiculation + logit1yrisk:diam.cat + logit1yrisk:any.growth -1, 
            data=data.screen.abn.pos, family=binomial(link='logit'))
summary(mod.vars.bsw)
  # get contrast amounts for diameter variable (for stepwise model in table)
lapply(coefficients(mod.vars.bsw)[9:14], function(x) x-coefficients(mod.vars.bsw)[10])
    # Get LRT p-value for diam.cat: p<0.000001
mod.vars.bsw.nodiamcat <- glm(case ~ logit1yrisk + logit1yrisk:any.lingula + logit1yrisk:any.mixed + logit1yrisk:any.right.mid + logit1yrisk:any.poor.def + logit1yrisk:any.upper +
             logit1yrisk:any.other.att + logit1yrisk:any.margin.unab + logit1yrisk:any.spiculation + logit1yrisk:any.growth -1, 
             data=data.screen.abn.pos, family=binomial(link='logit'))
lrtest(mod.vars.bsw.nodiamcat, mod.vars.bsw)
# Lasso with intermediate lambda - selects any.upper, any.spiculation, any.smooth, any.mixed, any.other.att, any.growth, 4 of 6 diam.cat categories
set.seed(61116)                                                        # why isn't this omitting the intercept? 
x = model.matrix(case ~ logit1yrisk + logit1yrisk:. -1, data = data.screen.abn.pos[,c("case","logit1yrisk",abnlist.pos)])
cv.lasso = cv.glmnet(x, data.screen.abn.pos$case, alpha=1, family="binomial")
out <- glmnet(x, data.screen.abn.pos$case, alpha=1, family="binomial")
predict(out, type="coefficients", s=(cv.lasso$lambda.min+cv.lasso$lambda.1se)/2)
    # Fit this model - log scale
mod.vars.lasso <- glm(case ~ log1yrisk:any.upper + log1yrisk:any.spiculation + log1yrisk:any.smooth + log1yrisk:any.mixed +
                        log1yrisk:any.other.att + log1yrisk:any.growth + log1yrisk:diam.cat -1, data=data.screen.abn.pos, family=binomial(link='log'))
summary(mod.vars.lasso)
    # get contrast amounts for diameter variable (for lasso model in table)
lapply(coefficients(mod.vars.lasso)[7:12], function(x) x-coefficients(mod.vars.lasso)[8])
    # Get LRT p-value for diam.cat: p<0.00001
mod.vars.lasso.nodiamcat <- glm(case ~ log1yrisk + log1yrisk:any.upper + log1yrisk:any.spiculation + log1yrisk:any.smooth + log1yrisk:any.mixed +
                        log1yrisk:any.other.att + log1yrisk:any.growth -1, data=data.screen.abn.pos, family=binomial(link='log'))
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
                      data=data.screen.abn.pos, family=binomial(link='log'))
summary(glm.screen.pos.abn.dev)
# Get LRT p-value for diam.cat
glm.screen.pos.abn.dev.nodiamcat <- glm(case ~ log1yrisk + log1yrisk:any.growth +
                                log1yrisk:any.upper + log1yrisk:I(any.right.mid==1|any.lingula==1) +
                                log1yrisk:I(any.new.nodule==T & diam.cat %in% c("4-5","6-7")) +
                                log1yrisk:any.mixed + log1yrisk:any.spiculation + log1yrisk:I(any.poor.def==1|any.margin.unab==1) -1, 
                              data=data.screen.abn.pos, family=binomial(link='log'))
lrtest(glm.screen.pos.abn.dev.nodiamcat, glm.screen.pos.abn.dev)
# Compare this model to one that ignores pre-screening risk
glm.screen.pos.abn.nopsr <-  glm(case ~ diam.cat + any.growth + any.upper + 
                      I(any.right.mid==1|any.lingula==1) + any.mixed + any.spiculation + 
                      I(any.new.nodule==T & diam.cat %in% c("4-5","6-7")) +
                      I(any.poor.def==1|any.margin.unab==1) -1, data=data.screen.abn.pos, family=binomial(link='log'))
summary(glm.screen.pos.abn.nopsr)
# Check for residual effects of LCRAT variables.
titles <- c("var","null model # param", "extended model # param", "LRT p-value", "check same # obs")
mat.out.ns <- matrix(rep(NA),nrow=length(varlist),ncol=length(titles))
mod.without <- glm(case ~ log1yrisk:diam.cat + log1yrisk:any.growth + log1yrisk:any.upper + log1yrisk:I(any.right.mid==1|any.lingula==1) + log1yrisk:any.mixed + log1yrisk:any.spiculation + log1yrisk:I(any.poor.def==1|any.margin.unab==1) + log1yrisk:I(any.new.nodule==T & diam.cat %in% c("4-5","6-7")) -1, data=data.screen.abn.pos, family=binomial(link='log'))
for (x in seq_along(varlist)) {
  mod.with <- glm(substitute(case ~ log1yrisk:diam.cat + log1yrisk:any.growth + log1yrisk:any.upper + log1yrisk:I(any.right.mid==1|any.lingula==1) + log1yrisk:any.mixed + log1yrisk:any.spiculation + log1yrisk:I(any.poor.def==1|any.margin.unab==1) + log1yrisk:i -1, list(i=as.name(varlist[x]))), data=data.screen.abn.pos, family=binomial(link='log'))
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
                id=pid, data=data.screen.abn.pos, family=binomial(link='log'), corstr="exchangeable", waves=interval))
summary(glm.screen.pos.abn.log)

### Calibration and validation analyses ###
# 10-fold cross-validated calibration
set.seed(61116)
data.screen.abn.pos$randgrp <- base::sample(1:10, nrow(data.screen.abn.pos), replace=T)
data.screen.abn.pos$cvpred <- NA
for (i in 1:10) {
  fit <- glm(formula = case ~ log1yrisk:diam.cat + log1yrisk:any.growth + log1yrisk:any.upper + 
               log1yrisk:I(any.right.mid==1|any.lingula==1) + log1yrisk:any.mixed + 
               log1yrisk:any.spiculation + log1yrisk:I(any.poor.def==1|any.margin.unab==1) -1, 
             family = binomial(link = "log"), data = dplyr::filter(data.screen.abn.pos, randgrp!=i)) 
  data.screen.abn.pos[data.screen.abn.pos$randgrp==i,]$cvpred <- predict(fit, newdata=data.screen.abn.pos[data.screen.abn.pos$randgrp==i,], type="response")
}
data.screen.abn.pos <- mutate(data.screen.abn.pos, cvpred.ntile = ntile(cvpred, 5))
data.screen.abn.pos %>% group_by(cvpred.ntile) %>% summarise(pred.cases= sum(cvpred), obs.cases = sum(case))
c(sum(data.screen.abn.pos$cvpred), sum(data.screen.abn.pos$case))
poisson.test(round(sum(data.screen.abn.pos$cvpred),0), sum(data.screen.abn.pos$case), alternative="two.sided") # p-value, requires rounding

### Calculate AUCs ###
## Optimism-corrected AUCs - have to use logistic models for this, and have to actually add the interaction terms to the dataset.
  # http://thestatsgeek.com/2014/10/04/adjusting-for-optimismoverfitting-in-measures-of-predictive-ability-using-bootstrapping/
# install.packages('rms')
library(rms)
data.screen.abn.pos <- mutate(data.screen.abn.pos, logit1yriskdiamcat0 = logit1yrisk*I(diam.cat=="0"), 
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
                    logit1yriskanypoordefanymarginunab -1, x=T, y=T, data=data.screen.abn.pos)
set.seed(61116)
validate(mod.ns.pos, B=1000)
c(0.5*(0.5938+1), 0.5*(0.5767+1)) # AUC = 0.5(Dxy+1). 
  # the above gives naive (index.orig) and optimism-corrected (index.corrected) AUCs. 0.79 is OC-AUC



#### Figure 1 - risk density #### 

prq <- with(dplyr::filter(data.screen.abn.pos, interval==2), quantile(post.risk.abn.pos, probs=c(0.5, 0.75, 0.95)))
prq_ar_end <- c(prq[1]-0.003, prq[2]-0.003, prq[3]+0.003)
prq

png(file=here("abnormals_density.png"),width=1200,height=850)
ggplot(data=dplyr::filter(data.screen.abn.pos, interval==2)) + theme_bw() +
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
with(dplyr::filter(data.screen.abn.pos, interval==2), quantile(post.risk.abn.pos, probs=seq(0,1,0.1)))
with(dplyr::filter(data.screen.abn.pos, interval==2), CrossTable(I(post.risk.abn.pos>0.002), I(post.risk.abn.pos<0.03)))

# Look at the 5% with highest risk (>7.1%)
View(dplyr::filter(data.screen.abn.pos, interval==2 & post.risk.abn.pos>0.071))
with(dplyr::filter(data.screen.abn.pos, interval==2 & post.risk.abn.pos>0.071), CrossTable(diam.cat))
with(dplyr::filter(data.screen.abn.pos, interval==2 & post.risk.abn.pos>0.071), CrossTable(any.upper))
with(dplyr::filter(data.screen.abn.pos, interval==2 & post.risk.abn.pos>0.071), CrossTable(I(any.right.mid == 1 | any.lingula == 1)))
with(dplyr::filter(data.screen.abn.pos, interval==2 & post.risk.abn.pos>0.071), CrossTable(any.mixed))
with(dplyr::filter(data.screen.abn.pos, interval==2 & post.risk.abn.pos>0.071), CrossTable(any.spiculation))
with(dplyr::filter(data.screen.abn.pos, interval==2 & post.risk.abn.pos>0.071), CrossTable(I(any.poor.def == 1 | any.margin.unab == 1)))
with(dplyr::filter(data.screen.abn.pos, interval==2 & post.risk.abn.pos>0.071), CrossTable(any.growth))

# Look at the 8% with a risk <0.20% (e.g. potential extended interval)
View(dplyr::filter(data.screen.abn.pos, interval==2 & post.risk.abn.pos<0.002))
      # Among these 8% look at 5-year pre-screening risk (1.9% threshold is proposed by LCRAT)
with(dplyr::filter(data.screen.abn.pos, interval==2 & post.risk.abn.pos<0.002), quantile(I(5*prescr.1yrisk), probs=seq(0,1,.1)))
with(dplyr::filter(data.screen.abn.pos, interval==2), quantile(I(5*prescr.1yrisk), probs=seq(0,1,.1)))
with(dplyr::filter(data.screen.abn.pos, interval==2 & post.risk.abn.pos<0.002), CrossTable(diam.cat))

