

#### Setup ####

library(here)
# Read in my Kovalchik prediction function
source("kovalchik.R")


# install coxph.risk tar file
# install lcmodels package
devtools::install_github('marskar/coxph_risk')
devtools::install_github('marskar/lcmodels')

# install.packages('pROC')
packages <- c("dplyr","ggplot2","survival","gmodels","coxph.risk","geepack","MESS","psych","Hmisc","glmnet","boot")
lapply(packages, require, c = T)
# Load PLCO and NLST data sets. New versions provided by Li Cheung on 11 July 2016, 14 July 2016, 20 July 2016, 8 Aug 2016.
load("hilary.RData")  # Load NLST and PLCO data
# In the PLCO dataset, impute missing family history values to 0
plco$fam.lung.trend <- ifelse(is.na(plco$fam.lung.trend), 0, plco$fam.lung.trend)
plco.control <- subset(plco, control.group==1) # control arm of PLCO who had no chest xray
# Remove people with <30 pack-years and age<55 or age>74 from NLST
nlst <- filter(nlst, pkyears.cat!="[0,30)" & age>=55 & age<=74)
# Make a new pack-years variable to get rid of the [0,30) level
nlst <- mutate(nlst, pkyears.cat.clone=ifelse(pkyears.cat=="[30,40)", "[30,40)", ifelse(pkyears.cat=="[40,50)", "[40,50)",
                      ifelse(pkyears.cat=="[50,Inf]", "[50,Inf]", NA))))  # I checked this with a table
nlst$pkyears.cat <- as.factor(nlst$pkyears.cat.clone)
# Make a variable for days to diagnosis
nlst$days_to_dx <- ifelse(nlst$case==1, 365*nlst$incidence.years, NA)
# Make a subset of NLST data with the LCRAT variables that we will need later to merge back with at-risk datasets
varlist <- c("female","race","edu6","fam.lung.trend","emp","bmi","cpd","pkyears.cat","age","qtyears","smkyears")
nlst.sub <- as.data.frame(cbind(nlst[,varlist], pid=nlst$pid, lss=as.numeric(nlst$lss)))

# To later calculate pre-screening risk, we must first fit the incidence model and other-cause death models in PLCO.
LCRAT <- coxph(Surv(incidence.years, case) ~ 
                 female+race+edu6+fam.lung.trend+emp+I(bmi<=18.5)+I(cpd>20)+as.factor(pkyears.cat)+
                 I(log(age))+I(log(bmi))+I(log(qtyears+1))+smkyears,data = plco.control)
cox.death <- coxph(Surv(years.followed, other.cause.death) ~ 
                     female+race+edu6+emp+I(bmi <= 18.5)+I(cpd>20)+as.factor(pkyears.cat)+I((age)^2)+I((bmi-25)^2)+
                     I(log(qtyears+1))+smkyears, data = plco.control)  

# Subset to CT arm in NLST and make a pos/neg variable for the first, second, and third screens
nlst$T0posneg <- ifelse(nlst$truefalse_scrnres_ly0 %in% c(4,5,6), 0, NA)
nlst$T0posneg <- ifelse(nlst$truefalse_scrnres_ly0 %in% c(1,2,3), 1, nlst$T0posneg)
nlst$T1posneg <- ifelse(nlst$truefalse_scrnres_ly1 %in% c(4,5,6), 0, NA)
nlst$T1posneg <- ifelse(nlst$truefalse_scrnres_ly1 %in% c(1,2,3), 1, nlst$T1posneg)
nlst$T2posneg <- ifelse(nlst$truefalse_scrnres_ly2 %in% c(4,5,6), 0, NA)
nlst$T2posneg <- ifelse(nlst$truefalse_scrnres_ly2 %in% c(1,2,3), 1, nlst$T2posneg)
nlst$prescr.1yrisk.T0 <- risk.kovalchik(0, 1, nlst, LCRAT, cox.death)  # add 1y risk to NLST dataset for descriptive stats

# Subset to CT arm and create screening history variables
nlst.CT <- subset(nlst, screen_group=="CT")
nlst.CT <- mutate(nlst.CT, hist.T0.T1 = 1*(T0posneg==0 & T1posneg==0) + 2*(T0posneg==0 & T1posneg==1) + 3*(T0posneg==1 & T1posneg==0) + 4*(T0posneg==1 & T1posneg==1))
nlst.CT$hist.T0.T1 <- factor(nlst.CT$hist.T0.T1, levels=c(1,2,3,4), labels=c("Neg-Neg","Neg-Pos","Pos-Neg","Pos-Pos"))
nlst.CT <- mutate(nlst.CT, hist.T1.T2 = 1*(T1posneg==0 & T2posneg==0) + 2*(T1posneg==0 & T2posneg==1) + 3*(T1posneg==1 & T2posneg==0) + 4*(T1posneg==1 & T2posneg==1))
nlst.CT$hist.T1.T2 <- factor(nlst.CT$hist.T1.T2, levels=c(1,2,3,4), labels=c("Neg-Neg","Neg-Pos","Pos-Neg","Pos-Pos"))

### Create datasets for risk from T0 to T1. One for interval, one for screen-detected.
    ## .neg intended for analysis of interval cancer risk; .scrisk intended for analysis of screen-detected cancers
  # 0 inadeq, 1 true-pos, 2 poss true-pos, 3 false-pos, 4 true-neg, 5 poss false-neg, 6 false-neg
# At risk in interval: T0 negatives only. Case status: false-negative at T0.
nlst.CT.T1.neg <- filter(nlst.CT, truefalse_scrnres_ly0 %in% c(4,5,6))  
nlst.CT.T1.neg$case_T1_interval <- ifelse(nlst.CT.T1.neg$truefalse_scrnres_ly0==6, 1, 0)
# At risk for screen-detected at T1: either false-positive or true-negative at T0, and did not have any of the following at T1:
  # inadequate image, left study, refused, wrong screen, erroneous report of LC, form not submitted (no missing values of scr_res0).
  # Case status: case=1 AND either of (true-pos at T1 or T1 is coded as "not expected: cancer/death in screening window")
nlst.CT.T1.scrisk <- filter(nlst.CT, truefalse_scrnres_ly0 %in% c(2,3,4,5) & scr_res1 %!in% c(10,11,15,17,95,97))
nlst.CT.T1.scrisk$case_T1_screen <- ifelse(nlst.CT.T1.scrisk$case==1 & 
                          (nlst.CT.T1.scrisk$truefalse_scrnres_ly1==1 | nlst.CT.T1.scrisk$scr_res1 %in% c(23,24)), 1, 0)
### Create datasets for risk from T1 to T2. One for interval, one for screen-detected.
# At risk in interval: T1 negatives only. Case status: false-negative at T1.
nlst.CT.T2.neg <- filter(nlst.CT, truefalse_scrnres_ly1 %in% c(4,5,6))  
nlst.CT.T2.neg$case_T2_interval <- ifelse(nlst.CT.T2.neg$truefalse_scrnres_ly1==6, 1, 0)
# At risk for screen-detected at T2: either false-positive or true-negative at T1, and did not have any of the following at T2:
  # inadequate image, left study, refused, wrong screen, erroneous report of LC, form not submitted (no missing values of scr_res0).
  # Case status: case=1 AND either of (true-pos at T2 or T2 is coded as "not expected: cancer/death in screening window")
nlst.CT.T2.scrisk <- filter(nlst.CT, truefalse_scrnres_ly1 %in% c(2,3,4,5) & scr_res2 %!in% c(10,11,15,17,95,97))
nlst.CT.T2.scrisk$case_T2_screen <- ifelse(nlst.CT.T2.scrisk$case==1 &
                          (nlst.CT.T2.scrisk$truefalse_scrnres_ly2==1 | nlst.CT.T2.scrisk$scr_res2 %in% c(23,24)), 1, 0)
### Create a dataset for risk during "interval" after T2 (within 1 year)
nlst.CT.T3.neg <- filter(nlst.CT, truefalse_scrnres_ly2 %in% c(4,5,6))
nlst.CT.T3.neg$case_T3_interval <- ifelse(nlst.CT.T3.neg$truefalse_scrnres_ly2==6, 1, 0)

### Construct dataset to model risk of ALL interval cancers across all 3 screens
  # Interval variable in data.interval datasets is 1 for T0-T1, 2 for T1-T2, and 3 for post-T2 intervals
data.interval <- data.frame(pid=c(nlst.CT.T1.neg$pid, nlst.CT.T2.neg$pid, nlst.CT.T3.neg$pid),
                            case=c(nlst.CT.T1.neg$case_T1_interval, nlst.CT.T2.neg$case_T2_interval, nlst.CT.T3.neg$case_T3_interval),
                            interval=c(rep(1,times=nrow(nlst.CT.T1.neg)), rep(2, times=nrow(nlst.CT.T2.neg)), rep(3, times=nrow(nlst.CT.T3.neg))))
# Merge this back with covariates from NLST. Add screening history variable.
data.interval <- merge(data.interval, nlst.sub, by="pid", all.x=T)
data.interval <- merge(data.interval, select(nlst.CT, pid, hist.T0.T1, hist.T1.T2), by="pid", all.x=T, all.y=F)  ## delete??
data.interval$screen.hist <- ifelse(data.interval$interval==2, data.interval$hist.T0.T1, ifelse(data.interval$interval==3, data.interval$hist.T1.T2, NA))
data.interval$screen.hist <- factor(data.interval$screen.hist, levels=c(1,2,3,4), labels=c("Neg-Neg","Neg-Pos","Pos-Neg","Pos-Pos"))  
data.interval <- select(data.interval, -c(hist.T0.T1,hist.T1.T2))                        ## delete??
data.interval <- arrange(data.interval, pid, interval)
# Update age, quit-years, and smoke-years by adding a year for T1 and T2
data.interval <- mutate(data.interval, age=ifelse(interval==2, age+1, ifelse(interval==3, age+2, age)),
                smkyears=ifelse(interval==2 & qtyears==0, smkyears+1, ifelse(interval==3 & qtyears==0, smkyears+2, smkyears)),
                qtyears=ifelse(interval==2 & qtyears>0, qtyears+1, ifelse(interval==3 & qtyears>0, qtyears+2, qtyears)))
data.interval <- mutate(data.interval, pkyears.cont=cpd*smkyears/20)  # using new smoke-years, update pack-years, then re-categorize
data.interval <- mutate(data.interval, pkyears.cat=as.factor(ifelse(pkyears.cont>=30 & pkyears.cont<40, "[30,40)",
                  ifelse(pkyears.cont>=40 & pkyears.cont<50, "[40,50)", ifelse(pkyears.cont>=50 & pkyears.cont<999,"[50,Inf]",NA)))))

# Construct dataset to model risk of ALL screen-detected cancers (at T1 and T2)
  # Interval variable in data.screen datasets is 1 for risk at T1 and 2 for risk at T2
data.screen <- data.frame(pid=c(nlst.CT.T1.scrisk$pid, nlst.CT.T2.scrisk$pid),
                          case=c(nlst.CT.T1.scrisk$case_T1_screen, nlst.CT.T2.scrisk$case_T2_screen),
                          screen.result=c(nlst.CT.T1.scrisk$T0posneg, nlst.CT.T2.scrisk$T1posneg),
                          interval=c(rep(1,times=nrow(nlst.CT.T1.scrisk)), rep(2, times=nrow(nlst.CT.T2.scrisk))))
# Merge this back with covariates from NLST
data.screen <- merge(data.screen, nlst.sub, by="pid", all.x=T)
# Add a variable for lagged screen result & a 6-level variable for all combinations
data.screen <- data.screen %>% group_by(pid) %>% mutate(lag.screen = lag(screen.result, order_by=interval))
data.screen <- mutate(data.screen, screen.comb = 1*(interval==1 & screen.result==0) +
          2*(interval==1 & screen.result==1) + 3*(interval==2 & lag.screen==0 & screen.result==0) +
          4*(interval==2 & lag.screen==0 & screen.result==1) + 5*(interval==2 & lag.screen==1 & screen.result==0) +
          6*(interval==2 & lag.screen==1 & screen.result==1))
data.screen$screen.comb <- factor(data.screen$screen.comb, levels = c(1,2,3,4,5,6),
                                  labels = c("Neg","Pos","Neg-Neg","Neg-Pos","Pos-Neg","Pos-Pos"))
# Update age, quit-years, and smoke-years by adding a year for T1
data.screen <- mutate(data.screen, age=as.numeric(age), smkyears=as.numeric(smkyears), qtyears=as.numeric(qtyears))
data.screen <- mutate(data.screen, age=ifelse(interval==2, age+1, age),
                      smkyears=ifelse(interval==2 & qtyears==0, smkyears+1, smkyears),
                      qtyears=ifelse(interval==2 & qtyears>0, qtyears+1, qtyears))
data.screen <- mutate(data.screen, pkyears.cont=cpd*smkyears/20) # using new smoke-years, update pack-years, then re-categorize
data.screen <- mutate(data.screen, pkyears.cat=as.factor(ifelse(pkyears.cont>=30 & pkyears.cont<40, "[30,40)",
                      ifelse(pkyears.cont>=40 & pkyears.cont<50, "[40,50)", ifelse(pkyears.cont>=50 & pkyears.cont<999,"[50,Inf]",NA)))))

# Load abnormalities data (person-screen level) and merge with data.screen
  # This dataset was generated by the program prepare_abn_data_vX.R
load(here("abn.spl.20160810.rdata"))
data.screen.abn <- merge(data.screen, abn.pl.all, by=c("pid","interval"), all.x=T, all.y=F)
# Replace NAs with 0 (not present) for appropriate variables
replacevars  <- names(abn.pl.all)[!names(abn.pl.all) %in% c("pid","interval","LRcat","LRcatcol.neg","LRcatcol.pos")]
data.screen.abn[replacevars][is.na(data.screen.abn[replacevars])]  <- 0
# Make a variable for including observations in Lung-RADS analysis
data.screen.abn$LR.include <- (data.screen.abn$LRcat %in% c("1","2","3","4A","4B","4X"))
# Merge abnormalities data (person-screen level) with data.interval
data.interval.abn <- merge(data.interval, abn.pl.all, by=c("pid","interval"), all.x=T, all.y=F)
data.interval.abn[replacevars][is.na(data.interval.abn[replacevars])]  <- 0
# Will need this vector for exploratory analysis of abnormalities (24 Jan 2017 - also create interaction vectors)
  # abnlist.neg contains a list of abnormalities variables that are relevant to negative CTs. abnlist.pos is the list for positive CTs.
abnlist <- abnlist.pos <- names(abn.pl.all)[3:32]
abnlist.neg <- abnlist[c(4,5,6,7,8,9,10,11,12,13,14)]
abnlist.neg.int <- lapply(abnlist.neg, function(x) {substitute(logit1yrisk:i, list(i=as.name(x)))}) # these lists create interaction terms
abnlist.pos.int <- lapply(abnlist.pos, function(x) {substitute(logit1yrisk:i, list(i=as.name(x)))})
# Create variable for log(diameter)
data.screen.abn$log.diam <- log(data.screen.abn$longest.diam+1)
# Calculate pre-screening risk inside this dataset
data.interval.abn$prescr.1yrisk <- risk.kovalchik(0, 1, data.interval, LCRAT, cox.death)
data.interval.abn <- mutate(data.interval.abn, log1yrisk=log(prescr.1yrisk), logit1yrisk=log(prescr.1yrisk/(1-prescr.1yrisk)))
data.screen.abn$prescr.1yrisk <- risk.kovalchik(0, 1, data.screen.abn, LCRAT, cox.death)
data.screen.abn <- mutate(data.screen.abn, log1yrisk=log(prescr.1yrisk), logit1yrisk=log(prescr.1yrisk/(1-prescr.1yrisk)))
# These datasets are needed to separately model screen-detected cancers incorporating abnormalities for negatives and false-positives
data.screen.abn.neg <- filter(data.screen.abn, screen.result==0)
data.screen.abn.pos <- filter(data.screen.abn, screen.result==1)
# Make a categorical variable for diameter
data.screen.abn.pos <- mutate(data.screen.abn.pos, diam.cat = 1*(longest.diam==0)+2*(longest.diam>0 & longest.diam<=5)+
                        3*(longest.diam>5 & longest.diam<=7) + 4*(longest.diam>7 & longest.diam<=10) +
                          5*(longest.diam>10 & longest.diam<=13) + 6*(longest.diam>13 & longest.diam<100))
data.screen.abn.pos$diam.cat <- factor(data.screen.abn.pos$diam.cat, levels=c(1:6),labels=c("0","4-5","6-7","8-10","11-13","14+"))
# Create a variable for any.growth that reflects a group in which growth can't be assessed (i.e. screen=T0)
data.screen.abn.pos <- mutate(data.screen.abn.pos, growth.3l = 
                          1*(interval==1) + 2*(interval==2 & any.growth==0) + 3*(interval==2 & any.growth==1))
data.screen.abn.pos$growth.3l <- factor(data.screen.abn.pos$growth.3l, levels = c(1,2,3), labels=c("NA","No","Yes"))

# Make dataset of unique individuals for descriptive table of screen-negatives
all.subj.neg <- filter(nlst.CT, pid %in% data.interval.abn$pid | pid %in% data.screen.abn.neg$pid)
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

# Run the main models
# Overall effects (without specific CT features)
glm.interval <- glm(case ~ log1yrisk -1, data=data.interval.abn, family=binomial(link='log'))
data.interval.abn$post.risk.interv <- fitted.values(glm.interval)
glm.screen.neg <- glm(case ~ log1yrisk -1, data=data.screen.abn.neg, family=binomial(link='log'))
data.screen.abn.neg$post.risk.neg.overall <- fitted.values(glm.screen.neg)
# With specific CT findings
glm.int.abn <- glm(case ~ log1yrisk + log1yrisk:adenop.consol -1, data=data.interval.abn, family=binomial(link='log'))
data.interval.abn$post.risk.abn <- fitted.values(glm.int.abn)
glm.screen.neg.abn <- glm(case ~ log1yrisk + log1yrisk:consolidation + log1yrisk:emphysema -1, data=data.screen.abn.neg, family=binomial(link='log'), na.action=na.exclude)
data.screen.abn.neg$post.risk.abn <- fitted.values(glm.screen.neg.abn)



# --------------------------- run code to this line for data setup ----------------------------- #


#### Descriptive stats ####

# Descriptive table 1 for analysis after NEGATIVE screen
nrow(all.subj.neg)  # number of unique individuals
length(unique(filter(nlst.CT, T0posneg==0 | T1posneg==0 | T2posneg==0)$pid)) # confirm - this is the same # with at least one negative screen
length(unique(data.interval$pid)) # number of unique individuals in interval cancers analysis
CrossTable((data.interval %>% group_by(pid) %>% summarise(times.in.interv.analysis=n()))$times.in.interv.analysis) # times included in interval ca analysis
length(unique(data.screen.abn.neg$pid)) # number of unique individuals in next-screen analysis
CrossTable((data.screen.abn.neg %>% group_by(pid) %>% summarise(times.in.screen.analysis=n()))$times.in.screen.analysis) # times included in next-screen analysis
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
with(filter(data.screen.abn.neg, interval==2), CrossTable(screen.comb, case))   # numbers considered for Markov assumption test
c(sum(data.interval.abn$case), nrow(data.interval.abn), sum(data.interval.abn$case)/nrow(data.interval.abn)) # interval cancers: cases, # at risk, overall risk
c(sum(data.screen.abn.neg$case), nrow(data.screen.abn.neg), sum(data.screen.abn.neg$case)/nrow(data.screen.abn.neg)) # next-screen cancers after negative: cases, # at risk, overall risk
range_without_outliers(all.subj.neg$prescr.1yrisk)


#### Model development: Interval cancer among NEGATIVES ####

# Interval cancers: overall model (no abnormalities) - properties of screening
# Confirm that pre-screening risk improves the model
int.nopsr <- glm(case ~ 1, data=data.interval.abn, family=binomial(link='log'))
summary(int.nopsr)
int.psr <- glm(case ~ log1yrisk+1, data=data.interval.abn, family=binomial(link='log'))
summary(int.psr)
1-pchisq(int.nopsr$deviance - int.psr$deviance, length(int.psr$coefficients)-length(int.nopsr$coefficients))  # LRT
# Overall model results
  # glm.interval.abn <- glm(case ~ log1yrisk -1, data=data.interval.abn, family=binomial(link='log'))   # run above in setup
  # data.interval.abn$post.risk.interv <- fitted.values(glm.interval)                                   # run above in setup
summary(glm.interval)
confint(glm.interval)
# Does the risk coefficient differ by interval? No (LRT p=0.23). Steps below: fit model, estimate 3 exponents, get p-value, get counts
glm.interval.intervals <- glm(case ~ log1yrisk + log1yrisk:I(as.numeric(interval==2)) + log1yrisk:I(as.numeric(interval==3)) -1, data=data.interval.abn, family=binomial(link='log'))
c(coefficients(glm.interval.intervals)[1], coefficients(glm.interval.intervals)[1]+coefficients(glm.interval.intervals)[2], coefficients(glm.interval.intervals)[1]+coefficients(glm.interval.intervals)[3])
1-pchisq(glm.interval$deviance - glm.interval.intervals$deviance, length(glm.interval.intervals$coefficients)-length(glm.interval$coefficients))
with(data.interval.abn, table(interval))
# Do previous screens matter? No (LRT p=0.99)
glm.int.2levels <- glm(case ~ log1yrisk:as.factor(screen.hist) -1, data=filter(data.interval.abn, interval %in% c(2,3)), family=binomial(link='log'))
summary(glm.int.2levels)
confint(glm.int.2levels)
glm.int.1level <- glm(case ~ log1yrisk -1, data=filter(data.interval.abn, interval %in% c(2,3) & !is.na(screen.hist)), family=binomial(link='log'))
summary(glm.int.1level)
1-pchisq(glm.int.1level$deviance-glm.int.2levels$deviance, df=length(glm.int.2levels$coefficients-length(glm.int.1level$coefficients)))

# Interval cancers: effects of abnormalities
  # Following a negative screen, the relevant CT features are in abnlist.neg. The relevant p-value is for the interaction (i.e. risk differs between the 0 and 1 levels)
# Backwards stepwise selection: selects other.above, benign.nodule, consolidation, adenopathy
int.full <- glm(paste("case ~ logit1yrisk -1 +",paste(abnlist.neg.int, collapse="+"),sep=""), data=data.interval.abn, family=binomial(link='logit'))
bsw.int <- step(int.full, direction="backward", scope = list(lower = case ~ logit1yrisk -1, upper = int.full))
  # Look at a model including these 4 effects
summary(glm(case ~ log1yrisk + log1yrisk:other.above + log1yrisk:benign.nodule + log1yrisk:consolidation + log1yrisk:adenopathy -1, data=data.interval.abn, family=binomial(link='log')))
# Lasso using intermediate lambda: selects adenopathy and consolidation
set.seed(61116)
x  <- model.matrix(case ~ logit1yrisk -1 + logit1yrisk:., data = data.interval.abn[,c("case","logit1yrisk",abnlist.neg)])
cv.lasso <- cv.glmnet(x, data.interval.abn$case, alpha=1, family="binomial")
out <- glmnet(x, data.interval.abn$case, alpha=1, family="binomial")
predict(out, type="coefficients", s=(cv.lasso$lambda.min+cv.lasso$lambda.1se)/2)
  # Look at a model including these two effects
summary(glm(case ~ log1yrisk + log1yrisk:consolidation + log1yrisk:adenopathy -1, data=data.interval.abn, family=binomial(link='log')))
# Based on discussion with Chris: include adenopathy and consolidation. Model as 1 variable (effect size is the same)
# Switch back to log scale for final model for interpretability (these models are run above in data setup section)
    # glm.int.abn.log <- glm(case ~ log1yrisk + log1yrisk:adenop.consol -1, data=data.interval.abn, family=binomial(link='log'))
    # data.interval.abn$post.risk.abn <- fitted.values(glm.int.abn.log)
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
  mod.without <- glm(case ~ logit1yrisk + logit1yrisk:adenop.consol -1, data=data.interval.abn, family=binomial(link='logit'))
  mod.with <- glm(substitute(case ~ logit1yrisk + logit1yrisk:adenop.consol + logit1yrisk:i -1, list(i=as.name(varlist[x]))), data=data.interval.abn, family=binomial(link='logit'))
  print(summary(mod.with))
  mat.out.interv[x,] <- c(varlist[x], length(mod.without$coefficients), sum(!is.na(mod.with$coefficients)), 1-pchisq(mod.without$deviance-mod.with$deviance, df=sum(!is.na(mod.with$coefficients))-length(mod.without$coefficients)), I(length(mod.without$residuals)==length(mod.with$residuals)))
}
rbind(titles, mat.out.interv)


#### Model development: Next-screen cancer among NEGATIVES ####

# Overall model for next-screen cancer among negatives (no abnormalities) - properties of screening
# Confirm that pre-screening risk improves the model
ns.nopsr <- glm(case ~ 1, data=data.screen.abn.neg, family=binomial(link='log'))
summary(ns.nopsr)
ns.psr <- glm(case ~ log1yrisk+1, data=data.screen.abn.neg, family=binomial(link='log'))
summary(ns.psr)
1-pchisq(ns.nopsr$deviance - ns.psr$deviance, length(ns.psr$coefficients)-length(ns.nopsr$coefficients))
# Overall model results
  # glm.screen.neg <- glm(case ~ log1yrisk -1, data=data.screen.abn.neg, family=binomial(link='log'))  # this is run above the line
  # data.screen.abn.neg$post.risk.neg.overall <- fitted.values(glm.screen.neg)
summary(glm.screen.neg)
confint(glm.screen.neg)
# Does the interval matter? no (p=0.38)
glm.screen.neg.by.int <-  glm(case ~ log1yrisk:as.factor(interval) -1, data=data.screen.abn.neg, family=binomial(link='log'))
summary(glm.screen.neg.by.int)
1-pchisq(glm.screen.neg$deviance - glm.screen.neg.by.int$deviance, length(glm.screen.neg.by.int$coefficients) - length(glm.screen.neg$coefficients))
with(data.screen.abn.neg, table(interval))
# Do previous screens matter? no (p=0.26)
glm.screen.neg.2levels <- glm(case ~ log1yrisk:screen.comb -1, data=filter(data.screen.abn.neg, interval==2 & !is.na(screen.comb)), family=binomial(link='log'))
summary(glm.screen.neg.2levels)
confint(glm.screen.neg.2levels)
glm.screen.neg.1level <- glm(case ~ log1yrisk -1, data=filter(data.screen.abn.neg, interval==2 & !is.na(screen.comb)), family=binomial(link='log'))
summary(glm.screen.neg.1level)
1-pchisq(glm.screen.neg.1level$deviance - glm.screen.neg.2levels$deviance, length(glm.screen.neg.2levels$coefficients) - length(glm.screen.neg.1level$coefficients))
# Do previous screens matter if we ignore pre-screening risk? p=0.14
glm.screen.neg.2levels.nopsr <- glm(case ~ screen.comb -1, data=filter(data.screen.abn.neg, interval==2 & !is.na(screen.comb)), family=binomial(link='log'))
summary(glm.screen.neg.2levels.nopsr)
glm.screen.neg.1level.nopsr <- glm(case ~ 1, data=filter(data.screen.abn.neg, interval==2 & !is.na(screen.comb)), family=binomial(link='log'))
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
scr.neg.full <- glm(paste("case ~ logit1yrisk -1 +",paste(abnlist.neg.int, collapse="+"),sep=""), data=data.screen.abn.neg, family=binomial(link='logit'))
bsw.scr.neg <- step(scr.neg.full, direction="backward", scope = list(lower = case ~ logit1yrisk -1, upper = scr.neg.full))
  # Look at a model including these 4 effects
summary(glm(case ~ log1yrisk + log1yrisk:opac.fibr + log1yrisk:nod6.not.susp + log1yrisk:consolidation + log1yrisk:emphysema -1, data=data.screen.abn.neg, family=binomial(link='log')))
# Lasso using intermediate lambda: selects ONLY logit1yrisk
set.seed(61116)
x  <-  model.matrix(case ~ logit1yrisk -1 + logit1yrisk:. , data = data.screen.abn.neg[,c("case","logit1yrisk",abnlist.neg)])
cv.lasso <- cv.glmnet(x, data.screen.abn.neg$case, alpha=1, family="binomial")
out <- glmnet(x, data.screen.abn.neg$case, alpha=1, family="binomial")
predict(out, type="coefficients", s=(cv.lasso$lambda.min+cv.lasso$lambda.1se)/2)
# From discussion with Chris: keep consolidation and emphysema.
# Switch back to log scale for final model for interpretability (this model is run above in data setup)
    # glm.screen.neg.abn <- glm(case ~ log1yrisk + log1yrisk:consolidation + log1yrisk:emphysema -1, data=data.screen.abn.neg, family=binomial(link='log'))
    # data.screen.abn.neg$post.risk.abn <- fitted.values(glm.screen.neg.abn)
summary(glm.screen.neg.abn)
      # Get estimates and CIs for the exponents
mat <- c(1,0,0) # Use this matrix for "neither noted"
mat <- c(1,1,0) # Use this matrix for consolidation
mat <- c(1,0,1) # Use this matrix for emphysema
stder <- sqrt(c(t(mat) %*% vcov(glm.screen.neg.abn) %*% mat))
c(coefficients(glm.screen.neg.abn) %*% mat, (coefficients(glm.screen.neg.abn) %*% mat)-1.96*stder, (coefficients(glm.screen.neg.abn) %*% mat)+1.96*stder)
# Check for residual effects of LCRAT variables. All p>0.05
titles <- c("var","null model # param", "extended model # param", "LRT p-value", "check same # obs")
mat.out.ns <- matrix(rep(NA),nrow=length(varlist),ncol=length(titles))
for (x in seq_along(varlist)) {
  mod.without <- glm(case ~ logit1yrisk + logit1yrisk:consolidation + logit1yrisk:emphysema -1, data=data.screen.abn.neg, family=binomial(link='logit'))
  mod.with <- glm(substitute(case ~ logit1yrisk + logit1yrisk:consolidation + logit1yrisk:emphysema + logit1yrisk:i -1, list(i=as.name(varlist[x]))), data=data.screen.abn.neg, family=binomial(link='logit'))
  mat.out.ns[x,] <- c(varlist[x], length(mod.without$coefficients), sum(!is.na(mod.with$coefficients)), 1-pchisq(mod.without$deviance-mod.with$deviance, df=sum(!is.na(mod.with$coefficients))-length(mod.without$coefficients)), I(length(mod.without$residuals)==length(mod.with$residuals)))
}
rbind(titles, mat.out.ns)
# What if we account for screening history along with pre-screening risk, emphysema, and consolidation? p=0.34
m1 <- glm(case ~ log1yrisk + log1yrisk:consolidation + log1yrisk:emphysema -1, data=filter(data.screen.abn.neg, interval==2 & !is.na(screen.comb)), family=binomial(link='log'))
summary(m1)
m2 <- glm(case ~ log1yrisk:screen.comb + log1yrisk:consolidation + log1yrisk:emphysema -1, data=filter(data.screen.abn.neg, interval==2 & !is.na(screen.comb)), family=binomial(link='log'))
summary(m2)
1-pchisq(m1$deviance - m2$deviance, length(m2$coefficients) - length(m1$coefficients))






#### Additional analyses (GEE, AUCs, cross-validation, etc) ####

### Comparison with GEE - this impacts the SEs negligibly ###
# Interval cancer models
summary(geeglm(case ~ log1yrisk -1, id=pid, data=data.interval, family=binomial(link='log'), corstr="exchangeable", waves=interval))
summary(glm.interval)
summary(geeglm(case ~ log1yrisk + log1yrisk:adenop.consol -1, id=pid, data=data.interval.abn, family=binomial(link='log'), corstr="exchangeable", waves=interval))
summary(glm.int.abn.log)
# Next-screen among negatives model
summary(geeglm(case ~ log1yrisk -1, id=pid, data=data.screen.abn.neg, family=binomial(link='log'), corstr="exchangeable", waves=interval))
summary(glm.screen.neg)
summary(geeglm(case ~ log1yrisk + log1yrisk:consolidation + log1yrisk:emphysema -1, id=pid, data=data.screen.abn.neg, family=binomial(link='log'), corstr="exchangeable", waves=interval))
summary(glm.screen.neg.abn)


### Calibration and validation analyses ###
# 10-fold cross-validated calibration
# Interval cancers
set.seed(61116)
data.interval.abn$randgrp <- base::sample(1:10, nrow(data.interval.abn), replace=T)
data.interval.abn$cvpred <- NA
for (i in 1:10) {
  fit <- glm(formula = case ~ log1yrisk + log1yrisk:adenop.consol - 1, 
             family = binomial(link = "log"), data = filter(data.interval.abn, randgrp!=i))
  data.interval.abn[data.interval.abn$randgrp==i,]$cvpred <- predict(fit, newdata=data.interval.abn[data.interval.abn$randgrp==i,], type="response")
}
data.interval.abn <- mutate(data.interval.abn, cvpred.ntile = ntile(cvpred, 5))
data.interval.abn %>% group_by(cvpred.ntile) %>% summarise(pred.cases= sum(cvpred), obs.cases = sum(case))
c(sum(data.interval.abn$cvpred), sum(data.interval.abn$case)) # number obs and expected cases
poisson.test(round(sum(data.interval.abn$cvpred),0), sum(data.interval.abn$case), alternative="two.sided") # p-value, requires rounding

# Next-screen among screen-negatives
set.seed(61116)
data.screen.abn.neg$randgrp <- base::sample(1:10, nrow(data.screen.abn.neg), replace=T)
data.screen.abn.neg$cvpred <- NA
for (i in 1:10) {
  fit <- glm(formula = case ~ log1yrisk + log1yrisk:consolidation + log1yrisk:emphysema - 1, 
             family = binomial(link = "log"), data = filter(data.screen.abn.neg, randgrp!=i))
  data.screen.abn.neg[data.screen.abn.neg$randgrp==i,]$cvpred <- predict(fit, newdata=data.screen.abn.neg[data.screen.abn.neg$randgrp==i,], type="response")
}
data.screen.abn.neg <- mutate(data.screen.abn.neg, cvpred.ntile = ntile(cvpred, 5))
data.screen.abn.neg %>% group_by(cvpred.ntile) %>% summarise(pred.cases= sum(cvpred), obs.cases = sum(case))
c(sum(data.screen.abn.neg$cvpred), sum(data.screen.abn.neg$case)) # number obs and expected cases
poisson.test(round(sum(data.screen.abn.neg$cvpred)), sum(data.screen.abn.neg$case), alternative="two.sided")

# 10-fold cross validation to get CV error. The first delta is standard version; second is bias-corrected.
set.seed(61116)
cv.err.int <- cv.glm(data.interval.abn, glm.int.abn, K=10)
cv.err.int$delta
cv.err.screen.neg <- cv.glm(data.screen.abn.neg, glm.screen.neg.abn, K=10)
cv.err.screen.neg$delta


### Calculate AUCs ###
## Regular AUCs. By default, the 95% CI are computed with 2000 stratified bootstrap replicates.
library(pROC)
# Interval cancer model
with(filter(data.interval.abn, interval==1), roc(case, post.risk.abn, ci=T, plot=T))  # T0-T1 - change interval for T1-T2, post-T2
# Next-screen model among negatives
with(filter(data.screen.abn.neg, interval==1), roc(case, post.risk.abn, ci=T, plot=T))   # T1 - change interval for T2
## Optimism-corrected AUCs - have to use logistic models for this, and have to actually add the interaction terms to the dataset.
library(rms)
data.screen.abn.neg <- mutate(data.screen.abn.neg, logit1yriskconsolidation = logit1yrisk*consolidation, logit1yriskemphysema = logit1yrisk*emphysema)
data.interval.abn <- mutate(data.interval.abn, logit1yriskadenopconsol = logit1yrisk*adenop.consol)
# Interval cancer model
mod.int <- lrm(case ~ logit1yrisk + logit1yriskadenopconsol -1, x=T, y=T, data=data.interval.abn)
set.seed(61116)
validate(mod.int, B=1000)
c(0.5*(0.5072+1), 0.5*(0.5082+1)) # AUC = 0.5(Dxy+1). Naive, optimism-corrected AUCs - 0.75 is OC-AUC
# Next-screen model among negatives
mod.ns <- lrm(case ~ logit1yrisk + logit1yriskconsolidation + logit1yriskemphysema -1, x=T, y=T, data=data.screen.abn.neg)
set.seed(61116)
validate(mod.ns, B=1000)     # AUC = 0.5(Dxy+1)
c(0.5*(0.4760+1), 0.5*(0.4689+1)) # AUC = 0.5(Dxy+1). Naive, optimism-corrected AUCs - 0.73 is OC-AUC

#### Figures #### 

# Effect of screen findings on risk of INTERVAL ca among screen-negatives
med.risk.interv.prescr <- median(filter(data.interval.abn, interval==1)$prescr.1yrisk, na.rm=T)
med.risk.interv.post.noac <- median(filter(data.interval.abn, interval==1 & adenop.consol==0)$post.risk.abn, na.rm=T)
med.risk.interv.post.ac <- median(filter(data.interval.abn, interval==1 & adenop.consol==1)$post.risk.abn, na.rm=T)
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
  geom_boxplot(data=subset(data.interval.abn, interval==1), aes(x=0, y=prescr.1yrisk), lwd=1, width=0.4, outlier.shape=NA) +
  geom_boxplot(data=subset(data.interval.abn, interval==1 & adenop.consol==0), aes(x=0.9, y=post.risk.abn), lwd=1, width=0.98*.8, outlier.shape=NA) +
  geom_boxplot(data=subset(data.interval.abn, interval==1 & adenop.consol==1), aes(x=1.1, y=post.risk.abn), lwd=1, width=0.02*.8, outlier.shape=NA) +
  geom_segment(aes(x=0, y=med.risk.interv.prescr, xend=0.9, yend=med.risk.interv.post.noac), linetype="dashed", size=0.6) +
  geom_segment(aes(x=0, y=med.risk.interv.prescr, xend=1.1, yend=med.risk.interv.post.ac), linetype="dashed", size=0.6) +
  annotate(geom="text", x=0.6, y=0.0062, label = "Adenopathy or consolidation (2%)", angle=4, size=9) +
  annotate(geom="text", x=0.45, y=0.0028, label = "Neither noted (98%)", angle=-4, size=9) +
  annotate(geom="text", x=0, y=0.045, label="Pre-screening risk", size=9) +
  annotate(geom="text", x=1, y=0.045, label="Risk during 1-year interval", size=9)
dev.off()
# Numbers for the text
with(data.interval.abn, CrossTable(interval, adenop.consol)) # 1.8% have adenop or consol at T0
# Some percentiles for below the figure
neg.i.psr.q <- quantile(filter(data.interval.abn, interval==1)$prescr.1yrisk, probs=c(0.25, 0.5, 0.75))
neg.i.no.q <- quantile(filter(data.interval.abn, interval==1 & adenop.consol==0)$post.risk.abn, probs=c(0.25, 0.5, 0.75))
neg.i.adcon.q <- quantile(filter(data.interval.abn, interval==1 & adenop.consol==1)$post.risk.abn, probs=c(0.25, 0.5, 0.75))
rbind(neg.i.psr.q, neg.i.no.q, neg.i.adcon.q)  # print the quantiles for each group
c(neg.i.no.q[2]/neg.i.psr.q[2], neg.i.adcon.q[2]/neg.i.psr.q[2]) # median RRs for no, adenop.consol
c(neg.i.no.q[2]-neg.i.psr.q[2], neg.i.adcon.q[2]-neg.i.psr.q[2]) # median RDs for no, adenop.consol
range_without_outliers(filter(data.interval.abn, interval==1 & adenop.consol==0)$post.risk.abn) # this uses my function defined in hilary_functions.R

# Effect of screen findings on risk of SCREEN-DETECTED ca among screen negatives
    # Note: I am not making a boxplot for the N=36 with emphysema and consolidation.
med.risk.screen.neg.prescr <- median(filter(data.screen.abn.neg, interval==1)$prescr.1yrisk)
med.risk.screen.neg.neither <- median(filter(data.screen.abn.neg, interval==1 & emphysema==0 & consolidation==0)$post.risk.abn)
med.risk.screen.neg.emp <-  median(filter(data.screen.abn.neg, interval==1 & emphysema==1)$post.risk.abn)
med.risk.screen.neg.consol <-  median(filter(data.screen.abn.neg, interval==1 & consolidation==1)$post.risk.abn)
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
  geom_boxplot(data=subset(data.screen.abn.neg, interval==1), aes(x=0, y=prescr.1yrisk), lwd=1, width=0.4, outlier.shape=NA) +
  geom_boxplot(data=subset(data.screen.abn.neg, interval==1 & emphysema==0 & consolidation==0), aes(x=0.8, y=post.risk.abn), lwd=1, width=0.69*1.2, outlier.shape=NA) +
  geom_boxplot(data=subset(data.screen.abn.neg, interval==1 & emphysema==1), aes(x=1, y=post.risk.abn), lwd=1, width=.3*1.2, outlier.shape=NA) +
  geom_boxplot(data=subset(data.screen.abn.neg, interval==1 & consolidation==1), aes(x=1.2, y=post.risk.abn), lwd=1, width=.01*1.2, outlier.shape=NA) +
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
with(filter(data.screen.abn.neg, interval==1), CrossTable(emp))
with(filter(data.screen.abn.neg, interval==1), CrossTable(emphysema))
with(filter(data.screen.abn.neg, interval==1), CrossTable(consolidation))
with(filter(data.screen.abn.neg, interval==1), CrossTable(emphysema, consolidation))
with(data.screen.abn.neg, CrossTable(interval, I(emphysema==0 & consolidation==0))) # no emphysema NOR consolidation (70% at T0)
  # Thus .738 *.696 = 51% of participants would be screen-negative without emp or consol.
# Some percentiles for below the figure
quantile(filter(data.screen.abn.neg, interval==1)$post.risk.abn, probs=c(0.25, 0.5, 0.75, 0.8, 0.85, 0.9)) # overall quantiles among all negatives
neg.s.psr.q <- quantile(filter(data.screen.abn.neg, interval==1)$prescr.1yrisk, probs=c(0.25, 0.5, 0.75))
neg.s.no.no.q <- quantile(filter(data.screen.abn.neg, interval==1 & emphysema==0 & consolidation==0)$post.risk.abn, probs=c(0.25, 0.5, 0.75))
neg.s.emp.q <- quantile(filter(data.screen.abn.neg, interval==1 & emphysema==1)$post.risk.abn, probs=c(0.25, 0.5, 0.75))
neg.s.consol.q <- quantile(filter(data.screen.abn.neg, interval==1 & consolidation==1)$post.risk.abn, probs=c(0.25, 0.5, 0.75))
rbind(neg.s.psr.q, neg.s.no.no.q, neg.s.emp.q, neg.s.consol.q)   # print the quantiles for each group
rbind(neg.s.no.no.q[2]/neg.s.psr.q[2], neg.s.emp.q[2]/neg.s.psr.q[2], neg.s.consol.q[2]/neg.s.psr.q[2]) # median RRs for no-no, emp, consol
rbind(neg.s.no.no.q[2]-neg.s.psr.q[2], neg.s.emp.q[2]-neg.s.psr.q[2], neg.s.consol.q[2]-neg.s.psr.q[2]) # median RDs for no-no, emp, consol


# Potential risk THRESHOLDS for longer interval after negative screen. Use T0-negatives and cases at T1
  # Note: 64 cancers at T1; only 30 among the 70% with no emp or consol.
# Make a dataset for impact of different risk thresholds among all screen-negatives
num_T0_neg <- nrow(filter(data.screen.abn.neg, interval==1))        # change to interval==2 to look at T2
num_T0_all <- nrow(filter(data.screen.abn, interval==1))            # change to interval==2 to look at T2
cases_T1 <- sum(filter(data.screen.abn.neg, interval==1)$case)      # change to interval==2 to look at T2
n <- perc_negs <- perc_all <- ca_N <- perc_of_ca <- perc_w_ca <- perc_w_emp <- perc_w_consol <- perc_w_emp_or_consol <- vector()
thresholds <- seq(0,0.1,0.0001)
for (i in seq_along(thresholds)) {
  dat <- filter(data.screen.abn.neg, interval==1 & post.risk.abn<=thresholds[i])  # change to interval==2 to look at T2
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



