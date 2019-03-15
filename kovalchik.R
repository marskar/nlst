c.coxph.risk <- function(coxph,...){
    models <- list(coxph)
    if(!missing(..1))
        models <- c(models, list(...))
    models  
}

coxph.risk.initialize <- function(coxph,...){
    
    models <- c.coxph.risk(coxph,...)
    relrisk <- lapply(models, projection.relrisk)
    baseline <- lapply(models, survfit.risk)
    
    list(
        models = models,
        relrisk = relrisk,
        baseline = baseline
    )
}
projection.relrisk <- function(object, data){
    
    if (is.numeric(object))
        return(object)
    else if (class(object)=="coxph")
        if (missing(data))
            return(coxph.relrisk.uncentered(object)) # RETURNS RELRISKS FOR FULL DATASET
    else
        return(coxph.relrisk.uncentered(object, data))
    else
        stop(cat("No method for class",class(object)))
    
}
coxph.relrisk.uncentered <- function(coxph.object, newdata){
    
    # RETURN VECTOR OF INDIVIDUAL RELATIVE RISKS FORM COXPH.OBJECT
    # FASTER IMPLEMENTATION THAN THROUGH FULL DESIGN
    center <- coxph.object$means%*%coef(coxph.object)
    
    if(!missing(newdata))
        lp <- predict(coxph.object, newdata, type="lp")
    else
        lp <- coxph.object$linear.predictor
    
    exp(lp+center)	
}
coxph.risk.baseline <- function(begin,end,models){
    
    in.interval <- function(x, begin, end) 
        x >= begin & x <= end
    
    # TAKE THE SET OF EVENT TIMES AND RETURN THE CUMULATIVE HAZARDS FOR EACH MODEL
    # IF NO TIMES, THE HAZARD IS ZERO
    event.times <- survfit(models[[1]])
    event.times <- event.times$time[event.times$n.event>=1]
    which.event.times <- in.interval(event.times, begin, end)
    
    if(all(!which.event.times))
        NA
    else{
        
        event.times <- event.times[which.event.times]
        
        Hazards <- mapply(basehaz.coxph.risk.revised, object=models,
                          MoreArgs=list(times=event.times),SIMPLIFY=FALSE)
        Hazards
    }
}

basehaz.coxph.risk.revised <- function(object, times){
    
    sfit <- summary(survfit(object), time=times, extend=TRUE)
    H <- -log(sfit$surv)
    z0 <- object$means
    bz0 <- sum(z0 * coef(object))
    H <- H * exp(-bz0)
    S <- exp(-H) # BASELINE
    
    # WHEN START OF PROJECTION IS NOT BEFORE FIRST EVENT
    # 	ADJUST FIRST HAZARD JUMP; ONLY IMPACTS PRIMARY EVENT
    
    H0 <- -log(summary(survfit(object), time=times[1]*.95, extend=TRUE)$surv)
    h <- c(H[1]-H0,diff(H))
    
    data.frame(time=times,haz=h,surv=S/S[1])
}
risk.kovalchik <-
function (begin, end, newdata, coxph1, ...) 
{
    #calculate absolute risk given hazards/survival and relative risks
    risk.fixed.interval <- function(H, RR) {
        if (!is.list(H)) 
            0
        else {
            absrisk <- H[[1]]$surv^RR[1] * H[[1]]$haz * RR[1]
            for (i in 2:length(H)) {
                absrisk <- absrisk * (H[[i]]$surv^RR[i])
            }
            sum(absrisk)
        }
    }
    
    models <- c.coxph.risk(coxph1, ...)

    #subset data to all covariates given in the models
    AllVars <- unique(unlist(sapply(models, function(x) all.vars(x$formula))))
    if ("years.followed" %in% AllVars) { AllVars <- AllVars[-grep("years.followed",AllVars)] }
    if ("lung.cancer.death" %in% AllVars) { AllVars <- AllVars[-grep("lung.cancer.death",AllVars)] }
    if ("other.cause.death" %in% AllVars) { AllVars <- AllVars[-grep("other.cause.death",AllVars)] }   # HAR changed other.cancer.death to other.cause.death on 17 Nov 2017
    if ("incidence.years" %in% AllVars) { AllVars <- AllVars[-grep("incidence.years",AllVars)] }
    if ("case" %in% AllVars) { AllVars <- AllVars[-grep("case",AllVars)] }
    newdata <- subset(newdata, select = unique(AllVars))
    
    #check for missing
    which.kept <- complete.cases(newdata)
    if (!all(which.kept)) {
        warning("Missing cases excluded.")
        if (length(begin) > 1) {
            begin <- begin[which.kept]
            end <- end[which.kept]
        }
        newdata <- newdata[which.kept, ]
    }

    #calculate relative risk for each subject and store as list
    rr <- sapply(models, projection.relrisk, data = newdata)

    if (is.matrix(rr))
        rr.list <- lapply(1:nrow(rr), function(x) rr[x, ])
    else rr.list <- list(rr)

    #estimate risk
    if (length(begin) == 1) {
        H <- coxph.risk.baseline(begin, end, models)   #hazard for each timepoint for lung cancer death and other deaths
        risks <- mapply(risk.fixed.interval, RR = rr.list, MoreArgs = list(H = H))
    }
    else {
        risks <- mapply(risk, begin = begin, end = end, RR = rr.list, 
            MoreArgs = list(models = models))
    }
    risks
}
