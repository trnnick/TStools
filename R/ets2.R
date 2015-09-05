ets2 <- function(data, model="ZZZ", persistence=NULL, phi=NULL,
                 bounds=c("usual","admissible"),
                 initial=NULL, initial.season=NULL,
                 IC=c("AICc","AIC"), trace=FALSE, CF.type=c("TLV","GV","TV","hsteps"),
                 intervals=FALSE, int.w=0.95, xreg=NULL,
                 holdout=FALSE, h=10, silent=FALSE, legend=TRUE,
                 ...){

# Start measuring the time of calculations
    start.time <- Sys.time()
    
    bounds <- substring(bounds[1],1,1)
    IC <- IC[1]
    CF.type <- CF.type[1]

# Check if the data is vector
    if(!is.numeric(data) & !is.ts(data)){
        stop("The provided data is not a vector or ts object! Can't build any model!", call.=FALSE)
    }

# Check if CF.type is appropriate in the case of trace==TRUE
    if(trace==TRUE){
        if(CF.type!="GV" & CF.type!="TLV" & CF.type!="TV" & CF.type!="hsteps"){
            message(paste0("The strange Cost Function is chosen: ",CF.type))
            sowhat(CF.type)
            message("Switching to 'TLV'")
            CF.type <- "TLV"
        }
    }

# Check if "bounds" parameter makes any sense
    if(bounds!="u" & bounds!="a"){
        message("The strange bounds are defined. Switching to 'usual'.")
        bounds <- "u"
    }

# If chosen model is "AAdN" or anything like that, we are taking the appropriate values
    if(nchar(model)==4){
        Etype <- substring(model,1,1)
        Ttype <- substring(model,2,2)
        Stype <- substring(model,4,4)
        damped <- TRUE
        if(substring(model,3,3)!="d"){
            message(paste0("You have defined a strange model: ",model))
            sowhat(model)
            model <- paste0(Etype,Ttype,"d",Stype)
        }
    }
    else if(nchar(model)==3){
        Etype <- substring(model,1,1)
        Ttype <- substring(model,2,2)
        Stype <- substring(model,3,3)
        damped <- FALSE
    }
    else{
        message(paste0("You have defined a strange model: ",model))
        sowhat(model)
        message("Switching to 'ZZZ'")
        model <- "ZZZ"

        Etype <- "Z"
        Ttype <- "Z"
        Stype <- "Z"
        damped <- TRUE
    }

    if(any(is.na(data))){
        message("Data contains NAs. These observations will be excluded.")
        datanew <- data[!is.na(data)]
        if(is.ts(data)){
            datanew <- ts(datanew,start=start(data),frequency=frequency(data))
        }
        data <- datanew
    }

# Define obs.all, the overal number of observations (in-sample + holdout)
    obs.all <- length(data) + (1 - holdout)*h

# Define obs, the number of observations of in-sample
    obs <- length(data) - holdout*h

# Define the actual values
    y <- coredata(data)

# Check if the data is ts-object
    if(!is.ts(data) & Stype!="N"){
        message("The provided data is not ts object. Only non-seasonal models are available.")
        Stype <- "N"
    }
    seasfreq <- frequency(data)

# Check the length of the provided data
    if(Stype!="N" & (obs/seasfreq)<2){
        message("Not enough observations for the seasonal model. Only non-seasonal models are available.")
        Stype <- "N"
    }

# If model selection is chosen, forget about the initial values and persistence
    if(Etype=="Z" | Ttype=="Z" | Stype=="Z"){
        if(!is.null(initial) | !is.null(initial.season) | !is.null(persistence) | !is.null(phi)){
            message("Model selection doesn't go well with the predefined values.")
            message("Switching to the estimation of all the parameters.")
            initial <- NULL
            initial.season <- NULL
            persistence <- NULL
            phi <- NULL
        }
    }

### Check all the parameters for the possible errors.
    if(!is.null(persistence) & length(persistence)>3){
        message("The length of persistence vector is wrong! It should not be greater than 3.")
        message("Changing to the estimation of persistence vector values.")
        persistence <- NULL
    }

### Check the error type
    if(Etype!="Z" & Etype!="A" & Etype!="M"){
        message("Wrong error type! Should be 'Z', 'A' or 'M'.")
        message("Changing to 'Z'")
        Etype <- "Z"
    }

### Check the trend type
    if(Ttype!="Z" & Ttype!="N" & Ttype!="A" & Ttype!="M"){
        message("Wrong trend type! Should be 'Z', 'N', 'A' or 'M'.")
        message("Changing to 'Z'")
        Ttype <- "Z"
    }

### Check the seasonaity type
    if(Stype!="Z" & Stype!="N" & Stype!="A" & Stype!="M"){
        message("Wrong seasonality type! Should be 'Z', 'N', 'A' or 'M'.")
        if(seasfreq==1){
            if(silent==FALSE){
                message("Data is non-seasonal. Changing seasonal component to 'N'")
            }
            Stype <- "N"
        }
        else{
            message("Changing to 'Z'")
            Stype <- "Z"
        }
    }
    if(Stype!="N" & seasfreq==1){
        if(silent==FALSE){
            message("Cannot build the seasonal model on the data with the frequency 1.")
            message(paste0("Switching to non-seasonal model: ETS(",substring(model,1,nchar(model)-1),"N)"))
        }
        Stype <- "N"
    }


    if(any(y<=0)){
        if(Etype=="M"){
            message("Can't apply multiplicative model to non-positive data. Switching error to 'A'")
            Etype <- "A"
        }
        if(Ttype=="M"){
            message("Can't apply multiplicative model to non-positive data. Switching trend to 'A'")
            Ttype <- "A"
        }
        if(Stype=="M"){
            message("Can't apply multiplicative model to non-positive data. Switching seasonal to 'A'")
            Stype <- "A"
        }
    }


# Function fills in the initial values of xt estimating the series.
define.xt <- function(Ttype,Stype,seasfreq,n.components){
    xt <- matrix(NA,seasfreq,n.components)

    if(Ttype=="N"){
        xt[1:seasfreq,1] <- rep(mean(y[1:min(12,obs)]),times=seasfreq)
    }
    else if(Ttype=="A"){
        slope <- cov(y[1:min(12,obs)],c(1:min(12,obs)))/var(c(1:min(12,obs)))
        intercept <- mean(y[1:min(12,obs)]) - slope * (mean(c(1:min(12,obs))) - 1)
        xt[1:seasfreq,1] <- rep(intercept,seasfreq)
        xt[1:seasfreq,2] <- rep(slope,seasfreq)
    }
    else if(Ttype=="M"){
        xt[1:seasfreq,1] <- rep(mean(y[1:min(12,obs)]),seasfreq)
        xt[1:seasfreq,2] <- rep(1,seasfreq)
    }

    if(Stype=="A"){
        xt[1:seasfreq,n.components] <- decompose(ts(y[1:obs],frequency=seasfreq),type="additive")$seasonal[1:seasfreq]
    }
    else if(Stype=="M"){
        xt[1:seasfreq,n.components] <- decompose(ts(y[1:obs],frequency=seasfreq),type="multiplicative")$seasonal[1:seasfreq]
    }

    return(xt)
}

# Function defines the number of components, lags and seasfreq depending on the chosen model
define.param <- function(Ttype,Stype,damped,phi){
# The number of components of ETS
    n.components <- 1
# The lag of components (needed mainly for the seasonal models)
    lags <- 1
# The names of the state vector components
    component.names <- "level"

    if(Ttype!="N"){
        n.components <- n.components + 1
        lags <- c(lags,1)
        component.names <- c(component.names,"trend")
        trend.component <- TRUE
    }
    else{
        trend.component <- FALSE
    }

    if(Stype!="N"){
        n.components <- n.components + 1
        lags <- c(lags,frequency(data))
        component.names <- c(component.names,"seasonality")
        seasonal.component <- TRUE
        seasfreq <- frequency(data)
    }
    else{
        seasonal.component <- FALSE
# In the cases of non-seasonal models built on seasonal data seasonality is equal to 1
        seasfreq <- 1
    }

### Create the essential matrices and vectors
# Create the matrix of state vectors.
# First seasfreq rows are initial values, excluded from estimation.
    matxt <- matrix(NA,nrow=(obs+seasfreq),ncol=n.components)
    colnames(matxt) <- component.names

# If the persistence vector is provided, use it
    if(!is.null(persistence)){
        vecg <- persistence
        estimate.persistence <- FALSE
    }
    else{
        if(Ttype=="M" | Stype=="M"){
            vecg <- rep(0.05,n.components)
            estimate.persistence <- TRUE
        }
        else{
            vecg <- c(0.3,0.2,0.1)[1:n.components]
            estimate.persistence <- TRUE
        }
    }

# If phi is not provided, mark that
    if(damped==TRUE){
        if(is.null(phi)){
            estimate.phi <- TRUE
            phi <- 0.95
        }
        else{
            if(phi<0 | phi>2){
                message("Damping parameter should lie in (0, 2) region.")
                message("Changing to the estimation of phi.")
                phi <- 0.95
                estimate.phi <- TRUE
            }
            else{
                estimate.phi <- FALSE
            }
        }
    }
    else{
        phi <- 1
        estimate.phi <- FALSE
    }

# If initial values are provided, write them down
    if(!is.null(initial)){
        matxt[1:seasfreq,1:(n.components - seasonal.component)] <- rep(initial,each=seasfreq)
        estimate.initial <- FALSE
    }
    else{
        matxt[1:seasfreq,1:(n.components - seasonal.component)] <- define.xt(Ttype,Stype,seasfreq,n.components)[1:seasfreq,1:(n.components - seasonal.component)]
        estimate.initial <- TRUE
    }

# If the seasonal model is chosen and initials are provided, fill in the first "seasfreq" values of seasonal component.
    if(seasonal.component==TRUE & !is.null(initial.season)){
        matxt[1:seasfreq,n.components] <- initial.season
        estimate.initial.season <- FALSE
    }
    else if(seasonal.component==TRUE & is.null(initial.season)){
        matxt[1:seasfreq,n.components] <- define.xt(Ttype,Stype,seasfreq,n.components)[1:seasfreq,n.components]
        estimate.initial.season <- TRUE
    }
    else{
        estimate.initial.season <- FALSE
    }

    return(list(n.components=n.components,lags=lags,seasfreq=seasfreq,
                matxt=matxt,vecg=vecg,phi=phi,
                trend.component=trend.component,
                seasonal.component=seasonal.component,
                estimate.persistence=estimate.persistence,
                estimate.phi=estimate.phi,
                estimate.initial=estimate.initial,
                estimate.initial.season=estimate.initial.season))
}

# Function returns transition matrix F and measurement matrix w depending on chosen model
mat.ets <- function(trend.component,seasonal.component,phi){
    if(seasonal.component==FALSE){
        if(trend.component==FALSE){
            matF <- matrix(1,1,1)
            matw <- 1
        }
        else{
            matF <- matrix(c(1,0,phi,phi),2,2)
            matw <- c(1,phi)
        }
    }
    else{
        if(trend.component==FALSE){
            matF <- matrix(c(1,0,0,1),2,2)
            matw <- c(1,1)
        }
        else{
            matF <- matrix(c(1,0,0,phi,phi,0,0,0,1),3,3)
            matw <- c(1,phi,1)
        }
    }
    return(list(matF=matF,matw=matw))
}

# Function fills in the initial values of matxt using values of C
estim.values <- function(matxt,vecg,phi,C,n.components,seasfreq,seasonal.component,Stype){

# If the persistence vector is provided, use it
    if(estimate.persistence==TRUE){
        vecg <- C[1:n.components]
    }

# If phi was not provided, use it in estimate
    if(estimate.phi==TRUE){
        phi <- C[n.components*estimate.persistence + 1]
    }

# If initial values are not provided, write them in the matrix
    if(estimate.initial==TRUE){
        matxt[1:seasfreq,1:(n.components - seasonal.component)] <- rep(C[(n.components*estimate.persistence + estimate.phi + 1):(n.components*estimate.persistence + estimate.phi + n.components - seasonal.component)],each=seasfreq)
    }

# If the seasonal model is chosen and initials are provided, fill in the first "seasfreq" values of seasonal component.
    if(seasonal.component==TRUE){
        if(estimate.initial.season==TRUE){
            matxt[1:seasfreq,n.components] <- C[(n.components*estimate.persistence + estimate.phi + (n.components - seasonal.component)*estimate.initial + 1):(n.components*estimate.persistence + estimate.phi + (n.components - seasonal.component)*estimate.initial + seasfreq)]
            if(Stype=="A"){
                matxt[1:seasfreq,n.components] <- matxt[1:seasfreq,n.components]-mean(matxt[1:seasfreq,n.components])
            }
            else if(Stype=="M"){
                matxt[1:seasfreq,n.components] <- exp(log(matxt[1:seasfreq,n.components])-mean(log(matxt[1:seasfreq,n.components])))
            }
        }
    }

    return(list(vecg=vecg,phi=phi,matxt=matxt))
}

# Cost function for ETS
CF <- function(C){

    init.ets <- estim.values(matxt,vecg,phi,C,n.components,seasfreq,seasonal.component,Stype)
    vecg <- init.ets$vecg
    phi <- init.ets$phi
    matxt <- init.ets$matxt

    matrices <- mat.ets(trend.component,seasonal.component,phi)
    matF <- matrices$matF
    matw <- matrices$matw

    if(estimate.persistence==TRUE){
        if(bounds=="a" & trend.component==TRUE & seasonal.component==TRUE){
            Theta.func <- function(Theta){
                return(abs((phi*C[1]+phi+1)/(C[3]) + ((phi-1)*(1+cos(Theta)-cos(seasfreq*Theta))+cos((seasfreq-1)*Theta)-phi*cos((seasfreq+1)*Theta))/(2*(1+cos(Theta))*(1-cos(seasfreq*Theta)))))
            }
            Theta <- 0.1
            Theta <- suppressWarnings(optim(Theta,Theta.func,method="Brent",lower=0,upper=1)$par)
        }
        else{
            Theta <- 0
        }
        CF.res <- costfunc(matxt,matF,matrix(matw,1,length(matw)),as.matrix(y[1:obs]),matrix(vecg,length(vecg),1),h,Etype,Ttype,Stype,seasfreq,trace,CF.type,normalizer,bounds,phi,Theta)
    }
    else{
        CF.res <- optimizerwrap(matxt,matF,matrix(matw,1,length(matw)),as.matrix(y[1:obs]),matrix(vecg,length(vecg),1),h,Etype,Ttype,Stype,seasfreq,trace,CF.type,normalizer)        
    }

    if(is.nan(CF.res) | is.na(CF.res) | is.infinite(CF.res)){
        CF.res <- 1e100
    }

    return(CF.res)
}

# Function returns interval forecasts.
int.ets <- function(xt,matF,matw,vecg,h=1){
    y.lo <- rep(NA,h)
    y.up <- y.lo
    return(list(y.up=y.up,y.lo=y.lo))
}

MASE.lvl <- function(a,f,round=3){
    # This function calculates Mean Absolute Scaled Error using level of the actual series
    # a - actual values,
    # f - forecasted values.
    result <- round(sum(abs(a-f),na.rm=TRUE)/sum(a),round)
    return(result)
}

MASE <- function(a,f,scale,round=3){
    # This function calculates Mean Absolute Scaled Error as in Hyndman & Koehler, 2006
    # a - actual values,
    # f - forecasted values.
    # scale - the measure to scale errors with. Usually - MAE of in-sample.
    result <- round(mean(abs(a-f),na.rm=TRUE)/scale,round)
    return(result)
}

# Function constructs default bounds where C values should lie
C.values <- function(bounds,Ttype,Stype,vecg,matxt,phi,seasfreq,n.components,seasonal.component){
    if(bounds=="u"){
        C <- NA
        C.lower <- NA
        C.upper <- NA

        if(estimate.persistence==TRUE){
            C <- c(C,vecg)
            C.lower <- c(C.lower,rep(0,length(vecg)))
            C.upper <- c(C.upper,rep(1,length(vecg)))
        }
        if(estimate.phi==TRUE){
            C <- c(C,phi)
            C.lower <- c(C.lower,0)
            C.upper <- c(C.upper,1)
        }
        if(estimate.initial==TRUE){
            C <- c(C,matxt[seasfreq,1:(n.components - seasonal.component)])
            if(Ttype!="M"){
                C.lower <- c(C.lower,rep(-Inf,(n.components - seasonal.component)))
                C.upper <- c(C.upper,rep(Inf,(n.components - seasonal.component)))
            }
            else{
                C.lower <- c(C.lower,1,0.01)
                C.upper <- c(C.upper,Inf,3)
            }
        }
        if(Stype!="N"){
            if(estimate.initial.season==TRUE){
                C <- c(C,matxt[1:seasfreq,n.components])
                if(Stype=="A"){
                    C.lower <- c(C.lower,rep(-Inf,seasfreq))
                    C.upper <- c(C.upper,rep(Inf,seasfreq))
                }
                else{
                    C.lower <- c(C.lower,rep(0,seasfreq))
                    C.upper <- c(C.upper,rep(10,seasfreq))
                }
            }
        }
    }
    else{
        C <- NA
        C.lower <- NA
        C.upper <- NA

        if(estimate.persistence==TRUE){
            C <- c(C,vecg)
            C.lower <- c(C.lower,rep(-5,length(vecg)))
            C.upper <- c(C.upper,rep(5,length(vecg)))
        }
        if(estimate.phi==TRUE){
            C <- c(C,phi)
            C.lower <- c(C.lower,0)
            C.upper <- c(C.upper,1)
        }
        if(estimate.initial==TRUE){
            C <- c(C,matxt[seasfreq,1:(n.components - seasonal.component)])
            if(Ttype!="M"){
                C.lower <- c(C.lower,rep(-Inf,(n.components - seasonal.component)))
                C.upper <- c(C.upper,rep(Inf,(n.components - seasonal.component)))
            }
            else{
                C.lower <- c(C.lower,1,0.01)
                C.upper <- c(C.upper,Inf,3)
            }
        }
        if(Stype!="N"){
            if(estimate.initial.season==TRUE){
                C <- c(C,matxt[1:seasfreq,n.components])
                if(Stype=="A"){
                    C.lower <- c(C.lower,rep(-Inf,seasfreq))
                    C.upper <- c(C.upper,rep(Inf,seasfreq))
                }
                else{
                    C.lower <- c(C.lower,rep(-0.0001,seasfreq))
                    C.upper <- c(C.upper,rep(20,seasfreq))
                }
            }
        }
    }

    C <- C[!is.na(C)]
    C.lower <- C.lower[!is.na(C.lower)]
    C.upper <- C.upper[!is.na(C.upper)]
    
    return(list(C=C,C.lower=C.lower,C.upper=C.upper))
}

Likelihood.value <- function(C){
    if(trace==TRUE & (CF.type=="TLV" | CF.type=="GV")){
        return(-obs/2 *((h^trace)*log(2*pi*exp(1)) + CF(C)))
    }
    else{
        return(-obs/2 *((h^trace)*log(2*pi*exp(1)) + log(CF(C))))
    }
}

## Function calculates ICs
IC.calc <- function(CF.objective=CF.objective,n.param=n.param,C,Etype=Etype){
# Information criteria are calculated with the constant part "log(2*pi*exp(1)*h)*obs".
# And it is based on the mean of the sum squared residuals either than sum.
# Hyndman likelihood is: llikelihood <- obs*log(obs*CF.objective)
    llikelihood <- Likelihood.value(C)


    AIC.coef <- 2*n.param*h^trace - 2*llikelihood
    AICc.coef <- AIC.coef + 2 * n.param * (n.param + 1) / (obs - n.param - 1)
    BIC.coef <- log(obs)*n.param*h^trace - 2*llikelihood

    ICs <- c(AIC.coef, AICc.coef, BIC.coef)
    names(ICs) <- c("AIC", "AICc", "BIC")

    return(list(llikelihood=llikelihood,ICs=ICs))
}

ets2.auto <- function(Etype,Ttype,Stype,IC="AICc",CF.type="none"){
# Script for the automatic model selection based on chosen IC.

# Define the pool of models to select from
    if(any(y<=0)){
        if(silent==FALSE){
            message("Only additive models are allowed with the negative data.")
        }
        errors.pool <- c("A")
        trends.pool <- c("N","A","Ad")
        season.pool <- c("N","A")
    }
    else{
        errors.pool <- c("A","M")
        trends.pool <- c("N","A","Ad","M","Md")
        season.pool <- c("N","A","M")
    }

    if(Etype!="Z"){
        errors.pool <- Etype
    }

    if(Ttype!="Z"){
        if(damped==TRUE){
            trends.pool <- paste0(Ttype,"d")
        }
        else{
            trends.pool <- Ttype
        }
    }

    if(Stype!="Z"){
        season.pool <- Stype
    }

# Number of observations in the mat.error matrix excluding NAs.
    errors.mat.obs <- obs - h + 1

    models.number <- (length(trends.pool)*length(season.pool)*length(errors.pool))
    models.pool <- array(c(1:models.number),
                         c(length(trends.pool),length(season.pool),length(errors.pool)),
                         dimnames=list(trends.pool,season.pool,errors.pool))

    results <- as.list(c(1:models.number))

    if(silent==FALSE){
        cat("Building model: ")
    }    
# Start cycle of models
    for(j in 1:models.number){

        Ttype <- dimnames(models.pool)[[1]][which(models.pool==j,arr.ind=TRUE)[1]]
        if(nchar(Ttype)==2){
            Ttype <- substring(Ttype,1,1)
            damped <- TRUE
            phi <- NULL
        }
        else{
            damped <- FALSE
            phi <- 1
        }
        
        Stype <- dimnames(models.pool)[[2]][which(models.pool==j,arr.ind=TRUE)[2]]
        Etype <- dimnames(models.pool)[[3]][which(models.pool==j,arr.ind=TRUE)[3]]

        Etype <<- Etype
        Ttype <<- Ttype
        Stype <<- Stype
        damped <<- damped
        phi <<- phi

        if(damped==TRUE){
            current.model <- paste0(Etype,Ttype,"d",Stype)
        }
        else{
            current.model <- paste0(Etype,Ttype,Stype)
        }
        if(silent==FALSE){
            cat(paste0(current.model," "))
        }        
        
        param.values <- define.param(Ttype,Stype,damped,phi)
        n.components <<- param.values$n.components
        lags <<- param.values$lags
        seasfreq <<- param.values$seasfreq
        matxt <<- param.values$matxt
        vecg <<- param.values$vecg
        phi <- param.values$phi
        trend.component <<- param.values$trend.component
        seasonal.component <<- param.values$seasonal.component
        estimate.persistence <<- param.values$estimate.persistence
        estimate.phi <<- param.values$estimate.phi
        estimate.initial <<- param.values$estimate.initial
        estimate.initial.season <<- param.values$estimate.initial.season

        Cs <- C.values(bounds,Ttype,Stype,vecg,matxt,phi,seasfreq,n.components,seasonal.component)
        C <- Cs$C
        C.upper <- Cs$C.upper
        C.lower <- Cs$C.lower
#        res <- nloptr::cobyla(C, CF, hin=hin.constrains, lower=C.lower, upper=C.upper)
#        CF.objective <- res$value
        res <- nloptr::nloptr(C, CF, lb=C.lower, ub=C.upper,
                              opts=list("algorithm"="NLOPT_LN_BOBYQA", "xtol_rel"=1e-8, "maxeval"=1000))
        CF.objective <- res$objective

        init.ets <- estim.values(matxt,vecg,phi,res$solution,n.components,seasfreq,seasonal.component,Stype)
        vecg <<- init.ets$vecg
        phi <- init.ets$phi
        matxt <<- init.ets$matxt

        matrices <- mat.ets(trend.component,seasonal.component,phi)
        matF <- matrices$matF
        matw <- matrices$matw

        n.param <- n.components*estimate.persistence + estimate.phi + (n.components - seasonal.component)*estimate.initial + seasfreq*estimate.initial.season

        IC.values <- IC.calc(CF.objective=CF.objective,n.param=n.param,C=res$solution,Etype=Etype)
        ICs <- IC.values$ICs

        results[[j]] <- c(ICs,Etype,Ttype,Stype,damped,CF.objective,res$solution)
    }
    if(silent==FALSE){
        cat("... Done! \n")
    }
    IC.selection <- rep(NA,length(models.pool))
    for(i in 1:length(models.pool)){
        IC.selection[i] <- as.numeric(eval(parse(text=paste0("results[[",i,"]]['",IC,"']"))))
    }

    i <- which(IC.selection==min(IC.selection))[1]

    return(results[[i]])
}

#########################################

    param.values <- define.param(Ttype=Ttype,Stype=Stype,damped=damped,phi=phi)
    n.components <- param.values$n.components
    lags <- param.values$lags
    seasfreq <- param.values$seasfreq
    matxt <- param.values$matxt
    vecg <- param.values$vecg
    phi <- param.values$phi
    trend.component <- param.values$trend.component
    seasonal.component <- param.values$seasonal.component
    estimate.persistence <- param.values$estimate.persistence
    estimate.phi <- param.values$estimate.phi
    estimate.initial <- param.values$estimate.initial
    estimate.initial.season <- param.values$estimate.initial.season

    if(trace==TRUE & CF.type=="GV"){
        normalizer <- mean(abs(diff(y[1:obs])))
    }
    else{
        normalizer <- 0;
    }

#########################################

### Check the length of initials and persistence vectors
# Check the persistence vector length
    if(!is.null(persistence)){
        if(n.components != length(persistence)){
            message("The length of persistence vector does not correspond to the chosen model!")
            message("Values will be estimated")
            persistence <- NULL
        }
    }

# Check the inital vector length
    if(!is.null(initial)){
        if(length(initial)>2){
            message("The length of the initial value is wrong! It should not be greater than 2.")
            message("Values of initial vector will be estimated.")
            initial <- NULL
        }
        if((n.components - seasonal.component)!=length(initial)){
            message("The length of initial state vector does not correspond to the chosen model!")
            message("Values of initial vector will be estimated.")
            initial <- NULL
        }
    }

# Check the seasonal inital vector length
    if(!is.null(initial.season)){
        if(frequency(data)!=length(initial.season)){
            message("The length of seasonal initial states does not correspond to the frequency of the data!")
            message("Values of initial seasonals will be estimated.")
            initial.season <- NULL
        }
    }

# Vectors of fitted data and errors
    y.fit <- rep(NA,obs)
    errors <- rep(NA,obs)

# If we use trace, define matrix of errors.
    if(trace==TRUE){
        mat.error <- matrix(NA,nrow=obs,ncol=h)
    }
    else{
        mat.error <- matrix(NA,nrow=obs,ncol=1)
    }

# Fill in the vector of initial values and vector of constrains used in estimation
    if(estimate.persistence==TRUE | estimate.phi==TRUE | estimate.initial==TRUE | estimate.initial.season==TRUE){

# Number of observations in the mat.error matrix excluding NAs.
        errors.mat.obs <- obs - h + 1

        if(Etype=="Z" | Ttype=="Z" | Stype=="Z"){
############ ETS2.auto should be rewritten in C++!!! ############
            results <- ets2.auto(Etype,Ttype,Stype,IC=IC,CF.type=CF.type)

            Etype <- results[4]
            Ttype <- results[5]
            Stype <- results[6]
            damped <- as.logical(results[7])
            CF.objective <- as.numeric(results[8])
            C <- as.numeric(results[-c(1:8)])

            param.values <- define.param(Ttype=Ttype,Stype=Stype,damped=damped,phi=phi)
            n.components <- param.values$n.components
            lags <- param.values$lags
            seasfreq <- param.values$seasfreq
            matxt <- param.values$matxt
            vecg <- param.values$vecg
            phi <- param.values$phi
            trend.component <- param.values$trend.component
            seasonal.component <- param.values$seasonal.component
            estimate.persistence <- param.values$estimate.persistence
            estimate.phi <- param.values$estimate.phi
            estimate.initial <- param.values$estimate.initial
            estimate.initial.season <- param.values$estimate.initial.season

            init.ets <- estim.values(matxt,vecg,phi,C,n.components,seasfreq,seasonal.component,Stype)
            vecg <- init.ets$vecg
            phi <- init.ets$phi
            matxt <- init.ets$matxt

        }
        else{
            Cs <- C.values(bounds,Ttype,Stype,vecg,matxt,phi,seasfreq,n.components,seasonal.component)
            C <- Cs$C
            C.upper <- Cs$C.upper
            C.lower <- Cs$C.lower
#            res <- nloptr::cobyla(C, CF, hin=hin.constrains, lower=C.lower, upper=C.upper)
#            CF.objective <- res$value
#            C <- res$par
#            eval_g_ineq=hin.constrains,
#   
############ Can we introduce the constraints in the cost function (returning Inf if constrains are violated)? ############
            res <- nloptr::nloptr(C, CF, lb=C.lower, ub=C.upper,
                                  opts=list("algorithm"="NLOPT_LN_BOBYQA", "xtol_rel"=1e-8, "maxeval"=1000))
            CF.objective <- res$objective
            C <- res$solution

            init.ets <- estim.values(matxt,vecg,phi,C,n.components,seasfreq,seasonal.component,Stype)
            vecg <- init.ets$vecg
            phi <- init.ets$phi
            matxt <- init.ets$matxt
        }
    }

    if(damped==TRUE){
        model <- paste0(Etype,Ttype,"d",Stype)
    }
    else{
        model <- paste0(Etype,Ttype,Stype)
    }

    matrices <- mat.ets(trend.component,seasonal.component,phi)
    matF <- matrices$matF
    matw <- matrices$matw

    fitting <- fitterwrap(matxt,matF,matrix(matw,1,length(matw)),as.matrix(y[1:obs]),matrix(vecg,length(vecg),1),Etype,Ttype,Stype,seasfreq)
    matxt <- ts(fitting$matxt,start=(time(data)[1] - deltat(data)*seasfreq),frequency=frequency(data))
    y.fit <- ts(fitting$yfit,start=start(data),frequency=frequency(data))

    errors.mat <- ts(errorerwrap(matxt,matF,matrix(matw,1,length(matw)),as.matrix(y[1:obs]),h,Etype,Ttype,Stype,seasfreq,TRUE),start=start(data),frequency=frequency(data))
    colnames(errors.mat) <- paste0("Error",c(1:h))
    errors <- ts(errors.mat[,1],start=start(data),frequency=frequency(data))

    y.for <- ts(forecasterwrap(matrix(matxt[(obs+1):(obs+seasfreq),],nrow=seasfreq),matF,matrix(matw,nrow=1),h,Ttype,Stype,seasfreq),start=time(data)[obs]+deltat(data),frequency=frequency(data))

    y <- data

    if(estimate.persistence==FALSE & estimate.phi==FALSE & estimate.initial==FALSE & estimate.initial.season==FALSE){
        C <- c(vecg,phi,initial,initial.season)
        errors.mat.obs <- obs - h + 1
        CF.objective <- CF(C)
        n.param <- 0
    }
    else{
        n.param <- n.components*estimate.persistence + estimate.phi + (n.components - seasonal.component)*estimate.initial + seasfreq*estimate.initial.season
    }

    FI <- numDeriv::hessian(Likelihood.value,C)

# Calculate IC values
    IC.values <- IC.calc(CF.objective=CF.objective,n.param=n.param,C=C,Etype=Etype)
    llikelihood <- IC.values$llikelihood
    ICs <- IC.values$ICs

# Convert bounds to ts
    if(intervals==T){
        y.low <- ts(y.low,start=start(y.for),frequency=frequency(data))
        y.high <- ts(y.high,start=start(y.for),frequency=frequency(data))
    }

if(silent==FALSE){
# Define plot.range for plot
    if(intervals==T){
        plot.range <- range(min(data,y.fit,y.for,y.low),max(data,y.fit,y.for,y.high))
    }
    else{
        plot.range <- range(min(data,y.fit,y.for),max(data,y.fit,y.for))
    }
    
# Print time elapsed on the construction
    print(paste0("Time elapsed: ",round(as.numeric(Sys.time() - start.time,units="secs"),2)," seconds"))
    print(paste0("Model constructed: ",model))
    print(paste0("Persistence vector: ", paste(round(vecg,3),collapse=", ")))
    if(damped==TRUE){
        print(paste0("Damping parameter: ", round(phi,3)))
    }
    print(paste0("Initial components: ", paste(round(matxt[seasfreq,1:(n.components - seasonal.component)],3),collapse=", ")))
    if(seasonal.component==TRUE){
        print(paste0("Initial seasonal components: ", paste(round(matxt[1:seasfreq,n.components],3),collapse=", ")))
    }
    print(paste0("Residuals sigma: ",round(sqrt(mean(errors^2)),3)))
    if(trace==TRUE){
        print(paste0("CF type: trace with ",CF.type))
    }
    else{
        print(paste0("CF type: one step ahead"))
    }
    print(paste0("CF value is: ",round(CF.objective,0)))
    print(paste0("Biased log-likelihood: ",round((llikelihood - n.param*h^trace),0)))
    print(paste0("AIC: ",round(ICs["AIC"],3)," AICc: ", round(ICs["AICc"],3)))
    if(holdout==T){
        print(paste0("MASE: ",MASE(coredata(data)[(obs+1):obs.all],coredata(y.for),mean(abs(diff(coredata(data)[1:obs]))),round=3)))
        print(paste0("MASE.lvl: ",MASE.lvl(coredata(data)[(obs+1):obs.all],coredata(y.for),round=5)*100,"%"))
    }

    par(mfrow=c(1,1), mar=c(5,3,2,1))
    plot(data,type="l",xlim=range(time(data)[1],time(y.for)[h]),
         ylim=plot.range,xlab="Time", ylab="")
    lines(y.fit,col="purple",lwd=2,lty=2)

    if(intervals==T){
        if(h>1){
            lines(y.low,col="darkgrey",lwd=3,lty=2)
            lines(y.high,col="darkgrey",lwd=3,lty=2)
# Draw the nice areas between the borders
            polygon(c(seq(deltat(y.high)*(start(y.high)[2]-1)+start(y.high)[1],deltat(y.high)*(end(y.high)[2]-1)+end(y.high)[1],deltat(y.high)),
                  rev(seq(deltat(y.low)*(start(y.low)[2]-1)+start(y.low)[1],deltat(y.low)*(end(y.low)[2]-1)+end(y.low)[1],deltat(y.low)))),
                c(coredata(y.high), rev(coredata(y.low))), col = "lightgray", border=NA, density=10)

            lines(y.for,col="blue",lwd=2)

            if(legend==TRUE){
# Define where to position the legend
                if(mean(c(y.fit,y.for)[1:round(obs.all/3,0)])<(plot.range[2]+plot.range[1])/2){
                    leg.place = "topleft"
                }
                else{
                    leg.place = "bottomleft"
                }

                legend(x=leg.place,
                   legend=c("Series","Fitted values","Point forecast",paste0(int.w*100,"% prediction interval"),"Forecast origin"),
                   col=c("black","purple","blue","darkgrey","red"),
                   lwd=c(1,2,2,3,2),
                   lty=c(1,2,1,2,1))
            }
        }
        else{
            points(y.low,col="darkgrey",lwd=3,pch=4)
            points(y.high,col="darkgrey",lwd=3,pch=4)
            points(y.for,col="blue",lwd=2,pch=4)

            if(legend==TRUE){
# Define where to position the legend
                if(mean(c(y.fit,y.for)[1:round(obs.all/3,0)])<(plot.range[2]+plot.range[1])/2){
                    leg.place = "topleft"
                }
                else{
                    leg.place = "bottomleft"
                }

                legend(x=leg.place,
                   legend=c("Series","Fitted values","Point forecast",paste0(int.w*100,"% prediction interval"),"Forecast origin"),
                   col=c("black","purple","blue","darkgrey","red"),
                   lwd=c(1,2,2,3,2),
                   lty=c(1,2,NA,NA,1),
                   pch=c(NA,NA,4,4,NA))
            }
        }
    }
    else{
        if(h>1){
            lines(y.for,col="blue",lwd=2)

            if(legend==TRUE){
# Define where to position the legend
                if(mean(c(y.fit,y.for)[1:round(obs.all/3,0)])<(plot.range[2]+plot.range[1])/2){
                    leg.place = "topleft"
                }
                else{
                    leg.place = "bottomleft"
                }

                legend(x=leg.place,
                   legend=c("Series","Fitted values","Point forecast","Forecast origin"),
                   col=c("black","purple","blue","red"),
                   lwd=c(1,2,2,2),
                   lty=c(1,2,1,1))
            }
        }
        else{
            points(y.for,col="blue",lwd=2,pch=4)
            if(legend==TRUE){
# Define where to position the legend
                if(mean(c(y.fit,y.for)[1:round(obs.all/3,0)])<(plot.range[2]+plot.range[1])/2){
                    leg.place = "topleft"
                }
                else{
                    leg.place = "bottomleft"
                }

                legend(x=leg.place,
                   legend=c("Series","Fitted values","Point forecast","Forecast origin"),
                   col=c("black","purple","blue","red"),
                   lwd=c(1,2,2,2),
                   lty=c(1,2,NA,1),
                   pch=c(NA,NA,4,NA))
            }
        }
    }

# Draw the line that divides the series into the "in-sample" and "holdout"
    abline(v=deltat(y.fit)*(end(y.fit)[2]-1)+end(y.fit)[1],col="red",lwd=2)

# Revert par to default parameters
    par(mfrow=c(1,1), mar=c(5,4,4,2))
}

return(list(persistence=vecg,phi=phi,states=matxt,fitted=y.fit,forecast=y.for,residuals=errors,errors=errors.mat,x=data,ICs=ICs,CF=CF.objective,FI=FI))
}