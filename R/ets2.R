ets2 <- function(data, model="ZZZ", persistence=NULL, phi=NULL,
                 bounds=c("usual","admissible"),
                 initial=NULL, initial.season=NULL,
                 IC=c("AIC","AICc","BIC"), trace=FALSE, CF.type=c("GV","TLV","TV"),
                 intervals=FALSE, int.w=0.95, xreg=NULL,
                 holdout=FALSE, h=10, silent=FALSE, legend=TRUE,
                 ...){

# Start measuring the time of calculations
    start.time <- Sys.time();

    bounds <- bounds[1];
    IC <- IC[1];
    CF.type <- CF.type[1];

    if(CF.type!="GV" & CF.type!="TLV" & CF.type!="TV"){
        message(paste0("The strange Cost Function is chosen: ",CF.type));
        message("Switching to 'GV'");
        CF.type <- "GV";
    }

# If chosen model is "AAdN" or anything like that, we are taking the appropriate values
    if(nchar(model)==4){
        error.type <- substring(model,1,1);
        trend.type <- substring(model,2,2);
        season.type <- substring(model,4,4);
        damped <- TRUE;
        if(substring(model,3,3)!="d"){
            message(paste0("You have defined a strange model: ",model));
            sowhat(model);
            model <- paste0(error.type,trend.type,"d",season.type);
        }
    }
    else if(nchar(model)==3){
        error.type <- substring(model,1,1);
        trend.type <- substring(model,2,2);
        season.type <- substring(model,3,3);
        phi <- 1;
        damped <- FALSE;
    }
    else{
        message(paste0("You have defined a strange model: ",model));
        sowhat(model);
        model <- "ZZZ";

        error.type <- "Z";
        trend.type <- "Z";
        season.type <- "Z";
        damped <- TRUE;
    }

# Define obs.all, the overal number of observations (in-sample + holdout)
    obs.all <- length(data) + (1 - holdout)*h;

# Define obs, the number of observations of in-sample
    obs <- length(data) - holdout*h;

# Define the actual values
    y <- coredata(data);

# Check if the data is vector
    if(!is.numeric(data) & !is.ts(data)){
        stop("The provided data is not a vector or ts object! Can't build any model!", call.=FALSE);
    }

# Check if the data is ts-object
    if(!is.ts(data)){
        message("The provided data is not ts object. Only non-seasonal models are available.");
        season.type <- "N";
    }
    seas.freq <- frequency(data);

### Check all the parameters for the possible errors.
    if(!is.null(persistence) & length(persistence)>3){
        message("The length of persistence vector is wrong! It should not be greater than 3.");
        message("Changing to the estimation of persistence vector values.");
        persistence <- NULL;
    }

### Check the error type
    if(error.type!="Z" & error.type!="A" & error.type!="M"){
        message("Wrong error type! Should be 'Z', 'A' or 'M'.");
        message("Changing to 'Z'");
        error.type <- "Z";
    }
# The number of components of ETS
        n.components <- 1;
# The lag of components (needed mainly for the seasonal models)
        lags <- 1;
# The names of the state vector components
        component.names <- "level";

### Check the trend type
    if(trend.type!="Z" & trend.type!="N" & trend.type!="A" & trend.type!="M"){
        message("Wrong trend type! Should be 'Z', 'N', 'A' or 'M'.");
        message("Changing to 'Z'");
        trend.type <- "Z";
    }

    if(trend.type!="N"){
        n.components <- n.components + 1;
        lags <- c(lags,1);
        component.names <- c(component.names,"trend");
        trend.component=TRUE;
    }
    else{
        trend.component=FALSE;
    }

### Check the seasonaity type
    if((season.type=="Z" | season.type=="A" | season.type=="M" ) & seas.freq==1){
        message("Cannot build the seasonal model on the data with the frequency 1.");
        message(paste0("Switching to non-seasonal model: ETS(",error.type,",",trend.type,",N)"));
        season.type <- "N";
    }
    if(season.type!="Z" & season.type!="N" & season.type!="A" & season.type!="M"){
        message("Wrong seasonality type! Should be 'Z', 'N', 'A' or 'M'.");
        if(seas.freq==1){
            message("Data is non-seasonal. Changing seasonal component to 'N'");
            season.type <- "N";
        }
        else{
            message("Changing to 'Z'");
            season.type <- "Z";
        }
    }

    if(season.type!="N"){
        n.components <- n.components + 1;
        lags <- c(lags,seas.freq);
        component.names <- c(component.names,"seasonality");
        seasonal.component <- TRUE;
    }
    else{
        seasonal.component <- FALSE;
# In the cases of non-seasonal models built on seasonal data seasonality is equal to 1
        seas.freq <- 1;
    }

    if(any(y<=0)){
        if(error.type=="M"){
            message("Can't apply multiplicative model to non-positive data. Switching error to 'A'");
            error.type <- "A";
        }
        if(trend.type=="M"){
            message("Can't apply multiplicative model to non-positive data. Switching trend to 'A'");
            trend.type <- "A";
        }
        if(season.type=="M"){
            message("Can't apply multiplicative model to non-positive data. Switching seasonal to 'A'");
            season.type <- "A";
        }
    }

### Check the length of initial and persistence vectors
# Check the persistence vector length
    if(!is.null(persistence)){
        if(n.components != length(persistence)){
            message("The length of persistence vector does not correspond to the chosen model!");
            message("Values will be estimated");
            persistence <- NULL;
        }
    }

# Check the inital vector length
    if(!is.null(initial)){
        if(length(initial)>2){
            message("The length of the initial value is wrong! It should not be greater than 2.");
            message("Values of initial vector will be estimated.");
            initial <- NULL;
        }
        if((n.components - seasonal.component)!=length(initial)){
            message("The length of initial state vector does not correspond to the chosen model!");
            message("Values of initial vector will be estimated.");
            initial <- NULL;
        }
    }

# Check the seasonal inital vector length
    if(!is.null(initial.season)){
        if(seas.freq!=length(initial.season)){
            message("The length of seasonal initial states does not correspond to the chosen frequency!");
            message("Values of initial vector will be estimated.");
            initial.season <- NULL;
        }
    }

### Create the essential matrices and vectors
# Create the matrix of state vectors.
# First seas.freq rows are initial values, excluded from estimation.
    mat.xt <- matrix(NA,nrow=(obs+seas.freq),ncol=n.components);
    colnames(mat.xt) <- component.names;

# If the persistence vector is provided, use it
    if(!is.null(persistence)){
        vec.g <- persistence;
        estimate.persistence <- FALSE;
    }
    else{
        vec.g <- rep(0.3,n.components);
        estimate.persistence <- TRUE;
    }

# If phi is not provided, mark that
    if(damped==TRUE){
        if(is.null(phi)){
            estimate.phi = TRUE;
            phi <- 0.95;
        }
        else{
            if(phi<0 | phi>2){
                message("Damping parameter should lie in (0, 2) region.");
                message("Changing to the estimation of phi.");
                phi <- 0.95;
                estimate.phi = TRUE;
            }
            else{
                estimate.phi = FALSE;
            }
        }
    }
    else{
        estimate.phi = FALSE;
    }

# If initial values are provided, write them down
    if(!is.null(initial)){
        mat.xt[1:seas.freq,1:(n.components - seasonal.component)] <- rep(initial,each=seas.freq);
        estimate.initial <- FALSE;
    }
    else{
        if(trend.type=="N"){
            mat.xt[1:seas.freq,1] <- rep(mean(y[1:min(12,obs)]),times=seas.freq);
        }
        else if(trend.type=="A"){
            slope <- cov(y[1:min(12,obs)],c(1:min(12,obs)))/var(c(1:min(12,obs)));
            intercept <- mean(y[1:min(12,obs)]) - slope * (mean(c(1:min(12,obs))) - 1);
            mat.xt[1:seas.freq,1] <- rep(intercept,seas.freq);
            mat.xt[1:seas.freq,2] <- rep(slope,seas.freq);
        }
        else if(trend.type=="M"){
            mat.xt[1:seas.freq,1] <- rep(mean(data[1:min(12,obs)]),seas.freq);
            mat.xt[1:seas.freq,2] <- rep(1,seas.freq);
        }
        estimate.initial <- TRUE;
    }

# If the seasonal model is chosen and initials are provided, fill in the first "seas.freq" values of seasonal component.
    if(seasonal.component==TRUE & !is.null(initial.season)){
        mat.xt[1:seas.freq,n.components] <- initial.season;
        estimate.initial.season <- FALSE;
    }
    else if(seasonal.component==TRUE & is.null(initial.season)){
        estimate.initial.season <- TRUE;
        if(season.type=="A"){
            mat.xt[1:seas.freq,n.components] <- decompose(ts(data[1:obs],frequency=seas.freq),type="additive")$seasonal[1:seas.freq];
        }
        else if(season.type=="M"){
            mat.xt[1:seas.freq,n.components] <- decompose(ts(data[1:obs],frequency=seas.freq),type="multiplicative")$seasonal[1:seas.freq];
        }
    }
    else{
        estimate.initial.season <- FALSE;
    }

# Vectors of fitted data and errors
    y.fit <- rep(NA,obs);
    errors <- rep(NA,obs);

# If we use trace, define matrix of errors.
    if(trace==TRUE){
        mat.error <- matrix(NA,nrow=obs,ncol=h);
    }
    else{
        mat.error <- matrix(NA,nrow=obs,ncol=1);
    }

# Function calculates the power of matrix
matrix.power <- function(A, n) {
    if(n==0){
        B <- diag(nrow(A));
    }
    else{
        B <- diag(nrow(A));
        for(i in 1:n){
            B <- A %*% B;
        }
    }
    return(B);
}

# Function returns transition matrix F and measurement matrix w depending on chosen model
mat.ets <- function(trend.component,seasonal.component,phi){
    if(seasonal.component==FALSE){
        if(trend.component==FALSE){
            mat.F <- matrix(1,1,1);
            mat.w <- 1;
        }
        else{
            mat.F <- matrix(c(1,0,phi,phi),2,2);
            mat.w <- c(1,phi);
        }
    }
    else{
        if(trend.component==FALSE){
            mat.F <- matrix(c(1,0,0,1),2,2);
            mat.w <- c(1,1);
        }
        else{
            mat.F <- matrix(c(1,0,0,phi,phi,0,0,0,1),3,3);
            mat.w <- c(1,phi,1);
        }
    }
    return(list(mat.F=mat.F,mat.w=mat.w));
}

# Function fills in the initial values of mat.xt
estim.values <- function(mat.xt,vec.g,phi,C){

# If the persistence vector is provided, use it
    if(estimate.persistence==TRUE){
        vec.g <- C[1:n.components];
    }

# If phi was not provided, use it in estimate
    if(estimate.phi==TRUE){
        phi <- C[n.components*estimate.persistence + 1];
    }
    
# If initial values are not provided, write them in the matrix
    if(estimate.initial==TRUE){
        mat.xt[1:seas.freq,1:(n.components - seasonal.component)] <- rep(C[(n.components*estimate.persistence + estimate.phi + 1):(n.components*estimate.persistence + estimate.phi + n.components - seasonal.component)],each=seas.freq);
    }
    
# If the seasonal model is chosen and initials are provided, fill in the first "seas.freq" values of seasonal component.
    if(estimate.initial.season==TRUE){
        mat.xt[1:seas.freq,n.components] <- C[(n.components*estimate.persistence + estimate.phi + (n.components - seasonal.component)*estimate.initial + 1):(n.components*estimate.persistence + estimate.phi + (n.components - seasonal.component)*estimate.initial + seas.freq)];
    }

    return(list(vec.g=vec.g,phi=phi,mat.xt=mat.xt));
}

# Function returns the value of r based on the type of the model
r.value <- function(xt,mat.w){
    if(error.type=="A"){
        if(trend.type=="N"){
            if(season.type!="M"){
                r <- rep(1,n.components);
            }
            else{
                r <- c(1/xt[2], 1/xt[1]);
            }
        }
        else if(trend.type=="A"){
            if(season.type!="M"){
                r <- rep(1,n.components);
            }
            else{
                r <- c(1/xt[3], 1/xt[3], 1/(xt[1] + mat.w[2] * xt[2]));
            }
        }
        else if(trend.type=="M"){
            if(season.type=="N"){
                r <- c(1, 1/xt[1]);
            }
            else if(season.type=="A"){
                r <- c(1, 1/xt[1], 1);
            }
            else if(season.type=="M"){
                r <- c(1/xt[3], 1/(xt[1]*xt[3]), 1/(xt[1]*xt[2]^mat.w[2]));
            }
        }
    }
    else{
        if(trend.type=="N"){
            if(season.type=="N"){
                r <- xt[1];
            }
            else if(season.type=="A"){
                r <- rep(xt[1] + xt[2],n.components);
            }
            else if(season.type=="M"){
                r <- c(x[1], x[2]);
            }
        }
        else if(trend.type=="A"){
            if(season.type!="M"){
                r <- rep(mat.w %*% xt,n.components);
            }
            else if(season.type=="M"){
                r <- c(xt[1] + xt[2]*mat.w[2], xt[1] + xt[2]*mat.w[2], xt[3]);
            }
        }
        else if(trend.type=="M"){
            if(season.type=="N"){
                r <- c(xt[1]*xt[2]^mat.w[2],xt[2]^mat.w[2]);
            }
            else if(season.type=="A"){
                r <- c(xt[1]*xt[2]^mat.w[2] + xt[3], (xt[1]*xt[2]^mat.w[2] + xt[3])/xt[1], xt[1]*xt[2]^mat.w[2] + xt[3])
            }
            else if(season.type=="M"){
                r <- c(xt[1]*xt[2]^mat.w[2], xt[2]^mat.w[2], xt[3]);
            }
        }
    }
    return(r);
}

# Function returns error depending on error type
error <- function(y.act,y.est){
    if(error.type=="A"){
        error <- y.act - y.est;
    }
    else{
        error <- (y.act - y.est) / y.est;
    }
    return(error);
}

# Function constructs ETS with given values
fit.ets <- function(y,mat.xt,mat.F,mat.w,vec.g){
## Cases of ZNN, ZAN, ZNA, ZAA
    if(trend.type!="M" & season.type!="M"){
        for(j in (seas.freq+1):(obs+seas.freq)){
            vec.r <- r.value(mat.xt[cbind((j-lags),c(1:n.components))],mat.w);
            y.fit[j-seas.freq] <- mat.w %*% mat.xt[cbind((j-lags),c(1:n.components))];
            errors[j-seas.freq] <- error(y[j-seas.freq],y.fit[j-seas.freq]);
            mat.xt[j,] <- mat.F %*% mat.xt[cbind((j-lags),c(1:n.components))] + vec.g * errors[j-seas.freq] * vec.r;
        }
    }
## Cases of ZMN, ZNM, ZMM
    else if(trend.type!="A" & season.type!="A"){
        if(trend.type!="N"){
            for(j in (seas.freq+1):(obs+seas.freq)){
                vec.r <- r.value(mat.xt[cbind((j-lags),c(1:n.components))],mat.w);
                y.fit[j-seas.freq] <- exp(mat.w %*% log(mat.xt[cbind((j-lags),c(1:n.components))]));
                errors[j-seas.freq] <- error(y[j-seas.freq],y.fit[j-seas.freq]);
                mat.xt[j,] <- exp(mat.F %*% log(mat.xt[cbind((j-lags),c(1:n.components))])) + vec.g * errors[j-seas.freq] * vec.r;
            }
        }
    }
## Case of ZAM
    else if(trend.type=="A" & season.type=="M"){
        for(j in (seas.freq+1):(obs+seas.freq)){
            vec.r <- r.value(mat.xt[cbind((j-lags),c(1:n.components))],mat.w);
            y.fit[j-seas.freq] <- mat.w[1:(n.components-1)] %*% mat.xt[cbind((j-lags[1:(n.components-1)]),c(1:(n.components-1)))] * mat.xt[j-lags[n.components],n.components];
            errors[j-seas.freq] <- error(y[j-seas.freq],y.fit[j-seas.freq]);
            mat.xt[j,] <- mat.F %*% mat.xt[cbind((j-lags),c(1:n.components))] + vec.g * errors[j-seas.freq] * vec.r;
        }
    }
## Case of ZMA
    else if(trend.type=="M" & season.type=="A"){
        for(j in (seas.freq+1):(obs+seas.freq)){
            vec.r <- r.value(mat.xt[cbind((j-lags),c(1:n.components))],mat.w);
            y.fit[j-seas.freq] <- exp(mat.w[1:(n.components-1)] %*% log(mat.xt[cbind((j-lags[1:(n.components-1)]),c(1:(n.components-1)))])) + mat.xt[j-lags[n.components],n.components];
            errors[j-seas.freq] <- error(y[j-seas.freq],y.fit[j-seas.freq]);
            mat.xt[j,] <- Re(exp(mat.F %*% log(as.complex(mat.xt[cbind((j-lags),c(1:n.components))])))) + vec.g * errors[j-seas.freq] * vec.r;
        }
    }
    return(list(mat.xt=mat.xt, errors=errors, y.fit=y.fit));
}

# Function makes forecast from the specified point to the specified horizon
forec.ets <- function(xt,mat.F,mat.w,vec.g,h=1,j){
    
    y.for <- rep(NA,h);

## Form vector of non-seasonal and vector of seasonal components
    if(season.type!="N"){
	    season.xt <- rep(xt[,n.components],times=ceiling(h/seas.freq));
	    xt <- matrix(xt[nrow(xt),1:(n.components-1)],(n.components-1),1);
	    mat.w <- matrix(mat.w[1:(n.components-1)],1,(n.components-1));
	    mat.F <- matrix(mat.F[1:(n.components-1),1:(n.components-1)],(n.components-1),(n.components-1));
    }
    else{
	    xt <- matrix(xt,n.components,1);
	    mat.w <- matrix(mat.w,1,n.components);
	    mat.F <- matrix(mat.F,n.components,n.components);
    }

## Cases of ZNN, ZAN
    if(trend.type!="M"){
        for(i in 1:h){
            y.for[i] <- mat.w %*% matrix.power(mat.F,(i-1)) %*% xt;
        }
    }
## Cases of ZMN
    else if(trend.type=="M"){
        for(i in 1:h){
            y.for[i] <- exp(mat.w %*% matrix.power(mat.F,(i-1)) %*% log(xt));
        }
    }
## Case of ZZA
    if(season.type=="A"){
        y.for <- y.for + season.xt[1:h];
    }
## Case of ZZM
    else if(season.type=="M"){
        y.for <- y.for * season.xt[1:h];
    }

    return(y.for)
}

# The function constructs the forecasts from each point and
# estimates errors. Needed for trace. SHOULD CALL forec.ets
errors.ets <- function(mat.xt,mat.F,mat.w,vec.g,trace){
    if(trace==TRUE){
	    residuals <- matrix(NA,obs,h);
	    for(j in (1+seas.freq):(obs+seas.freq)){
		hh <- min(h,(obs+seas.freq)-j+1);
	        residuals[(j-seas.freq),1:hh] <- y[(j-seas.freq):(j-seas.freq+hh-1)] - forec.ets(xt=mat.xt[(j-seas.freq):(j-1),],mat.F,mat.w,vec.g,h=hh);
	    }
    }
    else{
	    residuals <- matrix(NA,obs,1);
	    for(j in (1+seas.freq):(obs+seas.freq)){
	        residuals[(j-seas.freq),1] <- y[j-seas.freq] - forec.ets(xt=mat.xt[(j-seas.freq):(j-1),],mat.F,mat.w,vec.g,h=1);
	    }
    }
    return(residuals);
}

# Cost function for ETS
CF <- function(C){
    init.ets <- estim.values(mat.xt,vec.g,phi,C);
    vec.g <- init.ets$vec.g;
    phi <- init.ets$phi;
    mat.xt <- init.ets$mat.xt;

    matrices <- mat.ets(trend.component,seasonal.component,phi);
    mat.F <- matrices$mat.F;
    mat.w <- matrices$mat.w;

    fitting <- fit.ets(y,mat.xt,mat.F,mat.w,vec.g);
    mat.xt <- fitting$mat.xt;
    y.fit <- fitting$y.fit;
    errors <- fitting$errors;

    errors.mat <- errors.ets(mat.xt,mat.F,mat.w,vec.g,trace);

    if(trace==TRUE){
        if(CF.type=="GV"){
            errors.mat <- errors.mat[!is.na(errors.mat[,h]),];
            CF.res <- det(t(errors.mat) %*% (errors.mat) / errors.mat.obs);
        }
        else if(CF.type=="TLV"){
            CF.res <- exp(sum(log(colMeans(errors.mat^2,na.rm=TRUE))));
        }
        else if(CF.type=="TV"){
            CF.res <- sum(colMeans(errors.mat^2,na.rm=TRUE));
        }
    }
    else{
        CF.res <- mean(errors.mat^2,na.rm=TRUE);
    }

    return(CF.res);
}

# Function returns interval forecasts.
int.ets <- function(xt,mat.F,mat.w,vec.g,h=1){
    y.lo <- rep(NA,h);
    y.up <- y.lo;
    return(list(y.up=y.up,y.lo=y.lo));
}

MASE.lvl <- function(a,f,round=3){
    # This function calculates Mean Absolute Scaled Error using level of the actual series
    # a - actual values,
    # f - forecasted values.
    result <- round(sum(abs(a-f),na.rm=TRUE)/sum(a),round);
    return(result);
}

MASE <- function(a,f,scale,round=3){
    # This function calculates Mean Absolute Scaled Error as in Hyndman & Koehler, 2006
    # a - actual values,
    # f - forecasted values.
    # scale - the measure to scale errors with. Usually - MAE of in-sample.
    result <- round(mean(abs(a-f),na.rm=TRUE)/scale,round);
    return(result);
}

hin.constrains.usual <- function(C){
    p <- NA;
    d <- NA;
    i <- NA;
    i.s <- NA;

### Constrains for persistence vector
    if(estimate.persistence==TRUE){
# smoothing parameters constrains
        p <- rep(NA,n.components*2);

# alpha restrictions (0, 1)
        p[1] <- C[1];
        p[2] <- 1 - C[1];
        if(trend.component==TRUE){
# beta restrictions (0, alpha)
            p[3] <- C[2];
            p[4] <- C[1] - C[2];
            if(seasonal.component==TRUE){
# gamma restrictions (0, 1 - alpha)
                p[5] <- C[3];
                p[6] <- 1 - C[1] - C[3];
            }
        }
        else{
            if(seasonal.component==TRUE){
# gamma restrictions (0, 1 - alpha)
                p[3] <- C[2];
                p[4] <- 1 - C[1] - C[2];
            }
        }
    }

### Constrains for damping parameter (0, 1)
    if(estimate.phi==TRUE){
        d <- rep(NA,2);
        d[1] <- C[n.components*estimate.persistence + 1];
        d[2] <- 1 - C[n.components*estimate.persistence + 1];
    }

### Constrains on initial parameters (-Inf, Inf)
    if(estimate.initial){
        i <- 1;
    }

### Constrains on initial seasonal parameters (-Inf, Inf) for additive and (0.1, 1.9) for multiplicative
    if(estimate.initial.season==TRUE){
        i.s <- rep(NA,(seas.freq*2));
        n <- n.components*estimate.persistence + estimate.phi + (n.components - seasonal.component)*estimate.initial;
            if(season.type=="A"){
                i.s <- 1;
            }
            else{
                i.s[1:seas.freq] <- C[(n + 1):(n + seas.freq)] - 0.1;
                i.s[(seas.freq+1):(2*seas.freq)] <- 1.9 - C[(n + 1):(n + seas.freq)];
            }
    }

    constrains <- c(p,d,i,i.s);
    constrains <- constrains[!is.na(constrains)];

    return(constrains);
}

# For the automatic model selection. If data has negative or zero values,
# exclude multiplicative models
#if(any(y<=0)){
#    errors.pool <- c("N","A");
#    trends.pool <- c("N","A");
#    season.pool <- c("N","A");
#}

# Fill in the vector of initial values used in estimation
    if(estimate.persistence==TRUE | estimate.phi==TRUE | estimate.initial==TRUE | estimate.initial.season==TRUE){
        C <- NA;
        C.lower <- NA;
        C.upper <- NA;

        if(estimate.persistence==TRUE){
            C <- c(C,vec.g);
            C.lower <- c(C.lower,rep(0,length(vec.g)));
            C.upper <- c(C.upper,rep(1,length(vec.g)));
        }
        if(estimate.phi==TRUE){
            C <- c(C,phi);
            C.lower <- c(C.lower,0);
            C.upper <- c(C.upper,1.2);
        }
        if(estimate.initial==TRUE){
            C <- c(C,mat.xt[seas.freq,1:(n.components - seasonal.component)]);
            if(trend.type!="M"){
                C.lower <- c(C.lower,rep(-Inf,(n.components - seasonal.component)));
                C.upper <- c(C.upper,rep(Inf,(n.components - seasonal.component)));
            }
            else{
                C.lower <- c(C.lower,1,0.01);
                C.upper <- c(C.upper,Inf,2);
            }
        }
        if(estimate.initial.season==TRUE){
            C <- c(C,mat.xt[1:seas.freq,n.components]);
            if(season.type=="A"){
                C.lower <- c(C.lower,rep(-Inf,seas.freq));
                C.upper <- c(C.upper,rep(Inf,seas.freq));
            }
            else{
                C.lower <- c(C.lower,rep(0,seas.freq));
                C.upper <- c(C.upper,rep(2,seas.freq));
            }
        }
        C <- C[!is.na(C)];
        C.lower <- C.lower[!is.na(C.lower)];
        C.upper <- C.upper[!is.na(C.upper)];

# Number of observations in the mat.error matrix excluding NAs.
        errors.mat.obs <- obs - h + 1;

#library(nloptr);
        res <- cobyla(C, CF, hin=hin.constrains.usual, lower=C.lower, upper=C.upper);
        CF.objective <- res$value;

        init.ets <- estim.values(mat.xt,vec.g,phi,res$par);
        vec.g <- init.ets$vec.g;
        phi <- init.ets$phi;
        mat.xt <- init.ets$mat.xt;
    }

    matrices <- mat.ets(trend.component,seasonal.component,phi);
    mat.F <- matrices$mat.F;
    mat.w <- matrices$mat.w;

    fitting <- fit.ets(y,mat.xt,mat.F,mat.w,vec.g);
    mat.xt <- ts(fitting$mat.xt,start=(time(data)[1] - deltat(data)*seas.freq),frequency=frequency(data));
    y.fit <- ts(fitting$y.fit,start=start(data),frequency=frequency(data));
    errors <- ts(fitting$errors,start=start(data),frequency=frequency(data));

    y.for <- ts(forec.ets(xt=mat.xt[((obs+seas.freq)-seas.freq+1):(obs+seas.freq),],mat.F,mat.w,vec.g,h=h),start=time(data)[obs]+deltat(data),frequency=frequency(data));

    y <- data;

    if(intervals==T){
        y.low <- ts(y.low,start=start(y.for),frequency=frequency(data));
        y.high <- ts(y.high,start=start(y.for),frequency=frequency(data));
        plot.range <- range(min(data,y.fit,y.for,y.low),max(data,y.fit,y.for,y.high));
    }
    else{
        plot.range <- range(min(data,y.fit,y.for),max(data,y.fit,y.for));
    }

    if(estimate.persistence==FALSE & estimate.phi==FALSE & estimate.initial==FALSE & estimate.initial.season==FALSE){
        C <- c(vec.g,phi,initial,initial.season);
        CF.objective <- CF(C);
#        n.param <- n.components + damped + (n.components - seasonal.component) + seas.freq*seasonal.component;
        n.param <- 0;
    }
    else{
        n.param <- n.components*estimate.persistence + estimate.phi + (n.components - seasonal.component)*estimate.initial + seas.freq*estimate.initial.season;
    }

# Information criteria are calculated with the constant part "log(2*pi*exp(1)*h)*obs".
# And it is based on the mean of the sum squared residuals either than sum.
# Hyndman likelihood is: llikelihood <- obs*log(obs*CF.objective);
    llikelihood <- -obs/2 *((h^trace)*log(2*pi*exp(1)) + log(CF.objective));
    AIC.coef <- 2*n.param - 2*llikelihood;
    AICc.coef <- AIC.coef + 2 * n.param * (n.param + 1) / (obs - n.param - 1);
    BIC.coef <- log(obs)*n.param - 2*llikelihood;

    ICs <- c(AIC.coef, AICc.coef, BIC.coef);
    names(ICs) <- c("AIC", "AICc", "BIC");

if(silent==FALSE){
# Print time elapsed on the construction
    print(paste0("Time elapsed: ",round(as.numeric(Sys.time() - start.time,units="secs"),2)," seconds"));
    print(paste0("Model constructed: ",model));
    print(paste0("Persistence vector: ", paste(round(vec.g,3),collapse=", ")));
    if(damped==TRUE){
        print(paste0("Damping parameter: ", round(phi,3)));
    }
    print(paste0("Initial components: ", paste(round(mat.xt[seas.freq,1:(n.components - seasonal.component)],3),collapse=", ")));
    if(seasonal.component==TRUE){
        print(paste0("Initial seasonal components: ", paste(round(mat.xt[1:seas.freq,n.components],3),collapse=", ")));
    }
    print(paste0("Residuals sigma: ",round(sqrt(var(errors)),3)));
    if(trace==TRUE){
        print(paste0("CF type: ",CF.type));
    }
    print(paste0("CF value is: ",round(CF.objective,0)));
    print(paste0("Unbiased log-likelihood: ",round((llikelihood - n.param),0)));
    print(paste0("AIC: ",round(AIC.coef,3),"; AICc: ", round(AICc.coef,3), "; BIC: ", round(BIC.coef,3)));
    if(holdout==T){
        print(paste0("MASE: ",MASE(coredata(data)[(obs+1):obs.all],coredata(y.for),mean(abs(diff(coredata(data)[1:obs]))),round=3)));
        print(paste0("MASE.lvl: ",MASE.lvl(coredata(data)[(obs+1):obs.all],coredata(y.for),round=5)*100,"%"));
    }

    par(mfrow=c(1,1), mar=c(5,3,2,1))
    plot(data,type="l",xlim=range(time(data)[1],time(y.for)[h]),
         ylim=plot.range,xlab="Time", ylab="")
    
    if(intervals==T){
        lines(y.low,col="darkgrey",lwd=3,lty=2);
        lines(y.high,col="darkgrey",lwd=3,lty=2);
# Draw the nice areas between the borders
        polygon(c(seq(deltat(y.high)*(start(y.high)[2]-1)+start(y.high)[1],deltat(y.high)*(end(y.high)[2]-1)+end(y.high)[1],deltat(y.high)),
                  rev(seq(deltat(y.low)*(start(y.low)[2]-1)+start(y.low)[1],deltat(y.low)*(end(y.low)[2]-1)+end(y.low)[1],deltat(y.low)))),
                c(coredata(y.high), rev(coredata(y.low))), col = "lightgray", border=NA, density=10);
        
        lines(y.fit,col="purple",lwd=2,lty=2);
        lines(y.for,col="blue",lwd=2);
        if(legend==TRUE){
# Define where to position the legend
            if(mean(c(y.fit,y.for)[1:round(obs.all/3,0)])<(plot.range[2]+plot.range[1])/2){
                leg.place = "topleft";
            }
            else{
                leg.place = "bottomleft";
            }

            legend(x=leg.place,
                   legend=c("Series","Fitted values","Point forecast",paste0(int.w*100,"% prediction interval"),"Forecast origin"),
                   col=c("black","purple","blue","darkgrey","red"),
                   lwd=c(1,2,2,3,2),
                   lty=c(1,2,1,2,1))
        }
    }
    else{
        lines(y.fit,col="purple",lwd=2,lty=2);
        lines(y.for,col="blue",lwd=2);
        
        if(legend==TRUE){
# Define where to position the legend
            if(mean(c(y.fit,y.for)[1:round(obs.all/3,0)])<(plot.range[2]+plot.range[1])/2){
                leg.place = "topleft";
            }
            else{
                leg.place = "bottomleft";
            }

            legend(x=leg.place,
                   legend=c("Series","Fitted values","Point forecast","Forecast origin"),
                   col=c("black","purple","blue","red"),
                   lwd=c(1,2,2,2),
                   lty=c(1,2,1,1))
        }
    }

# Draw the line that divides the series into the "in-sample" and "holdout"
    abline(v=deltat(y.fit)*(end(y.fit)[2]-1)+end(y.fit)[1],col="red",lwd=2);

# Revert par to default parameters
    par(mfrow=c(1,1), mar=c(5,4,4,2))
}

return(list(persistence=vec.g,phi=phi,states=mat.xt,fitted=y.fit,forecast=y.for,residuals=errors,x=data,ICs=ICs,CF=CF.objective));
}
