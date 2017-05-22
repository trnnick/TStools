intervals.qr <- function(model, level=0.95){
    if(length(class(model))==1){
        if(class(model)!="smooth"){
            stop("Sorry, but we need smooth object for this stuff. And your model is not it...");
        }
    }
    else{
        stop("Sorry, but we need smooth object for this stuff. And your model is not it...");
    }
    
    QFunction <- function(C,alpha){
        values <- C[1] + C[2]*t + C[3]*t^2;
        CF <- (1-alpha)*sum(abs(errors[errors<values]-values[errors<values]))+alpha*sum(abs(errors[errors>=values]-values[errors>=values]));
        return(CF);
    }
    
    errors <- model$errors;
    h <- ncol(errors);
    errors <- errors[-c(1:(h-1)),];
    obs <- nrow(errors);
    
    if(level>1){
        if(level>100){
            warning("The specified level is meaningless. Swithing to 95%.");
            level <- 95;
        }
        level <- level/100;
    }
    t <- matrix(c(1:h),byrow=TRUE,nrow=obs,ncol=h);
    C <- rep(0.5,3);
    res1 <- nlminb(C,QFunction,alpha=(1-level)/2);
    res2 <- nlminb(C,QFunction,alpha=(1+level)/2);

    lower <- (model$forecast + res1$par[1] + res1$par[2]*t[1,] + res1$par[3]*t[1,]^2);
    upper <- (model$forecast + res2$par[1] + res2$par[2]*t[1,] + res2$par[3]*t[1,]^2);
    
    return(list(lower=lower,upper=upper));
}

intervals.qr2 <- function(model, level=0.95){
    if(length(class(model))==1){
        if(class(model)!="smooth"){
            stop("Sorry, but we need smooth object for this stuff. And your model is not it...");
        }
    }
    else{
        stop("Sorry, but we need smooth object for this stuff. And your model is not it...");
    }
    
    QFunction <- function(C,alpha){
        values <- C[1]*t^C[2];
        CF <- (1-alpha)*sum(abs(errors[errors<values]-values[errors<values]))+alpha*sum(abs(errors[errors>=values]-values[errors>=values]));
        return(CF);
    }
    
    errors <- model$errors;
    h <- ncol(errors);
    errors <- errors[-c(1:(h-1)),];
    obs <- nrow(errors);
    
    if(level>1){
        if(level>100){
            warning("The specified level is meaningless. Swithing to 95%.");
            level <- 95;
        }
        level <- level/100;
    }
    t <- matrix(c(1:h),byrow=TRUE,nrow=obs,ncol=h);
    C <- rep(1,2);
    res1 <- nlminb(C,QFunction,alpha=(1-level)/2);
    res2 <- nlminb(C,QFunction,alpha=(1+level)/2);
    
    lower <- (model$forecast + res1$par[1]*t[1,]^res1$par[2]);
    upper <- (model$forecast + res2$par[1]*t[1,]^res2$par[2]);
    
    return(list(lower=lower,upper=upper));
}