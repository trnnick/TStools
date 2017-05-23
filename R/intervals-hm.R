intervals.hm <- function(model, level=0.95, centre=TRUE){
    if(length(class(model))==1){
        if(class(model)!="smooth"){
            stop("Sorry, but we need smooth object for this stuff. And your model is not it...");
        }
    }
    else{
        stop("Sorry, but we need smooth object for this stuff. And your model is not it...");
    }
    
    errors <- model$errors;
    h <- ncol(errors);
    errors <- errors[-c(1:(h-1)),];
    if(centre){
        centre <- colMeans(errors);
        errors <- errors - matrix(centre,nrow=nrow(errors),ncol=h,byrow=TRUE);
    }
    else{
        centre <- 0;
    }
    
    if(level>1){
        if(level>100){
            warning("The specified level is meaningless. Swithing to 95%.");
            level <- 95;
        }
        level <- level/100;
    }
    
    hsmN <- gamma(0.75)*pi^(-0.5)*2^(-0.75);
    quantValues <- qnorm(c((1-level)/2,(1+level)/2));
    hmValues <- apply(errors,2,hm,C=0);
    
    lower <- (model$forecast + centre + quantValues[1] * Im(hmValues)^2 / hsmN^2);
    upper <- (model$forecast + centre + quantValues[2] * Re(hmValues)^2 / hsmN^2);

    names(hmValues) <- paste0("h",1:h);
    
    return(list(lower=lower,upper=upper,hmValues=hmValues));
}