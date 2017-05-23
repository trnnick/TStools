intervals.sd <- function(model, level=0.95, centre=TRUE){
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
    
    lower <- rep(NA,h);
    upper <- rep(NA,h);
    
    sigma <- rep(NA,h);
    
    quantValues <- qnorm(c((1-level)/2,(1+level)/2));
    
    for(i in 1:h){
        sigma[i] <- sqrt(mean((errors[,i])^2));
        
        lower[i] <- sigma[i] * quantValues[1];
        upper[i] <- sigma[i] * quantValues[2];
    }
    
    lower <- model$forecast + centre + lower;
    upper <- model$forecast + centre + upper;
    
    sigma <- sigma;
    names(sigma) <- paste0("h",1:h);
    
    return(list(lower=lower,upper=upper,sigma=sigma));
}