trendtest <- function(y,extract=c("FALSE","TRUE"),type=c("aicc","cs"),mta=c(FALSE,TRUE)){
    # Test for trend
    #
    # Inputs:
    #   y         Time series (must be ts class)
    #   extract   Extract trend using CMA. Test is un on CMA.
    #   type      Type of trend test. "aicc" uses AICc of ETS models. "cs" uses Cox-Stuart
    #   mta       Enhance test with Multiple Temporal Aggregation
    # 
    # Output:
    #   H         Result: 1 there is evidence of trend, 0 no evidence of trend
    #
    # 
    # Example:
    #   trend.test(AirPassengers,TRUE)
    #  
    # Nikolaos Kourentzes, 2017 <nikolaos@kourentzes.com>
    
    # Defaults
    extract <- extract[1]
    type <- match.arg(type,c("aicc","cs"))
    mta <- mta[1]
    
    if (class(y) != "ts"){
        stop("Input y must be of class ts")
    }
    
    # Extract trend
    if (extract == TRUE){
        y.in <- cmav(y) 
    } else {
        y.in <- y
    }
    
    # MTA if needed
    if (mta == TRUE){
        f <- frequency(y) 
        yaggr <- MAPA::tsaggr(y.in,1:f)$out
        # Remove levels with very few observations, but always keep dissagregate
        idx <- unlist(lapply(yaggr,function(x){length(x)>=5}))
        idx[1] <- TRUE
        yaggr <- yaggr[idx]
    } else {
        yaggr <- MAPA::tsaggr(y.in,1)$out
    }
        
    if (type == "cs"){
        H <- unlist(lapply(yaggr,test.cs))
    } else {
        H <- unlist(lapply(yaggr,test.aicc))
    }
    
    H <- median(H)>0.5
    
    return(H)
    
}

test.cs <- function(y){
    H <- coxstuart(y)$H
    return(H)
}

test.aicc <- function(y){
    # Construct and AICc test
    aicc1 <- tryCatch({ets(y,model="ANN")$aicc}, error = function(e){10^5})
    aicc2 <- tryCatch({ets(y,model="AAN",damped=FALSE)$aicc}, error = function(e){10^5})
    aicc3 <- tryCatch({ets(y,model="AAN",damped=TRUE)$aicc}, error = function(e){10^5})
    aicc <- c(aicc1,aicc2,aicc3)
    H.aicc <- which(aicc == min(aicc))[1] != 1
}