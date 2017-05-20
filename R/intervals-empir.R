intervals.empir <- function(model,level=0.95,type=c("multi","single"),quantiletype=7){
  # Calculate empirical PIs 
    
  # Arguments checks
  type <- match.arg(type,c("multi","single"))
  
  if (!any(class(model)=="smooth")){
    stop("Argument model should be of class `smooth'");
  }
  
  if (any(level>1) | length(level)>1){
    stop("Level should be a positive scalar < 1.")
  }
  
  # PI level
  b <- (1-level)/2 + c(0,1)*level
  
  # Get error matrix from model
  eh <- model$errors
  h <- dim(eh)[2]       # Forecast horizon
  eh <- eh[-(1:h),]
  
  # Calculate intervals
  lower <- upper <- vector("numeric",h)
  # Use only t+1 errors or up to t+h
  if (type == "single"){
    J <- 1
  } else {
    J <- h
  }
    
  for (j in 1:J){
    x <- as.numeric(quantile(eh[,j],b,type=quantiletype))
    lower[j] <- x[1]
    upper[j] <- x[2]
  }
  
  if (type == "single"){
    upper[] <- upper[1]*sqrt(1:h)
    lower[] <- lower[1]*sqrt(1:h)
  }
  
  # Add to forecasts
  lower <- model$forecast + lower
  upper <- model$forecast + upper

  return(list(lower=lower,upper=upper))
  
}