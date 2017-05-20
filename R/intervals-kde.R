intervals.kde <- function(model,level=0.95,type=c("multi","single"),
                          kdetype=c("diffusion","SJ","nrd0")){
  # Calculate PIs using KDE
    
  # Arguments checks
  type <- match.arg(type,c("multi","single"))
  kdetype <- match.arg(kdetype,c("diffusion","SJ","nrd0"))
  
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
    # Estimate KCDE
    ehat <- TStools::kdemode(eh[,j],type=kdetype)
    # Scale cdf - seems to be OK!
    kcde <- cumsum(ehat$fd)
    kcde <- kcde/max(kcde)
    # Get requested points
    k <- length(b)
    x <- vector("numeric",k)
    for (i in 1:k){
      idx <- order(abs(kcde-b[i]))[1:2]
      x[i] <- approx(kcde[idx],ehat$xd[idx],xout=b[i])$y
    }
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