ro <- function(data,h=1,origins=1,call,value=NULL,all=FALSE){
# Function makes Rolling Origin for the data using the call
#    Copyright (C) 2015  Yves Sagaert & Ivan Svetunkov
  if(!is.null(value)){
    if(substring(value,1,1)!="$"){
      value <- paste0("$",value)
    }
  }

# Check if the data is vector or ts
  if(!is.numeric(data) & !is.ts(data)){
    stop("The provided data is not a vector or ts object! Can't work with it!", call.=FALSE);
  }
  
# Write down the entire dataset y and the original horizon
  y <- data
  hh <- h
  obs <- length(y)
  in.sample <- obs - origins

# Make matrix to store the results origins in colnames,
#   different rows are entire horizons forecasted
  forecasts <- matrix(NA,nrow=h,ncol=origins)
  colnames(forecasts) <- paste0("origin",c(1:origins))
  rownames(forecasts) <- paste0("h",c(1:h))
  actuals <- forecasts

  cat(paste0("Origins done:  "))

  if(all==TRUE){
    for (i in 1:origins){
# Adjust forecasting horizon to not exeed the sample size
      h <- min(hh,obs - (in.sample+i-1))
# Make the in-sample
      data <- y[1:(in.sample+i-1)]
# Evaluate the call string and save to object o.m.
      o.m <- eval(parse(text=call))
# Save the forecast and the corresponding actuals in matrices
      forecasts[1:h,i] <- eval(parse(text=paste0("o.m",value)))
      actuals[1:h,i] <- y[(in.sample+i):(in.sample+i+h-1)]
      cat(paste(rep("\b",nchar(i)),collapse=""))
      cat(i)
    }
  }
  else{
    for (i in 1:origins){
# Make the in-sample
      data <- y[1:(in.sample+i-h)]
# Evaluate the call string and save to object o.m.
      o.m <- eval(parse(text=call))
# Save the forecast and the corresponding actuals in matrices
      forecasts[,i] <- eval(parse(text=paste0("o.m",value)))
      actuals[,i] <- y[(in.sample+i-h+1):(in.sample+i)]
      cat(paste(rep("\b",nchar(i)),collapse=""))
      cat(i)
    }
  }

  cat("\n")

  return(list(forecasts=forecasts,actuals=actuals))
}