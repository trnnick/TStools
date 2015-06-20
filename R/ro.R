ro <- function(data,h=1,origins=1,call,value=NULL,
               ci=FALSE,co=FALSE){
# Function makes Rolling Origin for the data using the call
#    Copyright (C) 2015  Yves Sagaert & Ivan Svetunkov
  l.value <- length(value)
  if(!is.null(value)){
    for(i in 1:l.value){
      if(substring(value[i],1,1)!="$"){
        value[i] <- paste0("$",value[i])
      }
    }
  }

# Check if the data is vector or ts
  if(!is.numeric(data) & !is.ts(data)){
    stop("The provided data is not a vector or ts object! Can't work with it!", call.=FALSE);
  }

# Write down the start and frequency
  data.start <- start(data)
  data.freq <- frequency(data)
  
# Write down the entire dataset y and the original horizon
  y <- data
  hh <- h
  obs <- length(y)
  in.sample <- obs - origins

# Make matrix to store the results origins in colnames,
#   different rows are entire horizons forecasted
#  forecasts <- array(NA,c(h,origins,length(value)))
#  dimnames(forecasts) <- list(paste0("h",c(1:h)),
#                              paste0("origin",c(1:origins)),
#                              value)
  actuals <- matrix(NA,nrow=h,ncol=origins)
  colnames(actuals) <- paste0("origin",c(1:origins))
  rownames(actuals) <- paste0("h",c(1:h))

#  forecasts <- actuals
  forecasts <- list(NA)
  cat(paste0("Origins done:  "))

  if(co==FALSE){
    for (i in 1:origins){
# Adjust forecasting horizon to not exeed the sample size
      h <- min(hh,obs - (in.sample+i-1))
# Make the in-sample
      if(ci==FALSE){
        data <- ts(y[1:(in.sample-h+i-1)],start=data.start,frequency=data.freq)
      }
      else{
        data <- ts(y[i:(in.sample-h+i-1)],start=data.start,frequency=data.freq)
      }
# Evaluate the call string and save to object o.m.
      o.m <- eval(parse(text=call))
# Save the forecast and the corresponding actuals in matrices
#      forecasts[1:h,i] <- eval(parse(text=paste0("o.m",value)))
      for(j in 1:l.value){
        forecasts[[(i-1)*l.value+j]] <- eval(parse(text=paste0("o.m",value[j])))
      }
      actuals[1:h,i] <- y[(in.sample+i):(in.sample+i+h-1)]
      cat(paste(rep("\b",nchar(i)),collapse=""))
      cat(i)
    }
  }
  else{
    for (i in 1:origins){
# Make the in-sample
      if(ci==FALSE){
        data <- ts(y[1:(in.sample-h+i-1)],start=data.start,frequency=data.freq)
      }
      else{
        data <- ts(y[i:(in.sample-h+i-1)],start=data.start,frequency=data.freq)
      }
# Evaluate the call string and save to object o.m.
      o.m <- eval(parse(text=call))
# Save the forecast and the corresponding actuals in matrices
#      forecasts[,i] <- eval(parse(text=paste0("o.m",value)))
      for(j in 1:l.value){
        forecasts[[(i-1)*l.value+j]] <- eval(parse(text=paste0("o.m",value[j])))
      }
      actuals[,i] <- y[(in.sample+i-h+1):(in.sample+i)]
      cat(paste(rep("\b",nchar(i)),collapse=""))
      cat(i)
    }
  }

  cat("\n")

  if(is.null(value)){
    value <- "$output"
  }
  value <- substring(value,2,nchar(value))

  returned.list <- list(actuals)

# Create several matrices from the list
  for(i in 1:l.value){
    stuff <- matrix(NA,nrow=length(forecasts[[i]]),ncol=origins)
    colnames(stuff) <- colnames(actuals)
    for(j in 1:origins){
      stuff[1:length(forecasts[[(j-1)*l.value+i]]),j] <- forecasts[[(j-1)*l.value+i]]
    }
    returned.list[[i+1]] <- stuff
  }

  names(returned.list) <- c("actuals",value)
  return(returned.list)
}