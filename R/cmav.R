cmav <- function(y,ma=NULL,fill=c(TRUE,FALSE),outplot=c(FALSE,TRUE),fast=c(TRUE,FALSE)){
# Calculate centred moving average
#
# Inputs:
#   y             Time series vector (can be ts object)
#   ma            Length of centred moving average. If y is a ts object 
#                 then the default is its frequency
#   fill          If TRUE then fill first and last ma/2 observations using ETS
#   outplot       If TRUE then produce plot
#   fast          If TRUE then only a limited set of models are evaluated for CMA extrapolation
#
# Outputs:
#   cma           Centred moving average. If y is a ts object, then cma has the
#                 same properties
#
# Example:
#   cmav(wineind,outplot=TRUE)
#
# Nikolaos Kourentzes, 2014 <nikolaos@kourentzes.com>
    
  # Set defaults
  fill <- fill[1]
  outplot <- outplot[1]
  fast <- fast[1]
  
  # Get MA length
  if (is.null(ma)){
    if (any(class(y) == "ts")){
      ma <- frequency(y)
    } else {
      stop("MA length not defined (y not ts object).")
    }
  }
  
  # Get bounds for MA and correct length
  n <- length(y)
  mlbounds <- c(floor(ma/2)+1, n-floor(ma/2))
  isodd = ma %% 2 
  if (isodd == 0){
    ml <- ma+1
  } else {
    ml <- ma
  }
  
  # Calculate MA
  # Loop across MA-order to speed up things
  mamat <- matrix(NA,nrow=ml,ncol=(n-ml+1))
  for (i in 1:ml){
    mamat[i,] <- y[i:(n-ml+i)]
  }
  if (isodd == 0){
    mamat[c(1,ml),] <- mamat[c(1,ml),]/2
  }
  mamat <- colSums(mamat)/ma
  
  cma <- y
  cma[] <- NA
  cma[mlbounds[1]:mlbounds[2]] <- mamat
  
  # Fill MA is requested
  if (fill == TRUE){
    if (fast == FALSE){
      if ((n-mlbounds[2]) >= 1){   
        cma[(mlbounds[2]+1):n] <- as.vector(forecast(ets(cma[(mlbounds[1]:mlbounds[2])], 
                                                       model="ZZN"),h=(n-mlbounds[2]))$mean)
      }
      if ((mlbounds[1]-1) >= 1){   
      cma[1:(mlbounds[1]-1)] <- rev(as.vector(forecast(ets(rev(cma[(mlbounds[1]:mlbounds[2])]),
                                                           model="ZZN"),h=(mlbounds[1]-1))$mean))
      }
    } else {
      fit <- ets(cma[(mlbounds[1]:mlbounds[2])],model="AZN")
      if ((n-mlbounds[2]) >= 1){   
        cma[(mlbounds[2]+1):n] <- as.vector(forecast(fit,h=(n-mlbounds[2]))$mean)
      }
      if ((mlbounds[1]-1) >= 1){   
        cma[1:(mlbounds[1]-1)] <- rev(as.vector(forecast(ets(rev(cma[mlbounds[1]:mlbounds[2]]),
                                                             fit,use.initial.values=FALSE),h=(mlbounds[1]-1))$mean))
      }
    }
  }
  
  if (outplot == TRUE){
    plot(y)
    lines(cma,col="red")
  }
  
  return(cma)
  
}