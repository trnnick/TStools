library(forecast)

init.ets <- function(y,level=NULL,trend=NULL,season=NULL,...){
# This function overrides the ETS initial values. The model is not reoptimised and therefore 
# should be used in conjuction with pre-defined smoothing parameters
#
# Inputs:
#   y           A numeric vector or time series.
#   level       Initial level.
#   trend       Initial trend.
#   season      Vector of initial seasonal indices.
#   ...         ets input arguments.
#
# Output:
#   fit         An object of class "ets".
#
# Example:
#   init.ets(AirPassengers,trend=1.1)
#
# Nikolaos Kourentzes, 2014 <nikolaos@kourentzes.com>
  
  # Fit model and get initial states
  fit <- ets(y,...)
  initstate <- fit$initstate
  ninit <- names(initstate)
  
  rerun <- FALSE
  
  # Override states
  # Level
  if (!is.null(level)){
    if (length(level) != 1){
      stop("Use a single value for the initial level.")
    } 
    initstate[ninit == "l"] <- level
    rerun <- TRUE
  }
  # Trend
  if (!is.null(trend) && sum(ninit == "b")>0){
    if (length(trend) != 1){
      stop("Use a single value for the initial trend.")
    }
    initstate[ninit == "b"] <- trend
    rerun <- TRUE
  }
  # Season
  if (!is.null(season) && sum(ninit == "s1")>0){
    s.idx <- grep("s",ninit)
    if (length(season) != length(s.idx)){
      stop("Length of initial values for season do not match with ETS. Check time series seasonality.")
    } 
    initstate[s.idx] <- season
    rerun <- TRUE
  }
  
  if (rerun==TRUE){
    # Rerun fit without re-optimisation
    fit$initstate <- initstate
    fit <- ets(y,fit,use.initial.values=TRUE,...)
  }
      
  # Return result
  return(fit)

}