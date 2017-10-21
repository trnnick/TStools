coxstuart <- function(y,type=c("trend","deviation","dispersion"),alpha=0.05){
# Perform Cox - Stuart test for location and dispersion  
#
# H0 --> There is no trend present in location/dispersion
# H1 --> There is trend (upwards or downwards)  
#
# Inputs
#   y         Vector of data.
#   type      Type of test:
#               trend - test for changes in trend [default];
#               deviation - test for changes in deviation;
#               dispersion - test for changes in dispersion (range).
#   alpha     Significance level.
#
# Outputs
#   H         Hypothesis (H0/H1).
#   p.value   P-value.
#   Htxt      Description of the result.
#
# Example
#   coxstuart(referrals)
#
# Nikolaos Kourentzes, 2014 <nikolaos@kourentzes.com>
  
  type <- match.arg(type,c("trend","deviation","dispersion"))
  
  switch(type,
         "trend" = {
           data <- y
         },
         "deviation"={
           data <- splitdata(y)
           n <- dim(data)[2]
           v.sd <- Vectorize(function(i){temp <- sd(data[,i])})
           data <- v.sd(1:n)
         },
         "dispersion"={
           data <- splitdata(y)
           n <- dim(data)[2]
           v.range <- Vectorize(function(i){temp <- range(data[,i])
                                            temp <- max(temp)-min(temp)})
           data <- v.range(1:n)
         })
  
  # Find number of pairs for comparison
  n <- length(data)
  C <- ceiling(n/2)
  
  # Create pairs
  idx1 <- 1:(n-C)
  idx2 <- (1+C):n
  pair <- data[idx1] - data[idx2]
  
  # Calculate statistic
  Nplus <- sum(pair>0)
  Nminus <- sum(pair<0)
  stat <- min(c(Nplus,Nminus))
  
  # P-value
  if (sum(pair)!=0){
    p <- pbinom(stat,C,0.5)
  } else {
    p <- 1
  }
  
  if (p <= alpha/2){
    H <- 1
    txt <- "H1: There is trend (upwards or downwards)"
  } else {
    H <- 0
    txt <- "H0: There is no trend present in location/dispersion"
  }
  
  return(list(H=H,p.value=p,Htxt=txt))
  
}
 
# ------------------------------------------
splitdata <- function(y){
# Helper function
  
  n <- length(y)
  
  # Find K according to Cox Stuart guidelines
  if (n < 48){
    k <- 2
  } else if(n < 64){
    k <- 3
  } else if(n < 90){
    k <- 4
  } else {
    k <- 5
  }
  
  # Split time series to subsamples
  srem <- n %% k
  ktimes <- floor(n/k)
  idx <- array(1:(n-srem),c(k,ktimes))
  idx[,(round(ktimes/2)+1):ktimes] <- idx[,(round(ktimes/2)+1):ktimes] + srem   
  data <- array(y[idx], c(k,ktimes))
    
}
