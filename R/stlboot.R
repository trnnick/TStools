stlboot  <- function(ts,k=1,test.season=c(TRUE,FALSE),outplot=c(FALSE,TRUE)){
# STL/Loess bootstrapping by Bergmeir, Hyndman and Benitez, 2014
# http://robjhyndman.com/working-papers/bagging-ets/
#
# Inputs:
#   ts          Time series to bootstrap. Must be ts object.
#   k           Number of bootstraps.
#   test.season If TRUE then test for presence of seasonality. If FALSE then all non-yearly 
#               are decomposed using STL, assuming seasonality.
#   outplot     If TRUE provide a plot of the bootstrapped series.
#
# Output:
#   ts.recon    Array of bootstrapped series. Each row is a time series.
#
# Example:
#   stlboot(referrals,k=20,outplot=TRUE)
#
# Nikolaos Kourentzes, 2015 <nikolaos@kourentzes.com>
  
  test.season <- test.season[1]
  outplot <- outplot[1]
      
  # Box-Cox transformation
  lambda <- BoxCox.lambda(ts,method="guerrero",lower=0,upper=1)
  ts.bc<- BoxCox(ts,lambda)
  
  # Check if series is seasonal
  m <- frequency(ts)
  if (m > 1){
    season.exist <- TRUE
    if (test.season == TRUE){
      season.exist <- seasplot(ts.bc,outplot=0)$season.exist      
    }
  } else {
    season.exist <- FALSE
  }
  
  # Decompose time series
  if (season.exist==TRUE){
    # Seasonal time series, apply STL
    ts.decomp <- stl(ts.bc, s.window = "periodic", s.degree = 0, t.degree = 1, l.degree=1)
    remainder <- ts.decomp$time.series[,"remainder"]
  } else {
    # Nonseasonal series, use Loess
    x <- 1:length(ts.bc)
    ts.decomp <- loess(ts.bc~ x, data.frame(x=x, y=ts.bc),degree=1,span=6/length(ts.bc))
    ts.decomp.trend <- predict(ts.decomp, data.frame(x=x))
    remainder <- ts.bc - ts.decomp.trend
  }    
  
  # MBB bootstrapping
  n <- length(remainder) 
  l <- m*2        # Find block size
  boot <- function(x){sapply(x,function(x){remainder[(x-l+1):x]})}  # Function that bootstraps blocks
  boot.sample <- array(NA, c(k,n)) 
  for (i in 1:k) {           
    endpoint <- sample(l:n,size=(floor(n/l)+2))                             # Block ends
    ts.bt <- as.vector(sapply(endpoint,function(x){remainder[(x-l+1):x]}))  # Bootstrap blocks
    ts.bt <- ts.bt[-(1:(sample(l,size=1)-1))]                               # Trim random number of starting obs
    ts.bt <- ts.bt[1:n]                                                     # Trim to correct size
    boot.sample[i,] <- ts.bt      
  }
  
  # Reconstruct time series
  if (season.exist==TRUE){
    ts.bc.recon <- boot.sample + 
      t(replicate(k, as.vector(ts.decomp$time.series[,"trend"]))) + 
      t(replicate(k, as.vector(ts.decomp$time.series[,"seasonal"])))
    ts.recon <- InvBoxCox(ts.bc.recon,lambda)
  } else {
    ts.bc.recon <- boot.sample + t(replicate(k, ts.decomp.trend))
    ts.recon <- InvBoxCox(ts.bc.recon,lambda)
  }
  rownames(ts.recon) <- paste0('Boot',1:k)
  
  # Produce plot
  if (outplot == TRUE){
    plot(1:n,ts,type="l",xlab="Period",ylab="",lwd=2)
    for (i in 1:k){
      lines(1:n,ts.recon[i,],col="red")
    }
    lines(1:n,ts,col="black",lwd=2)
  }
  
  # Convert output to ts object
  ts.recon <- t(ts.recon)
  ts.recon <- ts(ts.recon, frequency=m, start=start(ts))
  
  return(ts.recon)
  
}