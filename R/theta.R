theta <- function(y,m=NULL,h=10,outplot=0,sign.level=0.05,
                  cost0=c("MSE","MdSE","MAE","MdAE"),
                  cost2=c("MSE","MdSE","MAE","MdAE"),
                  costs=c("MSE","MdSE","MAE","MdAE"),
                  multiplicative=c(TRUE,FALSE),cma=NULL){
# Theta method 
# This implementation of Theta method tests automatically for seasonality and trend.
# Seasonal decomposition can be done either additively or multiplicatively and the seasonality
# is treated as a pure seasonal model. The various components of Theta can be optimised 
# using different cost functions. The originally proposed Theta method always assume 
# multiplicative seasonality and presence of trend, while all theta lines are optimised 
# using MSE. Seasonality is estimated using classical decomposition.
#
# Inputs
#   y               Time series to model. Can be either a vector or a ts object
#   m               Periods in a season of the time series. If insample is a ts object then 
#                   this is taken from its frequency, unless overriden. 
#   h               Forecast horizon. Default is 10.
#   outplot         Provide plot:
#                     0: No plot
#                     1: Series and forecast
#                     2: As above with theta lines
#   sign.level      Significance level for trend and seasoanlity statistical tests.
#   cost0           Cost function of theta0 line.
#   cost2           Cost function of theta2 line.
#   costs           Cost function of seasonal element. 
#                   Costs may be: MSE, MdSE, MAE, MdAE.
#   multiplicative  If TRUE then multiplicative decomposition is performed. 
#                   Otherwise additive is used.
#   cma             Input pre-calculated centred moving average. 
#                   Use NULL to calculate internally.
#
# Output
#   frc         Forecasts.
#   exist       exist[1] is the result for trend, exist[2] is for season.
#   theta0      Forecasted values of theta0 line.
#   theta2      Forecasted values of theta2 line.
#   season      Forecasted values of seasonal element.
#   a           SES parameters of theta2.
#   b           Regression parameters of theta0.
#   g           Pure seasonal exponential smoothing parameters of season.
#
# Example:
#   theta(referrals,outplot=2)
#  
# Nikolaos Kourentzes, 2014 <nikolaos@kourentzes.com>

  # Defaults
  cost0 <- cost0[1]
  cost2 <- cost2[1]
  costs <- costs[1]
  multiplicative <- multiplicative[1]
  
  n <- length(y)
  
  # Get m (seasonality)
  if (is.null(m)){
    if (class(y) == "ts"){
      m <- frequency(y)
    } else {
      stop("Seasonality not defined (y not ts object).")
    }
  }
  
  # Check if CMA is given
  if (!is.null(cma)){
    if (n != length(cma)){
      stop("Length of series and cma do not match.")
    }
  } else {
    # Calculate CMA
    cma <- cmav(y,ma=m,outplot=0,fast=TRUE)
  }
  
  # Test for trend
  trend.exist <- coxstuart(cma)$p.value <= sign.level/2
  
  # Get seasonal matrix and test for seasonality
  if (m>1){
    if (multiplicative == TRUE){
      ynt <- y / cma
    } else {
      ynt <- y - cma
    }
    k <- m - (n %% m)
    ynt <- c(as.vector(ynt),rep(NA,times=k))
    ns <- length(ynt)/m
    ynt <- matrix(ynt,nrow=ns,ncol=m,byrow=TRUE)
    season.exist <- friedman.test(ynt)$p.value <= sign.level 
  } else {
    season.exist <- FALSE
  }
  
  # If seasonality exist then decompose 
  if (season.exist == TRUE){
    # y.des <- y/rep(ynt,ceiling(n/m)+1)[1:n]
    y.des <- cma
  } else {
    y.des <- y
  }
  
  # Create theta lines
  X <- cbind(matrix(1,nrow=n,ncol=1),matrix(c(1:n), nrow=n, ncol=1))
  if (trend.exist == TRUE){
    b0 <- solve(t(X)%*%X)%*%t(X)%*%y.des # Initialise theta0 parameters
    b <- opt.trnd(y.des,X,cost0,b0)      # Optimise theta0
  } else {
    # If no trend then theta0 is just a mean
    b <- rbind(mean(y.des),0)
  }
  theta0 <- X%*%b + 0*y.des           # 0*y.des To take ts object properties
  a0 <- rbind(0.1,y.des[1])           # Initialise theta2 parameters
  theta2 <- 2*y.des - theta0          # Construct theta2 
  a <- opt.ses(theta2,cost2,a0,2)     # Optimise theta2 
  
  # Prediction
  frc.theta0 <- b[1] + b[2]*((n+1):(n+h))
  frc.theta2 <- fun.ses(theta2,a)$out * rep(1,h)
  
  frc <- (frc.theta0 + frc.theta2)/2
  
  # Convert to ts object
  if (class(y) == "ts"){
    s <- end(y)
    if (s[2]==m){
      s[1] <- s[1]+1
      s[2] <- 1
    } else {
      s[2] <- s[2]+1
    }
    frc <- ts(frc,start=s,frequency=m)  
  } 
  
  # Reseasonalise
  if (season.exist == TRUE){
    # Seasonality is modelled with a pure seasonal smoothing
    sout <- opt.sfit(ynt,costs,n,m)
    season <- sout$season
    # sstd <- sd(season)
    season <- rep(season, h %/% m + 1)[1:h]
    g <- sout$g
    if (multiplicative == TRUE){
      frc <- frc * season
    } else {
      frc <- frc + season
    }
  } else {
    g <- NULL
    season <- NULL
    # sstd <- NULL
  }
  
  if (outplot==1){
    # Simple in-sample and forecast
    if (class(y) == "ts"){
      ts.plot(y,frc,gpars=list(col=c("black","blue"),lwd=c(1,2)))
    } else {
      ymin <- min(min(y),min(frc))
      ymax <- max(min(y),max(frc))
      yminmax <- c(ymin-0.1*(ymax-ymin),ymax+0.1*(ymax-ymin))
      plot(1:n,y,type="l",xlim=c(1,(n+h)),ylab="",xlab="Time",ylim=yminmax)
      lines((n+1):(n+h),frc,col="blue",type="l",lwd=2)
    }
  } 
  
  if (outplot==2){
    # As previous, including theta lines
    if (class(y) == "ts"){
      s <- end(y)
      if (s[2]==m){
        s[1] <- s[1]+1
        s[2] <- 1
      } else {
        s[2] <- s[2]+1
      }
      frc.theta0 <- ts(frc.theta0,start=s,frequency=m)  
      frc.theta2 <- ts(frc.theta2,start=s,frequency=m)  
      ts.plot(y,theta0,frc.theta0,theta2,frc.theta2,frc,
              gpars=list(col=c("black","forestgreen","forestgreen","red","red","blue"),
                         lwd=c(1,1,1,1,1,2),lty=c(1,1,2,1,2,1)))
    } else {
    ymin <- min(min(y),min(theta0),min(theta2),min(frc),min(frc.theta0),min(frc.theta2))
    ymax <- max(min(y),max(theta0),max(theta2),max(frc),max(frc.theta0),max(frc.theta2))
    yminmax <- c(ymin-0.1*(ymax-ymin),ymax+0.1*(ymax-ymin))
    plot(1:n,y,type="l",xlim=c(1,(n+h)),ylab="",xlab="Time",ylim=yminmax)
    lines(1:n,theta0,col="forestgreen",type="l")
    lines(1:n,theta2,col="red",type="l")
    lines((n+1):(n+h),frc.theta0,col="forestgreen",type="l",lty=2)
    lines((n+1):(n+h),frc.theta2,col="red",type="l",lty=2)
    lines((n+1):(n+h),frc,col="blue",type="l",lwd=2)
    }
  }
  
  # Prepare output
  exist <- rbind(trend.exist,season.exist)
  rownames(exist) <- c("Trend","Season")
  rownames(a) <- c("Alpha","Initial level")
  rownames(b) <- c("Beta","Intercept")
  if (season.exist==TRUE){
    g <- matrix(g,nrow=1,ncol=m+1)
    colnames(g) <- c("Gamma",paste("s",1:m,sep=""))
  }
  costf <- rbind(cost0,cost2,costs)
  rownames(costf) <- c("Theta0","Theta2","Seasonal")
  return(list("frc"=frc,"exist"=exist,"theta0"=frc.theta0,"theta2"=frc.theta2,
              "season"=season,"cost"=costf,"a"=a,"b"=b, "g"=g)) # ,"std.season"=sstd))
  
}

opt.sfit <- function(ynt,costs,n,m){
  # Optimise pure seasonal model and predict out-of-sample seasonality
  g0 <- c(0.001,colMeans(ynt,na.rm=TRUE))       # Initialise seasonal model
  season.sample <- matrix(t(ynt),ncol=1)        # Transform back to vector
  season.sample <- season.sample[!is.na(season.sample)]
  opt <- optim(par=g0, cost.sfit, method = "Nelder-Mead", season.sample=season.sample, 
               cost=costs, n=n, m=m, control = list(maxit = 2000))
  g <- opt$par
  season <- fun.sfit(season.sample,g,n,m)$outs
  return(list("season"=season,"g"=g))
}

fun.sfit <- function(season.sample,g,n,m){
  # Fit pure seasonal model
  s.init <- g[2:(m+1)]
  season.fit <- c(s.init,rep(NA,n))
  for (i in 1:n){
    season.fit[i+m] <- season.fit[i] + g[1]*(season.sample[i] - season.fit[i])
  }
  return(list("ins"=season.fit[1:n],"outs"=season.fit[(n+1):(n+m)]))  
}

cost.sfit <- function(g,season.sample,cost,n,m){
  # Cost function of pure seasonal model
  err <- season.sample-fun.sfit(season.sample,g,n,m)$ins
  err <- cost.err(err,cost,NULL)
  if (g[1]<0 | g[1]>1){
    err <- 9*10^99
  }
  return(err)   
}

fun.ses <- function(line,a){
  # Fit SES model on theta line
  n <- length(line)
  ses <- matrix(NA,nrow=n+1,ncol=1)
  ses[1] <- a[2]
  for (i in 2:(n+1)){
    ses[i] <- a[1]*line[i-1] + (1-a[1])*ses[i-1]
  }
  return(list("ins"=ses[1:n],"outs"=ses[n+1]))
}

opt.ses <- function(line,cost,a0,theta){
  # Optimise SES on theta
  opt <- optim(par=a0, cost.ses, method = "Nelder-Mead", line=line, cost=cost, 
               theta=theta, control = list(maxit = 2000))
  a <- opt$par
  return(a)
}

cost.ses <- function(a,line,cost,theta=2){
  # Cost function for SES optimisation
  err <- line-fun.ses(line,a)$ins
  err <- cost.err(err,cost,theta)
  if (!a[1]<0.99 | !a[1]>0.01){
    err <- 9*10^99
  }
  return(err)
}

opt.trnd <- function(y,X,cost,b0){
  # Optimise theta line 0
  opt <- optim(par=b0, cost.trnd, method = "Nelder-Mead", y=y, X=X, 
               cost=cost, control = list(maxit = 2000))
  b <- opt$par
  return(b)
}

cost.trnd <- function(b,y,X,cost){
  # Theta 0 cost function
  err <- y-X%*%b
  err <- cost.err(err,cost,NULL)
  return(err)
}

cost.err <- function(err,cost,theta=NULL){
  # Cost calculation
  if (cost == "MAE"){
    err <- mean(abs(err))
  }
  if (cost == "MdAE"){
    err <- median(abs(err))
  }
  if (cost == "MSE"){
    err <- mean((err)^2)
  }
  if (cost == "MdSE"){
    err <- median((err)^2)
  }
  if (cost == "MTE"){
    err <- mean(abs((err)^theta))
  }
  if (cost == "MdTE"){
    err <- median(abs((err)^theta))
  }
  return(err)
}