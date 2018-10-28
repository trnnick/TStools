theta <- function(y,m=NULL,sign.level=0.05,
                  cost0=c("MSE","MdSE","MAE","MdAE"),
                  cost2=c("MSE","MdSE","MAE","MdAE"),
                  costs=c("MSE","MdSE","MAE","MdAE"),
                  multiplicative=c(TRUE,FALSE),cma=NULL,
                  outliers=NULL){
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
#   sign.level      Significance level for trend and seasoanlity statistical tests.
#   cost0           Cost function of theta0 line.
#   cost2           Cost function of theta2 line.
#   costs           Cost function of seasonal element. 
#                   Costs may be: MSE, MdSE, MAE, MdAE.
#   multiplicative  If TRUE then multiplicative decomposition is performed. 
#                   Otherwise additive is used.
#   cma             Input pre-calculated centred moving average. 
#                   Use NULL to calculate internally.
#   outliers        Optional. Provide vector of location of observations that are considered outliers.
#                   These will be included in theta0 estimation. To consider no outliers then use NULL.
#
# Output
#   method          Forecasting method.
#   y               Input time series.
#   m               Periods in a season of the time series.
#   exist           exist[1] is the result for trend, exist[2] is for season.
#   multiplicative  If TRUE seasonality is modelled multiplicatively.
#   theta0          Fitted theta0 line values.
#   theta2          Fitted theta2 line values.
#   season          Seasonal profile.
#   a               SES parameters of theta2.
#   b               Regression parameters of theta0.
#   p               Coefficients of outliers from theta0 and theta2 estimation.
#   g               Pure seasonal exponential smoothing parameters of season.
#   fitted          Fitted values.
#   redisuals       In-sample residuals.
#   MSE             In-sample Mean Squared Error.
#
# Example:
#   theta(referrals)
#  
# Nikolaos Kourentzes, 2014 <nikolaos@kourentzes.com>
# Updates 2017.2: Use S3method

  # Defaults
  cost0 <- match.arg(cost0,c("MSE","MdSE","MAE","MdAE"))
  cost2 <- match.arg(cost2,c("MSE","MdSE","MAE","MdAE"))
  costs <- match.arg(costs,c("MSE","MdSE","MAE","MdAE"))
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
  if (m>1 && (length(y)/m)>=2){
    k <- m - (n %% m)
    if (multiplicative == TRUE){
      ynt <- y / cma
    } else {
      ynt <- y - cma
    }
    k <- m - (n %% m)
    if (k == m){k <- 0}
    ynt <- c(as.vector(ynt),rep(NA,times=k))
    ns <- length(ynt)/m
    ynt <- matrix(ynt,nrow=ns,ncol=m,byrow=TRUE)
    if (diff(range(colSums(ynt))) <= .Machine$double.eps ^ 0.5){
        season.exist <- FALSE
    } else {
        season.exist <- friedman.test(ynt)$p.value <= sign.level 
    }
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
  
  # Include in theta0 outliers if any provided
  if (is.null(outliers)){
    X.out <- NULL
    n.out <- 0
  } else {
    n.out <- length(outliers)
#     X.out <- array(0,c(n,n.out))
#     X.out[outliers+n*(0:(n.out-1))] <- 1
    # Through the MA the outlier is spread across m observations, create dummies to account for that
    m.half <- floor((m + (m+1) %% 2)/2)
    X.out <- array(0,c(n,n.out))
    for (i in 1:n.out){
      X.out[max(1,outliers[i]-m.half):min(n,outliers[i]+m.half),i] <- 1
    }
  }
  # Include trend component in theta0
  if (trend.exist == TRUE){
    X.trend <- matrix(c(1:n), nrow=n, ncol=1)
  } else {
    X.trend <- NULL
  }
  X <- cbind(matrix(1,nrow=n,ncol=1),X.trend,X.out)
  
#   if (trend.exist == TRUE || !is.null(outliers)){
  if ((n.out + trend.exist)>=1){
    b0 <- solve(t(X)%*%X)%*%t(X)%*%y.des # Initialise theta0 parameters
    b <- opt.trnd(y.des,X,cost0,b0)      # Optimise theta0
  } else {
    # If no trend then theta0 is just a mean
    b <- mean(y.des)
  }

  # Create theta0 line without outliers
  # 0*y.des To take ts object properties
  theta0 <- X[,1:(1+trend.exist)]%*%matrix(b[1:(1+trend.exist)],ncol=1) + 0*y.des               
  
  # Estimate SES
  theta2 <- 2*y.des - theta0              # Construct theta2 
  a0 <- rbind(0.1,theta2[1],rep(0,n.out)) # Initialise theta2 parameters
  a <- opt.ses(theta2,cost2,a0,2,X.out)   # Optimise theta2 
  
  # In-sample fit
  in.theta0 <- theta0
  in.theta2 <- fun.ses(theta2,a,X.out)$ins
  if (!is.null(X.out)){
    # Remove outlier from fit - the complete outlier will be modelled afterwards
    in.theta2 <- in.theta2 - X.out %*% matrix(a[3:(2+n.out)])
  }
  in.fit <- (in.theta0 + in.theta2)/2
  
  # Separate theta0 and theta2 parameters into a, b and outliers (p[,1:2])
  if (n.out > 0){
    p <- matrix(b[(2+trend.exist):(1+trend.exist+n.out)],nrow=n.out)
    p <- cbind(p,a[3:(2+n.out)])
    colnames(p) <- c('Theta0','Theta2')
  } else {
    p <- NULL
  }
  if (trend.exist == FALSE){
    b <- rbind(b[1], 0)
  } else {
    b <- b[1:2]
  } 
  b <- matrix(b,ncol=1)
  a <- matrix(a[1:2],ncol=1)

  # Reseasonalise
  if (season.exist == TRUE){
    # Seasonality is modelled with a pure seasonal smoothing
    sout <- opt.sfit(ynt,costs,n,m,y,in.fit,multiplicative,outliers)
    season <- sout$season
    g <- sout$g
    if (n.out > 0){
      p <- cbind(p,matrix(sout$p,ncol=1,dimnames=list(NULL,'Season')))
    }
    # sout$in.season includes the outlier
    if (multiplicative == TRUE){
      in.fit <- in.fit * sout$in.season
    } else {
      in.fit <- in.fit + sout$in.season
    }
  } else {
    g <- NULL
    season <- NULL
  }
  
  # Prepare output
  exist <- rbind(trend.exist,season.exist)
  rownames(exist) <- c("Trend","Season")
  rownames(a) <- c("Alpha","Initial level")
  rownames(b) <- c("Intercept","Slope")
  if (season.exist==TRUE){
    g <- matrix(g,nrow=1,ncol=m+1)
    colnames(g) <- c("Gamma",paste("s",1:m,sep=""))
  }
  costf <- rbind(cost0,cost2,costs)
  rownames(costf) <- c("Theta0","Theta2","Seasonal")
  if (class(y) == "ts"){
    in.fit <- ts(in.fit,frequency=m,end=end(y))    
    theta0 <- ts(theta0,frequency=m,end=end(y))    
    theta2 <- ts(theta2,frequency=m,end=end(y))    
  }
  
  return(structure(list("method"="Theta","y"=y,"m"=m,"exist"=exist,
                        "multiplicative"=multiplicative,
                        "theta0"=theta0,"theta2"=theta2,
                        "season"=season,"x.out"=X.out,
                        "cost"=costf,"a"=a,"b"=b,"p"=p,"g"=g,
                        "fitted"=in.fit,"residuals"=y-in.fit,
                        "MSE"=mean((y-in.fit)^2))
                   ,class="theta")) 
  
}

forecast.theta <- function(object,h=NULL,...){
    # Produce forecasts with Theta
    
    a <- object$a
    b <- object$b
    p <- object$p
    m <- object$m
    n <- length(object$y)
    X.out <- object$x.out
    theta2 <- object$theta2
    y <- object$y
    
    if (is.null(h)){
        h <- m
    }
    
    # Forecast theta line
    frc.theta0 <- b[1] + b[2]*((n+1):(n+h))
    if (!is.null(X.out)){
        #### This is added because this is not defined anywhere in this function! ####
        n.out <- length(X.out)
        ####
        a.frc <- rbind(a,array(p[,2],c(n.out,1)))
    } else {
        a.frc <- a
    }
    frc.theta2 <- fun.ses(theta2,a.frc,X.out)$outs * rep(1,h)
    frc <- (frc.theta0 + frc.theta2)/2
    
    # Include seasonality
    if (object$exist[2] == TRUE){
        season <- rep(object$season, h %/% m + 1)[1:h]
        if (object$multiplicative == TRUE){
            frc <- frc * season
        } else {
            frc <- frc + season
        }
    } else {
        season <- NULL
    }
    
    # Convert to ts
    if (class(y) == "ts"){
        fstart <- c(end(y)[1],end(y)[2]+1)
        if (is.na(fstart[2])){  # If the second element of end(y) does not exist because it is fractional
            fstart <- end(y) + deltat(y)
        }
        frc <- ts(frc,start=fstart,frequency=m)  
        frc.theta0 <- ts(frc.theta0,start=fstart,frequency=m)  
        frc.theta2 <- ts(frc.theta2,start=fstart,frequency=m)  
    } 
    
    return(structure(list("method"=class(object),"mean"=frc,
                "frc.theta0"=frc.theta0,"frc.theta2"=frc.theta2,
                "frc.season"=season,"x"=y,
                "fitted"=object$fitted,"residuals"=object$residuals),
                class="forecast"))
    
}

plot.theta <- function(x,thetalines=c(TRUE,FALSE),...){
    # Produce in-sample fit plot
    
    thetalines = thetalines[1]
    is.ts <- class(x$y) == "ts"
    
    # Default limits of plot
    yy <- range(c(x$y,x$fitted,x$theta0,x$theta2))
    yy <- yy+c(-1,1)*0.04*diff(yy)
    if (is.ts){
        xx <- time(x$y)[c(1,length(x$y))]
    } else {
        xx <- 1:length(x$y)
    }
    
    # Allow user to override plot defaults
    args <- list(...)
    if (!("main" %in% names(args))){
        args$main <- "Theta method"
    }
    if (!("xlab" %in% names(args))){
        args$xlab <- "Time"
    }
    if (!("ylab" %in% names(args))){
        args$ylab <- ""
    }
    if (!("xlim" %in% names(args))){
        args$xlim <- xx
    }
    if (!("ylim" %in% names(args))){
        args$ylim <- yy
    }
    # Remaining defaults
    args$x <- args$y <- NA
    # Use do.call to use manipulated ellipsis (...)
    do.call(plot,args)
    # Plot the rest
    lines(x$y,col="black")
    lines(x$fitted,col="blue")
    if (thetalines == TRUE){
        lines(x$theta0, col="red")
        lines(x$theta2, col="forestgreen")
    }
    
}

summary.theta <- function(object,...){
    print(object)
}

print.theta <- function(x,...){
    
    if (all(x$exist == FALSE)){
        mdl <- ""
    } else {
        mdl <- "("
        if (x$exist[1]){
            mdl <- paste0(mdl,"trend")
        }
        if (all(x$exist)){
            mdl <- paste0(mdl,",")
        }
        if (x$exist[2]){
            if (x$multiplicative){
                stp <- "multiplicative "
            } else {
                stp <- "additive "
            }
            mdl <- paste0(mdl,stp,"season (",x$m,")")
        }
        mdl <- paste0(mdl,")")
    }
    
    writeLines(paste("Theta method",mdl))
    writeLines("")
    writeLines("Parameters:")
    writeLines(paste("    Theta 0, intercept =",round(x$b[1],3)))
    if (x$exist[1]){
        writeLines(paste("    Theta 0, slope     =",round(x$b[2],3)))
    }
    writeLines(paste("    Theta 2, alpha     =",round(x$a[1],3)))
    if (x$exist[2]){
        writeLines(paste("    Season, gamma      =",round(x$g[1],3)))
    }
    writeLines("")
    writeLines("Initial states:")
    writeLines(paste("    Theta 2            =",round(x$a[2],3)))
    if (x$exist[2]){
        writeLines(paste("    Season             =",
                         paste(round(x$g[2:(x$m+1)],3),collapse=", ")))
    }
    writeLines("")
    writeLines(paste0("RMSE: ",round(sqrt(x$MSE),3)))

}

theta.thief <- function(y,h=NULL,...){
    # This is a wrapper function to use Theta with THieF
    
    # Remove level input from ellipsis
    ellipsis.args <- list(...)
    ellipsis.args$level <- NULL
    ellipsis.args$y <- y
    
    # Fit network
    fit <- do.call(theta,ellipsis.args) 
    
    # Default h
    if (is.null(h)){
        h <- frequency(y)
    }
    
    # Forecast
    out <- forecast(fit,h)

    return(out)
    
}

opt.sfit <- function(ynt,costs,n,m,y,in.fit,multiplicative,outliers){
  # Optimise pure seasonal model and predict out-of-sample seasonality
  if (is.null(outliers)){
    g0 <- c(0.001,colMeans(ynt,na.rm=TRUE))       # Initialise seasonal model
    season.sample <- matrix(t(ynt),ncol=1)        # Transform back to vector
    season.sample <- season.sample[!is.na(season.sample)]
    X.out <- NULL
    n.out <- 0
  } else {
    n.out <- length(outliers)
    g0 <- c(0.001,colMeans(ynt,na.rm=TRUE),rep(0,n.out))       # Initialise seasonal model
    if (multiplicative == TRUE){
      season.sample <- as.numeric(y / in.fit)
    } else {
      season.sample <- as.numeric(y - in.fit)
    }
    X.out <- array(0,c(n,n.out))
    X.out[outliers+n*(0:(n.out-1))] <- 1
  }
  opt <- optim(par=g0, cost.sfit, method = "Nelder-Mead", season.sample=season.sample, 
               cost=costs, n=n, m=m, X.out = X.out, control = list(maxit = 2000))
  g <- opt$par
  sfit <- fun.sfit(season.sample,g,n,m,X.out)
  out.season <- sfit$outs
  in.season <- sfit$ins
  # Size of outliers
  if (n.out > 0){
    p <- g[(m+2):(m+1+n.out)]
  } else {
    p <- NULL
  }
  return(list("season"=out.season,"in.season"=in.season,"g"=g[1:(1+m)],"p"=p))
}

fun.sfit <- function(season.sample,g,n,m,X.out){
  # Fit pure seasonal model
  s.init <- g[2:(m+1)]
  season.fit <- c(s.init,rep(NA,n))
  for (i in 1:n){
    season.fit[i+m] <- season.fit[i] + g[1]*(season.sample[i] - season.fit[i])
  }
  if (!is.null(X.out)){
    n.out <- length(X.out[1,])
    g.out <- matrix(g[(m+2):(m+1+n.out)],ncol=1)
    X.out <- rbind(X.out, array(0,c(m,n.out)))
    season.fit <- season.fit + X.out %*% g.out
  }
  return(list("ins"=season.fit[1:n],"outs"=season.fit[(n+1):(n+m)]))  
}

cost.sfit <- function(g,season.sample,cost,n,m,X.out){
  # Cost function of pure seasonal model
  err <- season.sample-fun.sfit(season.sample,g,n,m,X.out)$ins
  err <- cost.err(err,cost,NULL)
  if (g[1]<0 | g[1]>1){
    err <- 9*10^99
  }
  return(err)   
}

fun.ses <- function(line,a,X.out=NULL){
  # Fit SES model on theta line
  n <- length(line)
  ses <- matrix(NA,nrow=n+1,ncol=1)
  ses[1] <- a[2]
  for (i in 2:(n+1)){
    ses[i] <- a[1]*line[i-1] + (1-a[1])*ses[i-1]
  }
  # If X.out is !null then model outliers
  if (!is.null(X.out)){
    n.out <- length(X.out[1,])
    a.out <- matrix(a[3:(n.out+2)],ncol=1)
    X.out <- rbind(X.out, array(0,c(1,n.out)))
    ses <- ses + X.out %*% a.out
  }
  return(list("ins"=ses[1:n],"outs"=ses[n+1]))
}

opt.ses <- function(line,cost,a0,theta,X.out=NULL){
  # Optimise SES on theta
  opt <- optim(par=a0, cost.ses, method = "Nelder-Mead", line=line, cost=cost, 
               theta=theta, X.out=X.out, control = list(maxit = 2000))
  a <- opt$par
  return(a)
}

cost.ses <- function(a,line,cost,theta=2,X.out=NULL){
  # Cost function for SES optimisation
  err <- line-fun.ses(line,a,X.out)$ins
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
