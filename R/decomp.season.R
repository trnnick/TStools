decomp.season <- function(y,m=NULL,s=NULL,trend=NULL,outplot=c(FALSE,TRUE),
                   decomposition=c("multiplicative","additive"),
                   h=0,type=c("mean","median","pure.seasonal"))
{
  # Seasonal plots and crude trend/season tests
  #
  # Inputs:
  #   y             Time series vector (can be ts object)
  #   m             Seasonal period. If y is a ts object then the default is its frequency  
  #   s             Starting period in the season. If y is a ts object then default is read
  #   trend         If TRUE then a trend is assumed and is removed using CMA
  #                 If FALSE then no trend is assumed
  #                 If NULL then trend is identified and removed if found
  #   colour        Single colour override for plots
  #   alpha         Significance level for statistical tests (kpss and friedman)
  #   outplot       Provide plot output:
  #                   0 - None
  #                   1 - Seasonal diagramme
  #                   2 - Seasonal boxplots
  #                   3 - Seasonal subseries
  #                   4 - Seasonal distribution
  #                   5 - Seasonal density
  #   decomposition typee of seasonal decompositio: "multiplicative" or "additive".
  #   cma           Input precalculated level/trend for the analysis. Overrides trend=NULL. 
  #
  # Outputs:
  #   List with the following elements:
  #     season        Matrix of (detrended) seasonal elements
  #     season.exist  TRUE/FALSE results of friedman test
  #     season.pval   Friedman test p-value
  #     trend         CMA estimate or NULL if trend == FALSE
  #     trend.exist   TRUE/FALSE result of kpss test
  #     trend.pval    kpss (approximate) p-value
  #
  # Nikolaos Kourentzes, 2014 <nikolaos@kourentzes.com>
  
  # Defaults
  outplot <- outplot[1]
  decomposition <- decomposition[1]
  type <- type[1]
  
  # Get m (seasonality)
  if (is.null(m)){
    if (class(y) == "ts"){
      m <- frequency(y)
    } else {
      stop("Seasonality not defined (y not ts object).")
    }
  }
  
  # Get starting period in seasonality if available
  if (is.null(s)){
    if (class(y) == "ts"){
      s <- start(y)
      s <- s[2]
    } else {
      s <- 1
    }
  } 
  
  n <- length(y)
  
  if ((decomposition == "multiplicative") && (min(y)<=0)){
    decomposition <- "additive"
  }
  
  # If trend is not given then calculate CMA
  if (is.null(trend)){
    trend <- cmav(y=y,ma=m,fill=TRUE,outplot=FALSE,fast=TRUE)
  } else {
    if (n != length(trend)){
      stop("Length of series and trend input do not match.")
    }
  }
  
  if (decomposition == "multiplicative"){
    ynt <- y/trend  
  } else {
    ynt <- y - trend
  }
    
  ymin <- min(ynt)
  ymax <- max(ynt)
  ymid <- median(ynt)
  
  # Fill with NA start and end of season
  k <- m - (n %% m)
  ks <- s-1
  ke <- k-ks
  ynt <- c(rep(NA,times=ks),as.vector(ynt),rep(NA,times=ke))
  ns <- length(ynt)/m
  ynt <- matrix(ynt,nrow=ns,ncol=m,byrow=TRUE)
  colnames(ynt) <- paste("p",1:m,sep="")
  rownames(ynt) <- paste("s",1:ns,sep="")
  
  # If h>0 then produce forecasts of seasonality 
  g <- NULL
  if (h == 0){
    f.season <- NULL
  } else {
    if (type=="mean"){
      # Calculate the seasonality as the overall mean
      season <- colMeans(ynt, na.rm=TRUE)
      f.season <- rep(season,h+m*2)
      f.season <- as.vector(f.season[(m-ke+1):(m-ke+h-1)])
      i.season <- rep(season,ns)
      i.season <- i.season[(ks+1):(ns*m-ke)]+y*0
    }
    if (type=="median"){
      # Calculate the seasonality as the overall median
      season <- array(NA,c(1,m))
      for (si in 1:m){
        season[si] <- median(ynt[,si], na.rm=TRUE)
      }
      f.season <- rep(season,h+m*2)
      f.season <- as.vector(f.season[(m-ke+1):(m-ke+h)])
      i.season <- rep(season,ns)
      i.season <- i.season[(ks+1):(ns*m-ke)] + y*0
    }
    if (type="pure.seasonal"){
      # Seasonality is modelled with a pure seasonal smoothing
      sout <- opt.sfit(ynt,"MSE",0.001,n,m)
      g <- sout$g
      season <- sout$season
      f.season <- rep(season,h+m*2)
      f.season <- as.vector(f.season[(m-k+1):(m-k+h)])
      f.season <- rep(season, h %/% m + 1)[1:h]
      season.sample <- matrix(t(ynt),ncol=1)        
      season.sample <- season.sample[!is.na(season.sample)]
      i.season <- as.vector(fun.sfit(season.sample,g,n,m)$ins) + y*0
    }
    
    # Convert f.season to ts object
    if (class(y) == "ts"){
      s <- end(y)
      if (s[2]==m){
        s[1] <- s[1]+1
        s[2] <- 1
      } else {
        s[2] <- s[2]+1
      }
      f.season <- ts(f.season,start=s,frequency=m)  
    } 
  } 
  
  if (decomposition == "multiplicative"){
    resid <- y - (trend*i.season)
  } else {
    resid <- y - (trend+i.season)
  }
  
  # Produce plots
  if (outplot == TRUE){
    ts.plot(i.season,f.season, gpars=list(col=c("black","blue"),
            lwd=c(1,2),lty=c(1,1)))
  }
  
  return(list(trend=trend,season=i.season,irregular=resid,
              f.season=f.season,g=g))
  
}


opt.sfit <- function(ynt,costs,g0,n,m){
  # Optimise pure seasonal model and predict out-of-sample seasonality
  g0 <- c(g0,colMeans(ynt,na.rm=TRUE))          # Initialise seasonal model
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
  err <- cost.err(err,cost)
  if (g[1]<0 | g[1]>1){
    err <- 9*10^99
  }
  return(err)   
}

cost.err <- function(err,cost){
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
  return(err)
}
