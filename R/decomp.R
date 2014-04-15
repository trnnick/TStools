decomp <- function(y,m=NULL,s=NULL,trend=NULL,outplot=c(FALSE,TRUE),
                   decomposition=c("multiplicative","additive"),
                   h=0,type=c("mean","median","pure.seasonal"))
{
# Decomposition of series
#
# Inputs:
#   y             Time series vector (can be ts object).
#   m             Seasonal period. If y is a ts object then the default is its frequency.
#   s             Starting period in the season. If y is a ts object then default is read.
#   trend         Vector of level/trend of the time series.
#                 If NULL then the level/trend is calculated using CMA.
#   outplot       If TRUE provide a plot of the decomposed components.
#   decomposition Type of seasonal decomposition: "multiplicative" or "additive".
#   h             Forecast horizon for seasonal component. 
#   type          Type of calculation for seasonal component:
#                   "mean"          - The mean of each seasonal period
#                   "median"        - The median of each seasonal period
#                   "pure.seasonal" - Model using a pure seasonal model
#
# Outputs:
#   List with the following elements:
#     trend         Trend component
#     season        Seasonal component
#     irregular     Irregular component
#     f.season      Forecasted seasonal component if h>0
#     g             Purse seasonal model parameters
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
    
  ymin.s <- min(ynt)
  ymax.s <- max(ynt)
  yminmax.s <- c(ymin.s-0.1*(ymax.s-ymin.s),ymax.s+0.1*(ymax.s-ymin.s))
  ymin.y <- min(y)
  ymax.y <- max(y)
  yminmax.y <- c(ymin.y-0.1*(ymax.y-ymin.y),ymax.y+0.1*(ymax.y-ymin.y))
  
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
  if (type=="mean"){
    # Calculate the seasonality as the overall mean
    season <- colMeans(ynt, na.rm=TRUE)
    if (h>0){
      f.season <- rep(season,h+m*2)
      f.season <- as.vector(f.season[(m-ke+1):(m-ke+h)])
    } else {
      f.season <- NULL
    }
    i.season <- rep(season,ns)
    i.season <- i.season[(ks+1):(ns*m-ke)]+y*0
  }
  if (type=="median"){
    # Calculate the seasonality as the overall median
    season <- array(NA,c(1,m))
    for (si in 1:m){
      season[si] <- median(ynt[,si], na.rm=TRUE)
    }
    if (h>0){
      f.season <- rep(season,h+m*2)
      f.season <- as.vector(f.season[(m-ke+1):(m-ke+h)])
    } else {
      f.season <- NULL
    }
    i.season <- rep(season,ns)
    i.season <- i.season[(ks+1):(ns*m-ke)] + y*0
  }
  if (type=="pure.seasonal"){
    # Seasonality is modelled with a pure seasonal smoothing
    sout <- opt.sfit(ynt,"MSE",0.001,n,m)
    g <- sout$g
    if (h>0){
      season <- sout$season
      f.season <- rep(season,h+m*2)
      f.season <- as.vector(f.season[(m-k+1):(m-k+h)])
      f.season <- rep(season, h %/% m + 1)[1:h]
    } else {
      f.season <- NULL
    }
    season.sample <- matrix(t(ynt),ncol=1)        
    season.sample <- season.sample[!is.na(season.sample)]
    i.season <- as.vector(fun.sfit(season.sample,g,n,m)$ins) + y*0
  }
  
  # Convert f.season to ts object
  if (class(y) == "ts" && h>0){
    s <- end(y)
    if (s[2]==m){
      s[1] <- s[1]+1
      s[2] <- 1
    } else {
      s[2] <- s[2]+1
    }
    f.season <- ts(f.season,start=s,frequency=m)  
  } 
  
  if (decomposition == "multiplicative"){
    resid <- y - (trend*i.season)
  } else {
    resid <- y - (trend+i.season)
  }
  
  # Produce plots
  if (outplot == TRUE){
    par(mfrow=c(4,1),mar=c(0,2,0,0),oma=c(2,2,2,2))
    
    # Series
    plot(1:n,as.vector(y),type="l",xlab="",xaxt='n',ylab="",yaxt='n',xlim=c(1,n+h),ylim=yminmax.y,xaxs = "i")
    if (decomposition == "multiplicative"){
      lines(1:n,trend*i.season,col="red",lty=1)
    } else {
      lines(1:n,trend+i.season,col="red",lty=1)
    }
    if (h>0){
      polygon(c(n+1,n+1,n+h,n+h),c(yminmax.y,yminmax.y[2:1]),border=NA,col="gray93")
    }
    mtext("Data",side=2,cex=0.8,padj=-2.5)
    axis(2,cex.axis=0.8)
    legend("topleft",c("Data","Reconstructed"),lty=c(1,1),lwd=c(1,1),col=c("black","red"),cex=0.8,bty="n",horiz=TRUE)
    
    # Trend
    plot(1:n,as.vector(trend),type="l",xlab="",xaxt='n',ylab="",yaxt='n',xlim=c(1,n+h),lty=1,ylim=yminmax.y,xaxs="i")
    mtext("Trend",side=2,cex=0.8,padj=-2.5)
    axis(2,cex.axis=0.8)
    if (h>0){
      polygon(c(n+1,n+1,n+h,n+h),c(yminmax.y,yminmax.y[2:1]),border=NA,col="gray93")
    }
    
    # Season
    plot(1:n,as.vector(i.season),type="l",xlab="",xaxt='n',ylab="",yaxt='n',xlim=c(1,n+h),lty=1,ylim=yminmax.s,xaxs="i")
    mtext("Season",side=2,cex=0.8,padj=-2.5)
    axis(2,cex.axis=0.8)
    if (h>0){
      polygon(c(n+1,n+1,n+h,n+h),c(yminmax.s,yminmax.s[2:1]),border=NA,col="gray93")
      lines((n+1):(n+h),f.season,col="blue")
    }
    if (decomposition == "multiplicative"){
      lines(1:(n+h),rep(1,n+h),lty=2,col="grey")
    } else {
      lines(1:(n+h),rep(0,n+h),lty=2,col="grey")
    }
    
    # Irregular
    yminmax.i = yminmax.y-mean(yminmax.y)
    plot(1:n,as.vector(resid),type="l",xlab="",xaxt='n',ylab="",yaxt='n',xlim=c(1,n+h),lty=1,ylim=yminmax.i,xaxs="i")
    mtext("Irregular",side=2,cex=0.8,padj=-2.5)
    text(1,yminmax.i[2],paste("RMSE:",round(sqrt(mean(resid^2)),2)),pos=4,cex=0.8)
    axis(2,cex.axis=0.8)
    if (h>0){
      polygon(c(n+1,n+1,n+h,n+h),c(yminmax.i,yminmax.i[2:1]),border=NA,col="gray93")
    }
    lines(1:n,rep(0,n),lty=2,col="grey")
    
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
