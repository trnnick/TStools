xyz <- function(x,m=NULL,prc=c(0.2,0.3,0.5),type=c("naive","ets","cv"),outplot=c(TRUE,FALSE)){
# XYZ analysis
#
# Inputs
#   x           This can either be an array, where each column is a series, or a vector of values.
#               If x is a vector of values forecastability is not calculated and the input is used as such.
#   m           Seasonal length for time series. Required when type is "naive" or "ets".
#   prc         A list of percentages to separate the items in. By default this is c(0.2,0.3,0.5),
#               but any set of percentage values can be used as long as 0<=prc[i]<=1 and sum(prc)==1.
#   type        The type of forecastability calculation
#                 "naive" - fit naive and seasonal naive and calculate fit RMSE/mean level
#                 "ets"   - fit ets and calculate fit RMSE/mean level
#                 "cv"    - use coefficient of variation as a proxy of forecastability
#   outplot     If TRUE provide a visualisation of the ABC analysis result.
#   
# Outputs
#   value       A vector containing the forecastability value of each series.
#   class       A vector containing the class membership of each series.
#   rank        A vector containing the rank of each series. Lowest forecastability is 1.
#   cum.error   The forecastability concentration of each class, as percentage of total value.
#   model       Fitted model for each series.
#
# Example
#   x <- abs(matrix(cumsum(rnorm(5400,0,1)),36,150))
#   xyz(x,m=12,outplot=TRUE)
#
# Nikolaos Kourentzes, 2014 <nikolaos@kourentzes.com>

  outplot <- outplot[1]
  type <- type[1]
  
  n <- dim(x)[2]              # Number of series total
  
  if (sum(dim(x)==1)>0){
    x.mean <- x
    x.model <- NULL
  } else {
    x.mean <- array(0,c(1,n))
    x.model <- array(NA,c(1,n))
    for (i in 1:n){
      if (type=="cv"){
        x.mean[i] <- sd(x[,i],na.rm=TRUE)/mean(x[,i],na.rm=TRUE)
        x.model[i] <- "cv"
      } else {
        if (is.null(m)){
          stop("Seasonality m not provided")
        }
        x.temp <- x[,i]
        x.temp <- x.temp[!is.na(x.temp)]
        if (type=="ets"){
          fit <- ets(ts(x.temp,frequency=m))
          x.mean[i] <- sqrt(fit$mse)/mean(x.temp)
          x.model[i] <- fit$method
        } else {
          l <- length(x.temp)
          n1.fit <- c(NA,x.temp[1:(l-1)])
          nm.fit <- c(rep(NA,m),x.temp[1:(l-m)])
          n.err <- c(mean((x.temp[(m+1):l] - n1.fit[(m+1):l])^2),
                     mean((x.temp[(m+1):l] - nm.fit[(m+1):l])^2))
          cmp.err <- (l-m)*log(n.err) + c(2,2*m)
          cmp <- cmp.err == min(cmp.err)
          x.mean[i] <- sqrt(n.err[cmp])/mean(x.temp)
          x.model[i] <- c("Naive","Seasonal Naive")[cmp]
        }
      }
    }
  }
  # Find rank and percentage contribution of each series
  x.rank <- order(x.mean, decreasing=TRUE)
  x.sort <- x.mean[x.rank]
  x.sort <- (x.sort/sum(x.sort))*100
  
  k <- length(prc)            # Number of classes
  p <- array(0,c(k,1))        # Number of series in each class
  x.ind <- array(k,c(1,n))    # Indicator for class of each series
  x.class <- array(NA,c(1,n)) # Class of each series
  x.imp <- array(0,c(k,1))    # Percentage importance of each class
  nam.xyz <- LETTERS[26:(26-k+1)]
  
  # Calculate classes
  for (i in 1:(k)){
    p[i] <- round(n*prc[i])
    if (i==1){
      x.ind[x.rank<=p[i]] <- i
      x.imp[i] <- sum(x.sort[1:sum(p[1:i])])
    } else if (i!=k) {
      x.ind[sum(p[1:(i-1)])<x.rank & x.rank<=sum(p[1:i])] <- i
      x.imp[i] <- sum(x.sort[1:sum(p[1:i])]) - sum(x.imp[1:(i-1)])
    } else {
      p[i] <- n - sum(p[1:(i-1)])
      x.imp[i] <- sum(x.sort[1:sum(p[1:i])]) - sum(x.imp[1:(i-1)])
    }
    x.class[x.ind==i] <- nam.xyz[i]
  }
  rownames(x.imp) <- nam.xyz
  
  # Produce plot
  if (outplot==TRUE){
    cmp <- rainbow(k,start=0,end=1/6,alpha=0.5)
    plot(((1:n)/n)*100,((1:n)/n)*100,xlim=c(1,100),ylim=c(0,100),type="l",xaxs="i",yaxs="i",
         xlab="",ylab="",lty=2)
    for (i in 1:k){
      yy <- sum(x.imp[1:i])
      xx <- (sum(p[1:i])/n)*100
      if (i == 1){
        polygon(c(0,xx,xx,0),c(0,0,yy,yy),col=cmp[i])
        text(xx/2,yy/2,nam.xyz[i],cex=1.2)
      } else {
        yy2 <- sum(x.imp[1:(i-1)])
        xx2 <- (sum(p[1:(i-1)])/n)*100
        polygon(c(xx2,xx,xx,0,0,xx2),c(0,0,yy,yy,yy2,yy2),col=cmp[i])
        text(xx2+(xx-xx2)/2,yy/2,nam.xyz[i],cex=1.2)
      }
      text(xx/2,yy,paste(format(round(x.imp[i],1),nsmall=1),"%",sep=""),cex=0.6,adj=c(0.5,1))
    }
    lines(((1:n)/n)*100,cumsum(x.sort))
  }
  
  return(list("value"=x.mean,"class"=x.class,"rank"=x.rank,"cum.error"=x.imp,"model"=x.model))

}