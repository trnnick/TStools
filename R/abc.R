abc <- function(x,prc=c(0.2,0.3,0.5),outplot=c(TRUE,FALSE),cex.prc=0.8,...){
# ABC analysis
#
# Inputs
#   x           This can either be an array, where each column is a series, or a vector of values.
#               If x is an array of time series, then the importance of each series is calculated
#               as the mean of the observations (sales volume). Otherwise, value can be whatever 
#               quantity is provided. 
#   prc         A list of percentages to separate the items in. By default this is c(0.2,0.3,0.5),
#               but any set of percentage values can be used as long as 0<=prc[i]<=1 and sum(prc)==1.
#   outplot     If TRUE provide a visualisation of the ABC analysis result.
#   cex.prc     Font size of percentage reported in plot.
#   ...         Additional arguments can be passed to the plot.
#   
# Outputs
#   value       A vector containing the importance value of each series.
#   class       A vector containing the class membership of each series.
#   rank        A vector containing the rank of each series. Highest value is 1.
#   importance  The importance concentration of each class, as percentage of total value.
#
# Example
#   x <- abs(matrix(cumsum(rnorm(5400,0,1)),36,150))
#   abc(x,outplot=TRUE)
#
# Nikolaos Kourentzes, 2014 <nikolaos@kourentzes.com>
# Update 2016.10: Change colour scheme, accept ellipsis

  outplot <- outplot[1]
  
  if (sum(dim(x)==1)>0 | class(x)=="numeric"){
    x.mean <- x
  } else {
    x.mean <- colMeans(x, na.rm = TRUE)
  }
  # Find rank and percentage contribution of each series
  x.rank <- order(x.mean, decreasing=TRUE)
  x.sort <- x.mean[x.rank]
  x.sort <- (x.sort/sum(x.sort))*100
  
  n <- length(x.mean)         # Number of series total
  k <- length(prc)            # Number of classes
  p <- array(0,c(k,1))        # Number of series in each class
  x.ind <- array(k,c(1,n))    # Indicator for class of each series
  x.class <- array(NA,c(1,n)) # Class of each series
  nam.abc <- LETTERS[1:k]
  x.imp <- array(0,c(k,1),dimnames=list(nam.abc,"Importance"))    # Percentage importance of each class
  
  # Calculate classes
  for (i in 1:k){
    p[i] <- round(n*prc[i])
    if (i == 1){
      x.ind[x.rank <= p[i]] <- i
      x.imp[i] <- sum(x.sort[1:sum(p[1:i])])
    } else if (i != k) {
      x.ind[sum(p[1:(i-1)])<x.rank & x.rank<=sum(p[1:i])] <- i
      x.imp[i] <- sum(x.sort[1:sum(p[1:i])]) - sum(x.imp[1:(i-1)])
    } else {
      p[i] <- n - sum(p[1:(i-1)])
      x.imp[i] <- sum(x.sort[1:sum(p[1:i])]) - sum(x.imp[1:(i-1)])
    }
    x.class[x.ind==i] <- nam.abc[i]
  }
  
  # Produce plot
  if (outplot==TRUE){
    
    # Get colours
    cmp <- colorRampPalette(brewer.pal(9,"Greens")[4:7])(k)
    
    # Allow user to override plot defaults
    args <- list(...)
    if ("main" %in% names(args)){
      tmain <- args$main
      args$main <- NULL
    } else {
      tmain <- "ABC analysis"
    }
    if ("xlab" %in% names(args)){
      txlab <- args$xlab
      args$xlab <- NULL
    } else {
      txlab <- "Cumulative number of items"
    }
    if ("ylab" %in% names(args)){
      tylab <- args$ylab
      args$ylab <- NULL
    } else {
      tylab <- "Cumulative importance"
    }
    if ("xaxs" %in% names(args)){
      txaxs <- args$xaxs
      args$xaxs <- NULL
    } else {
      txaxs <- "i"
    }
    if ("yaxs" %in% names(args)){
      tyaxs <- args$yaxs
      args$yaxs <- NULL
    } else {
      tyaxs <- "i"
    }
    # Remaining of plot inputs
    args.def = list(x=NA,y=NA,xlim=c(0,100),ylim=c(0,100),xaxs=txaxs,yaxs=tyaxs,xlab=txlab,ylab=tylab,main=tmain)
    # Use do.call to merge user defined ... inputs and defaults
    do.call(plot,c(args.def,args))
    # Plot the rest
    for (i in 1:k){
      yy <- sum(x.imp[1:i])
      xx <- (sum(p[1:i])/n)*100
      if (i == 1){
        polygon(c(0,xx,xx,0),c(0,0,yy,yy),col=cmp[i])
        text(xx/2,yy/2,nam.abc[i],cex=1.2)
      } else {
        yy2 <- sum(x.imp[1:(i-1)])
        xx2 <- (sum(p[1:(i-1)])/n)*100
        polygon(c(xx2,xx,xx,0,0,xx2),c(0,0,yy,yy,yy2,yy2),col=cmp[i])
        text(xx2+(xx-xx2)/2,yy/2,nam.abc[i],cex=1.2)
      }
      text(xx/2,yy,paste(format(round(x.imp[i],1),nsmall=1),"%",sep=""),cex=cex.prc,adj=c(0.5,1))
    }
    lines(((1:n)/n)*100,cumsum(x.sort))
    lines(((1:n)/n)*100,((1:n)/n)*100,lty=2)
  }
  
  return(list("value"=x.mean,"class"=x.class,"rank"=x.rank,"importance"=x.imp))

}