nzlength <- function(y,t=NULL,outplot=c(FALSE,TRUE)){
# Find dataset series non-zero length and identify series over or equal to threshold
#
# Inputs:
#   y             Dataset, each column is a series.
#   t             Threshold to compare length with. Use NULL if not needed.
#   outplot       If TRUE produce plot of the cumulative distribution of lengths.
#
# Outputs:
#   data.length   Vector with non-zero length of each time series.
#   index         Logical vector of series that are above or equal to threshold.
#
# Example:
#   x <- abs(matrix(cumsum(rnorm(5400,0,1)),36,150))
#   x[1:20,1:50] <- 0
#   nzlength(x,30,outplot=TRUE)
# 
# Nikolaos Kourentzes, 2014 <nikolaos@kourentzes.com>
  
  # Default
  outplot <- outplot[1]
  
  n <- ncol(y)
  data.length <- vector("numeric",n)
  for (i in 1:n){
    data.length[i] <- sum(y[,i]>0)
  }
  
  # Find series under the threshold
  if (!is.null(t)){
    idx <- data.length >= t
  } else {
    idx <- rep(TRUE,n)
  }
  
  if (outplot == TRUE){
    sl <- sort(data.length)
    sl.u <- sl[sl<t]
    sl.a <- sl[sl>=t]
    # Plot profile
    plot(1:n,sl,type="l",col="black",lwd=0,xaxs="i",yaxs="i",xlab="Number of series",ylab="Length",ylim=c(0,max(sl)))
    if (length(sl.u)>0){
      # Plot under threshold
      polygon(c(1:length(sl.u),length(sl.u),n,n,1),c(sl.u,min(sl.a)*c(1,1),0,0),col="grey",border="black")
      # Plot over threshold
      polygon(c((length(sl.u)+1):n,n,(length(sl.u)+1)),c(sl.a,(max(sl.u)+1)*c(1,1)),col="deepskyblue",border="black")
      mtext(c(paste0(round(100*sum(!idx)/n,2),'%'),paste0(round(100*sum(idx)/n,2),'%')),side=4,line=0.2,at=c(max(sl.u)/2,((max(sl)-max(sl.u))/2+max(sl.u))))
    } else {
      polygon(c(1:n,n,1),c(sl,0,0),col="deepskyblue",border="black")
    }
    lines(c(1,n),c(t,t),col="black",lwd=1,lty=1)
  }
  
  return(list(data.length=data.length,index=idx))

}