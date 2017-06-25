residout <- function(resid,t=2,outplot=c(TRUE,FALSE)){
# Create a control chart of residuals and identify outliers
#
# Inputs
#   resid     Vector of residuals
#   t         Threshold value over which standardised residuals are regarded as outliers
#   outplot   If TRUE then a control chart of the standardised residuals is plotted
#
# Outputs
#   location  Location of outliers
#   outliers  Value of outliers
#   residuals Standardised residuals
#
# Example
#   # Create some random normally distributed residuals
#   resid <- rnorm(50)
#   resid.out(resid, outplot=TRUE)
#
# Nikolaos Kourentzes, 2015 <nikolaos@kourentzes.com>

  outplot <- outplot[1]
  
  # Standardise residuals
  resid <- resid/sd(resid)
  
  if (outplot == TRUE || outplot == 1){
    n <- length(resid)
    rbar <- mean(resid)
    yylim <- max(t+1,max(abs(resid)))*1.1
    # Setup plot
    plot(-n,0,type="p",ylim=c(-1,1)*yylim,xlim=c(1,n),ylab="",xlab="")
    polygon(c(0,n+1,n+1,0),rbar+(yylim+1)*c(-1,-1,1,1),col="lightgrey",border=NA)
    polygon(c(0,n+1,n+1,0),rbar+(t+1)*c(-1,-1,1,1),col="lightcyan2",border=NA)
    polygon(c(0,n+1,n+1,0),rbar+t*c(-1,-1,1,1),col="lightblue",border=NA)
    # Seperate residuals to within and out of 2 st.dev.
    resid.in <- resid.out <- resid
    resid.in[abs(resid)>=t] <- NA
    resid.out[abs(resid)<t] <- NA
    # Plot mean and residuals
    lines(c(0,n+1),rbar*c(1,1),type="l",col="lightblue4",lwd=2)
    lines(1:n,resid.in,type="p",pch=20,cex=1.5)
    lines(1:n,resid.out,type="p",col="red",pch=20,cex=1.5)
    box()
  }
  
  # Find location of residuals > t and values
  idx <- which(abs(resid) >= t)
  if (length(idx) > 0){
    res <- resid[idx]
  } else {
    idx <- NULL
    res <- NULL
  }
  
  return(list(location=idx,outliers=res,residuals=resid))
  
}
