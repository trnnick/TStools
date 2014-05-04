abcxyz <- function(imp,frc,outplot=c(TRUE,FALSE)){
# ABC-XYZ visualisation
# 
# Inputs
#   imp       Output of abc.
#   frc       output of xyz.
#
# Output
#   abcxyz    Matrix containing the number of time series in each class.
#
# Example
#   x <- abs(matrix(cumsum(rnorm(5400,0,1)),36,150))
#   abcxyz(abc(x,outplot=FALSE),xyz(x,type="cv",outplot=FALSE))
#
# Nikolaos Kourentzes, 2014 <nikolaos@kourentzes.com>
  
  outplot <- outplot[1]
  
  k.abc <- length(imp$importance)
  nam.abc <- rownames(imp$importance)
  k.xyz <- length(frc$cum.error)
  nam.xyz <- rownames(frc$cum.error)
  
  n <- length(imp$value)
  matrix.abcxyz <- matrix(0,k.abc,k.xyz,dimnames=list(nam.abc,nam.xyz))
  
  for (i in 1:k.abc){
    for (j in 1:k.xyz){
      matrix.abcxyz[i,j] <- sum(imp$class == nam.abc[i] & frc$class == nam.xyz[j])
    }
  }
  
  if (outplot==TRUE){
    cmp = rainbow(k.abc*k.xyz,alpha=0.5,start=0,end=4/6)
    x.abc <- c(0,1:k.abc)
    y.xyz <- c(0,1:k.xyz)
    plot(c(1,1),c(0,k.xyz),type="l",xlim=c(0,k.abc),ylim=c(0,k.xyz),xaxs="i",
         yaxs="i",xlab="Importance",ylab="Forecastability",xaxt="n",yaxt="n")
    for (i in 1:k.abc){
      for (j in 1:k.xyz){
        polygon(c(x.abc[i],x.abc[i+1],x.abc[i+1],x.abc[i]),
                c(y.xyz[j],y.xyz[j],y.xyz[j+1],y.xyz[j+1]),col=cmp[i+k.abc*(j-1)])
        text((x.abc[i]+x.abc[i+1])/2,(y.xyz[j]+y.xyz[j+1])/2,
             paste(nam.abc[i],nam.xyz[j],":",matrix.abcxyz[i,j],sep=""))
      }
    }
    axis(side = 1, at = c(0.5, k.abc-0.5), labels = c("High","Low"))
    axis(side = 2, at = c(0.5, k.xyz-0.5), labels = c("Low","High"))
  }
  
  return(matrix.abcxyz)
  
}