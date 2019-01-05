bias.coeff <- function(mre,outplot=c(0,1,2)){
# Calculate the bias coefficient and optionally plot it  
#
# Inputs
#   mre       k Root Errors
#   outplot   Optional plot output: 0 - No; 1 - Histogram; 2 - Boxplot
#
# Output
#   bias      k Bias Coefficients
# 
# Example
#   # Create some random MRE
#   mre <- runif(10,0,5) + runif(10,0,5)*1i
#   bias.coeff(mre,outplot=2)
#
# Reference
#   Kourentzes N., Trapero J. R., Svetunkov I. Measuring the behaviour of experts on demand 
#   forecasting: a complex task. Working paper
#   http://kourentzes.com/forecasting/2014/12/17/measuring-the-behaviour-of-experts-on-demand-forecasting-a-complex-task/
#
# Nikolaos Kourentzes, 2014 <nikolaos@kourentzes.com>
  

  outplot <- outplot[1]
  
  gamma <- Arg(mre)
  bias <- 1 - 4*gamma/pi
  k <- length(bias)
  if (outplot == 1){
    if (k>1){
      x <- hist(bias,plot=FALSE)
      ymax <- ceiling(max(x$counts)*1.1)
      plot(x,xlim=c(-1,1),col="lightblue",xlab="Bias Coefficient",
          ylim=c(0,ymax),main=NULL)
    } else {
      ymax <- 1.1
      plot(x=bias,y=1,type="h",lwd=6,col="lightblue",xlab="Bias Coefficient",
           ylim=c(0,ymax),xlim=c(-1,1),main=NULL)
    }
    lines(c(0,0),c(0,ymax),col="black",lwd=3)
    text(-0.5,ymax,"Negative bias",cex=0.95)
    text(0.5,ymax,"Positive bias",cex=0.95)
  }
  if (outplot == 2){
    boxplot(bias,horizontal=TRUE,ylim=c(-1,1),col="lightblue",xlab="Bias Coefficient")
    lines(c(0,0),c(0,2),col="black",lwd=3)
    text(0.5,1.5,'--')
    text(-0.5,1.5,'--')
    text(-0.5,1.5,"\n\nNegative bias",cex=0.95)
    text(0.5,1.5,"\n\nPositive bias",cex=0.95)
    text(-0.75,1.5,"Strong",cex=0.95)
    text(-0.25,1.5,"Weak",cex=0.95)
    text(0.75,1.5,"Strong",cex=0.95)
    text(0.25,1.5,"Weak",cex=0.95)
  }
  
  return(bias)
  
}

mre <- function(e,op=c("mean","sum","gm")){
# Calculate Root Error  
# Given k errors mre calculates the Mean Root Error, the Sum Root Error 
# or the Mod(Geometric Squared Mean Root Error) = GRMSE
#
# Inputs
#   e       k forecast error(s)
#   op      Aggregation operator of the Root Error: 1) Mean; 2) Sum; 3) Magnitude of Geometric Mean
#
# Output
#   mre     Root Error
#
# Example
#   # Generate some random errors
#   e <- runif(10,-10,10)
#   mre(e)
#
# Reference
#   Kourentzes N., Trapero J. R., Svetunkov I. Measuring the behaviour of experts on demand 
#   forecasting: a complex task. Working paper
#   http://kourentzes.com/forecasting/2014/12/17/measuring-the-behaviour-of-experts-on-demand-forecasting-a-complex-task/
#
# Nikolaos Kourentzes, 2014 <nikolaos@kourentzes.com>

  op <- op[1]
  
  e <- as.complex(e)
  e <- sqrt(e)
  switch(op,
         "mean" = {mre <- mean(e)},
         "sum" = {mre <- sum(e)},
         "gm" = {mre <- Mod(prod(e)^(2/length(e)))})
  
  return(mre)
  
}

mre.plot <- function(mre,main=NULL,plot.legend=c(TRUE,FALSE),maxmag=NULL,ts=1,group=NULL){
# Mean Root Error Bias Plot
# Plots the `Bias Plot' for MRE.
#
# Inputs
#   mre         k Root or Mean Root Errors to be plotted. These must already be complex numbers.
#   main        Main title for plot
#   plot.legend Plot legend if k > 1. Default is TRUE.
#   maxmag      Maximum MRE magnitude.
#   ts          Text size cex.
#   group       Vector containing either characters or numerics to group errors.
#
# Example
#   # Create some random MRE
#   mre <- runif(10,0,5) + runif(10,0,5)*1i
#   plot.mre(mre)
#
# Reference
#   Kourentzes N., Trapero J. R., Svetunkov I. Measuring the behaviour of experts on demand 
#   forecasting: a complex task. Working paper
#   http://kourentzes.com/forecasting/2014/12/17/measuring-the-behaviour-of-experts-on-demand-forecasting-a-complex-task/
#
# Nikolaos Kourentzes, 2014 <nikolaos@kourentzes.com>

  # Default options
  plot.legend <- plot.legend[1]
  ts=1 
  
  # Initialise
  if (is.null(maxmag)){
    sz <- Mod(mre)  
  } else {
    sz <- maxmag
  }
  max.sz <- max(sz)*1.05
  rr <- seq(pi/4,3*pi/4,length.out=200)
  bc <- c(1,.5,0,-.5,-1)
  xx <- seq(0,max.sz,round(max.sz/10,1))
  xx <- xx[2:10]
  # Produce plot layout
  plot(NA,NA,ylim=c(0,sin(pi/2)*max.sz+max.sz*.125),xlim=c(cos(3*pi/4)*max.sz-max.sz*.05,cos(pi/4)*max.sz+max.sz*.05),
       axes=FALSE,xlab="",ylab="",xaxs="i",yaxs="i",main=main,cex.main=ts)
  for (i in 1:5){
    if (i == 1 || i == 5){
      cc <- "black"
    } else {
      cc <- "darkgrey"
    }
    lines(c(0,cos((pi/4)+((i-1)*pi/8)))*max.sz,c(0,sin((pi/4)+((i-1)*pi/8)))*max.sz,col=cc)  
    text(cos((pi/4)+((i-1)*pi/8))*max.sz,sin((pi/4)+((i-1)*pi/8))*max.sz,bc[i],pos=3,cex=ts)
  }
  for (i in seq(2,8,2)){
    lines(c(0,cos((pi/4)+((i-1)*pi/16)))*max.sz,c(0,sin((pi/4)+((i-1)*pi/16)))*max.sz,col="darkgrey",lty=2)  
  }
  lines(cos(rr)*max.sz,sin(rr)*max.sz,col="black")
  for (i in 1:9){
    lines(cos(rr)*xx[i],sin(rr)*xx[i],col="darkgrey",lty=2)
    text(cos(pi/4)*xx[i],sin(pi/4)*xx[i],xx[i],pos=4,cex=ts)
    text(cos(3*pi/4)*xx[i],sin(3*pi/4)*xx[i],xx[i],pos=2,cex=ts)
  } 
  text(cos(pi/4)*max.sz/2,sin(pi/4)*max.sz/2,"\n\n\n\nError",srt=45,cex=ts)
  text(cos(3*pi/4)*max.sz/2,sin(3*pi/4)*max.sz/2,"\n\n\n\nError",srt=-45,cex=ts)
  text(0,sin(pi/2)*max.sz,"Bias coefficient\n",pos=3,cex=ts)
  text(cos(3*pi/8)*max.sz*1.1,sin(3*pi/8)*max.sz*1.1,"Positive bias\n",cex=ts)
  text(cos(5*pi/8)*max.sz*1.1,sin(5*pi/8)*max.sz*1.1,"Negative bias\n",cex=ts)
  # Plot mre
  k <- length(mre)
  # Get colours
  if (is.null(group)){
    cmp <- colorRampPalette(RColorBrewer::brewer.pal(9,"YlGnBu")[2:8])(k)
    cmp2 <- RColorBrewer::brewer.pal(3,"Set1")[1]
  } else {
    uc <- table(group)
    un <- as.vector(uc)
    uc <- names(uc)
    u <- length(un)
    plts <- c("Blues","Greens","Reds","Purples","RdPu")
    if (u > 5){
      plts <- rep(plts,ceiling(u/5))
    }
    cmp <- vector("character",sum(un))
    for (i in 1:u){
      cmp[group==uc[i]] <- colorRampPalette(RColorBrewer::brewer.pal(9,plts[i])[3:7])(un[i])
    }
    cmpl <- unlist(lapply(plts[1:u],function(x){RColorBrewer::brewer.pal(5,x)[5]}))
    cmpo <- unlist(lapply(plts[1:u],function(x){RColorBrewer::brewer.pal(3,x)[2]}))
    cmp2 <- "black"
  }
  # cmp <- rainbow(k,start=3.4/6,end=4.4/6)
  for (i in 1:k){
    temp.gamma <- Arg(mre[i])+pi/4
    temp.r <- Mod(mre[i])
    lines(c(0,cos(temp.gamma)*temp.r),c(0,sin(temp.gamma)*temp.r),col=cmp[i],lwd=2)
    lines(cos(temp.gamma)*temp.r,sin(temp.gamma)*temp.r,type="o",col=cmp[i],pch=20,cex=1.5)
  }
  # Plot overall behaviour
  if (k>1){
    if (!is.null(group)){
      for (i in 1:u){
        mean.mre <- mean(mre[group==uc[i]])
        temp.gamma <- Arg(mean.mre)+pi/4
        temp.r <- Mod(mean.mre)
        lines(c(0,cos(temp.gamma)*temp.r),c(0,sin(temp.gamma)*temp.r),col=cmpo,lwd=3)
        lines(cos(temp.gamma)*temp.r,sin(temp.gamma)*temp.r,type="o",col=cmp2,pch=20,cex=2.5)
      }
    }
    mean.mre <- mean(mre)
    temp.gamma <- Arg(mean.mre)+pi/4
    temp.r <- Mod(mean.mre)
    lines(c(0,cos(temp.gamma)*temp.r),c(0,sin(temp.gamma)*temp.r),col=cmp2,lwd=3)
    lines(cos(temp.gamma)*temp.r,sin(temp.gamma)*temp.r,type="o",col=cmp2,pch=20,cex=2.5)
    if (plot.legend == TRUE){
      if (is.null(group)){
        legend("bottomleft",c("MRE","Overall"),lty=1,col=c(cmp[round(k/2)],cmp2),bty="n",pch=20,lwd=2,cex=ts)
      } else {
        legend("bottomleft",c(uc,paste(uc,"overall"),"Overall"),lty=1,col=c(cmpl,cmpo,cmp2),pt.bg=c(cmpl,rep(cmp2,u),cmp2),bty="n",pch=21,lwd=2,cex=ts)
      }
    }
  }
}
