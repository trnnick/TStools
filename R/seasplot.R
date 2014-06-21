seasplot <- function(y,m=NULL,s=NULL,trend=NULL,colour=NULL,alpha=0.05,
                     outplot=c(1,0,2,3,4,5),decomposition=c("multiplicative","additive"),
                     cma=NULL,labels=NULL)
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
#   decomposition type of seasonal decomposition: "multiplicative" or "additive".
#   cma           Input precalculated level/trend for the analysis. Overrides trend=NULL.
#   labels        External labels for the seasonal periods. Use NULL for default. 
#                 If length(labels) < m, then this input is ignored.
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
  decomposition <- decomposition[1]
  outplot <- outplot[1]
  
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
  
  # Make sure that labels input is fine
  if (!is.null(labels)){
    if (length(labels) < m){
      labels  <- NULL
    } else {
      labels <- labels[1:m]
    }
  }
  
  if (is.null(labels)){
    labels <- paste(1:m)
  }
  
  n <- length(y)
  
  if ((decomposition == "multiplicative") && (min(y)<=0)){
    decomposition <- "additive"
  }
  
  # Override trend if cma is given
  if (!is.null(cma)){
    trend <- NULL
    if (n != length(cma)){
      stop("Length of series and cma do not match.")
    }
  }
  
  # Calculate CMA
  if ((is.null(trend) || trend == TRUE) && (is.null(cma))){
      cma <- cmav(y=y,ma=m,fill=TRUE,outplot=FALSE)
  }
  
  # Test for changes in the CMA (trend)
  if (is.null(trend)){
    trend.pval <- coxstuart(cma)$p.value 
    trend.exist <- trend.pval <= alpha/2
    trend <- trend.exist
  } else {
    trend.exist <- NULL
    trend.pval <- NULL
  }
  
  if (trend == TRUE){
    if (decomposition == "multiplicative"){
      ynt <- y/cma  
    } else {
      ynt <- y - cma
    }
    title.trend <- "(Detrended)"
  } else {
    ynt <- y
    title.trend <- ""
    cma <- NULL
  }
  
  ymin <- min(ynt)
  ymax <- max(ynt)
  ymid <- median(ynt)
  
  # Fill with NA start and end of season
  k <- m - (n %% m)
  ks <- s-1
  ke <- m - ((n+ks) %% m)
  ynt <- c(rep(NA,times=ks),as.vector(ynt),rep(NA,times=ke))
  ns <- length(ynt)/m
  ynt <- matrix(ynt,nrow=ns,ncol=m,byrow=TRUE)
  colnames(ynt) <- labels
  rownames(ynt) <- paste("s",1:ns,sep="")
  
  # Check seasonality with Friedman
  if (m>1 && (length(y)/m)>=2){
    season.pval <- friedman.test(ynt)$p.value
    season.exist <- season.pval <= alpha
    if (season.exist==TRUE){
      title.season <- "Seasonal"
    } else {
      title.season <- "Nonseasonal"
    }
  } else {
    season.pval <- NULL
    season.exist <- NULL
    title.season <- "Nonseasonal"
  }
  
  # Produce plots
  if (outplot != 0){
    yminmax <- c(ymin - 0.1*(ymax-ymin),ymax + 0.1*(ymax-ymin))
    if (is.null(season.pval)){
      plottitle <- paste(title.trend, "\n", title.season,sep="")
    } else {
      plottitle <- paste(title.trend, "\n", title.season,
                   " (p-val: ",round(season.pval,3),")",sep="")
    }
  }
  
  if (outplot == 1){
    # Conventional seasonal diagramme
    if (is.null(colour)){
      cmp <- rainbow(ns,start=2.5/6,end=4.5/6)
      cmp <- rev(cmp)
    } else {
      cmp <- rep(colour,times=ns)
    }
    plot(ynt[1,],type="l",col=cmp[1],xlab="Period",ylab="",ylim=yminmax,xlim=c(1,m),xaxt="n")
    for (i in 2:ns){
      lines(ynt[i,],type="l",col=cmp[i])
    }
    lines(c(0,m+1),c(ymid,ymid),col="black",lty=2)
    title(paste("Seasonal diagramme ",main=plottitle,sep=""))
    legend("topleft",c("Oldest","Newest"),col=c(cmp[1],cmp[ns]),lty=1,bty="n",lwd=2,cex=0.7)
    axis(1,at=1:m,labels=labels)
  } 
  if (outplot == 2){
    # Seasonal boxplots
    if (is.null(colour)){
      cmp <- "cyan"
    }
    boxplot(ynt,col=cmp,ylim=yminmax,xlim=c(1,m),xlab="Period")
    lines(c(0,m+1),c(ymid,ymid),col="black",lty=2)
    title(paste("Seasonal boxplot ",main=plottitle,sep=""))
  }
  if (outplot == 3){
    # Subseries plots
    if (is.null(colour)){
      cmp <- array(NA,c(1,2))
      cmp[1] <- "blue"
      cmp[2] <- "red"
    }
    plot(1:ns,ynt[,1],type="o",col=cmp[1],xlab="Period",ylab="",ylim=yminmax,
         xlim=c(1,m*ns),pch=20,cex=0.75,xaxt="n")
    lines(c(1,ns),median(ynt[,1],na.rm=TRUE)*c(1,1),col=cmp[2],lwd=2)
    lines(c(1,1)*ns+0.5,yminmax,col="gray")
    for (i in 2:m){
      lines((1+(i-1)*ns):(ns+(i-1)*ns),ynt[,i],type="o",col=cmp[1],pch=20,cex=0.75)
      lines(c((1+(i-1)*ns),(ns+(i-1)*ns)),median(ynt[,i],na.rm=TRUE)*c(1,1),col=cmp[2],lwd=2)
      if (i < m){
        lines(c(1,1)*i*ns+0.5,yminmax,col="gray")
      }
    }
    lines(c(0,m*ns+1),c(ymid,ymid),col="black",lty=2)
    title(paste("Seasonal subseries ",main=plottitle,sep=""))
    axis(1,at=seq(0.5+(ns/2),0.5+m*ns-ns/2,ns),labels=labels)
  }
  if (outplot == 4){
    # Seasonal distribution
    if (is.null(colour)){
      cmp <- "royalblue"
    }
    qntl <- matrix(NA,nrow=9,ncol=m)
    for (i in 1:m){
      qntl[,i] <- quantile(ynt[!is.na(ynt[,i]),i], c(0, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 1))
    }
    plot(x=NULL,y=NULL,xlim=c(1,m),ylim=yminmax,ylab="",xlab="Period",xaxt="n")
    polygon(c(1:m,rev(1:m)),c(qntl[7,],rev(qntl[1,])),col=gray(0.8),border=NA)
    polygon(c(1:m,rev(1:m)),c(qntl[6,],rev(qntl[2,])),col="lightblue",border=NA)
    polygon(c(1:m,rev(1:m)),c(qntl[5,],rev(qntl[3,])),col="skyblue",border=NA)
    lines(1:m,qntl[4,],col=cmp,lwd=2)
    lines(c(0,m*ns+1),c(ymid,ymid),col="black",lty=2)
    legend("topleft",c("Median","25%-75%","10%-90%","MinMax"),col=c(cmp,"skyblue","lightblue",gray(0.8)),lty=1,bty="n",lwd=2,cex=0.7)
    title(paste("Seasonal distribution ",main=plottitle,sep=""))
    axis(1,at=1:m,labels=labels)
  }
  if (outplot == 5){
    dnst <- matrix(NA,nrow=m,ncol=512)
    for (i in 1:m){
       tmp <- density(ynt[!is.na(ynt[,i]),i], bw = "SJ", n = 512, from = yminmax[1], to = yminmax[2])
       dnst[i,] <- tmp$y
       if (i == 1){
        llc <- tmp$x
       }
    }
    image(1:m,llc,dnst,ylim=yminmax,xlab="Period",ylab="",
          col=c(rgb(0,0,0,0),colorRampPalette(c(rgb(1,1,1,0), rgb(0,0,1,0), rgb(0,1/3,1,0), rgb(0,2/3,1,0), rgb(0,1,1,0)))(100)),
          xaxt="n")
    # filled.contour(1:m,llc,dnst,ylim=yminmax,xlab="Period",ylab="",nlevels=50,color.palette=rainbow)
    lines(colMeans(ynt,na.rm=TRUE),type="o",lty=1,bg="skyblue",pch=21,cex=0.7)
    lines(c(0,m*ns+1),c(ymid,ymid),col="black",lty=2)
    title(paste("Seasonal density ",main=plottitle,sep=""))
    axis(1,at=1:m,labels=labels)
  }
  
  return(list(season=ynt,season.exist=season.exist,season.pval=season.pval,
              trend=cma,trend.exist=trend.exist,trend.pval=trend.pval))

}