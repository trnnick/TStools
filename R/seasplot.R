seasplot <- function(y,m=NULL,s=NULL,trend=NULL,colour=NULL,alpha=0.05,
                     outplot=c(1,0,2,3,4,5),decomposition=c("multiplicative","additive"),
                     cma=NULL,labels=NULL,...)
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
#   decomposition Type of seasonal decomposition: "multiplicative" or "additive".
#   cma           Input precalculated level/trend for the analysis. Overrides trend=NULL.
#   labels        External labels for the seasonal periods. Use NULL for default. 
#                 If length(labels) < m, then this input is ignored.
#   ...           Additional arguments can be passed to plotting functions. For example use 
#                 main="" to replace the title.
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
  decomposition <- match.arg(decomposition,c("multiplicative","additive"))
  outplot <- outplot[1]
  
  # Get m (seasonality)
  if (is.null(m)){
    if (any(class(y) == "ts")){
      m <- frequency(y)
    } else {
      stop("Seasonality not defined (y not ts object).")
    }
  }
  
  # Get starting period in seasonality if available
  if (is.null(s)){
    if (any(class(y) == "ts")){
      s <- start(y)
      s <- s[2]
      # Temporal aggregation can mess-up s, so override if needed
      if (is.na(s)){s<-1}
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
    yminmax <- c(-1,1)*max(abs(ymid-yminmax))+ymid
    if (is.null(season.pval)){
      plottitle <- paste(title.trend, "\n", title.season,sep="")
    } else {
      plottitle <- paste(title.trend, "\n", title.season,
                   " (p-val: ",round(season.pval,3),")",sep="")
    }
    # Allow user to override plot defaults
    args <- list(...)
    if (!("main" %in% names(args))){
      addtitle <- TRUE
    } else {
      addtitle <- FALSE
    }
    if (!("xlab" %in% names(args))){
      args$xlab <- "Period"
    }
    if (!("ylab" %in% names(args))){
      args$ylab <- ""
    }
    if (!("yaxs" %in% names(args))){
      args$yaxs <- "i"
    }
    if (!("ylim" %in% names(args))){
      args$ylim <- yminmax
    }
    # Remaining defaults
    args$x <- args$y <- NA
   }
  
  if (outplot == 1){
    # Conventional seasonal diagramme
    if (is.null(colour)){
      cmp <- colorRampPalette(brewer.pal(9,"YlGnBu")[4:8])(ns)
    } else {
      cmp <- rep(colour,times=ns)
    }
    # Additional plot options
    if (!("xlim" %in% names(args))){
      args$xlim <- c(1,m)
    }
    if (!("xaxs" %in% names(args))){
      args$xaxs <- "i"
    }
    if (addtitle){
      args$main <- paste0("Seasonal plot ",main=plottitle)
    }
    args$xaxt <- "n"
    # Produce plot
    do.call(plot,args)
    for (i in 1:ns){
      lines(ynt[i,],type="l",col=cmp[i])
    }
    lines(c(0,m+1),c(ymid,ymid),col="black",lty=2)
    legend("topleft",c("Oldest","Newest"),col=c(cmp[1],cmp[ns]),lty=1,bty="n",lwd=2,cex=0.7)
    axis(1,at=1:m,labels=labels)
  } 
  if (outplot == 2){
    # Seasonal boxplots
    if (is.null(colour)){
      cmp <- brewer.pal(3,"Set3")[1]
    } else {
      cmp <- colour
    }
    # Additional plot options
    if (!("xlim" %in% names(args))){
      args$xlim <- c(1,m)
    }
    if (addtitle){
      args$main <- paste0("Seasonal boxplot ",main=plottitle)
    }
    args$x <- ynt
    args$col <- cmp
    # Produce plot
    do.call(boxplot,args)
    lines(c(0,m+1),c(ymid,ymid),col="black",lty=2)
  }
  if (outplot == 3){
    # Subseries plots
    if (is.null(colour)){
      cmp <- brewer.pal(3,"Set1")
    } else {
      cmp <- colour
    }
    # Additional plot options
    if (!("xlim" %in% names(args))){
      args$xlim <- c(1,m*ns)
    }
    if (addtitle){
      args$main <- paste0("Seasonal subseries ",main=plottitle)
    }
    if (!("xaxs" %in% names(args))){
      args$xaxs <- "i"
    }
    args$xaxt <- "n"
    # Produce plot
    do.call(plot,args)
    lines(c(1,ns),median(ynt[,1],na.rm=TRUE)*c(1,1),col=cmp[1],lwd=2)
    for (i in 1:m){
      lines((1+(i-1)*ns):(ns+(i-1)*ns),ynt[,i],type="o",col=cmp[2],pch=20,cex=0.75)
      lines(c((1+(i-1)*ns),(ns+(i-1)*ns)),median(ynt[,i],na.rm=TRUE)*c(1,1),col=cmp[1],lwd=2)
      if (i < m){
        lines(c(1,1)*i*ns+0.5,yminmax,col="gray")
      }
    }
    lines(c(0,m*ns+1),c(ymid,ymid),col="black",lty=2)
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
    # Additional plot options
    if (!("xlim" %in% names(args))){
      args$xlim <- c(1,m)
    }
    if (addtitle){
      args$main <- paste0("Seasonal distribution ",main=plottitle)
    }
    if (!("xaxs" %in% names(args))){
      args$xaxs <- "i"
    }
    args$xaxt <- "n"
    # Produce plot
    do.call(plot,args)
    polygon(c(1:m,rev(1:m)),c(qntl[7,],rev(qntl[1,])),col=gray(0.8),border=NA)
    polygon(c(1:m,rev(1:m)),c(qntl[6,],rev(qntl[2,])),col="lightblue",border=NA)
    polygon(c(1:m,rev(1:m)),c(qntl[5,],rev(qntl[3,])),col="skyblue",border=NA)
    lines(1:m,qntl[4,],col=cmp,lwd=2)
    lines(c(0,m*ns+1),c(ymid,ymid),col="black",lty=2)
    legend("topleft",c("Median","25%-75%","10%-90%","MinMax"),col=c(cmp,"skyblue","lightblue",gray(0.8)),lty=1,bty="n",lwd=2,cex=0.7)
    axis(1,at=1:m,labels=labels)
    box()
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
    cmp <- c(rgb(0,0,0,0), colorRampPalette((brewer.pal(9,"Blues")[3:9]))(100))
    # Additional plot options
    if (addtitle){
      args$main <- paste0("Seasonal density ",main=plottitle)
    }
    args$xaxt <- "n"
    args$col <- cmp
    args$x <- 1:m
    args$y <- llc
    args$z <- dnst
    # Produce plot
    do.call(image,args)
    lines(colMeans(ynt,na.rm=TRUE),type="o",lty=1,bg=brewer.pal(3,"Set1")[1],pch=21,cex=1.1,lwd=2)
    lines(c(0,m*ns+1),c(ymid,ymid),col="black",lty=2)
    box()
    axis(1,at=1:m,labels=labels)
  }
  
  out <- list(season=ynt,season.exist=season.exist,season.pval=season.pval,
              trend=cma,trend.exist=trend.exist,trend.pval=trend.pval)
  out <- structure(out,class="classseas")
  return(out)

}

summary.classseas <- function(object,...){
  print(object)
}

print.classseas <- function(x,...){
  cat("Results of statistical testing\n")
  if (!is.null(x$trend.exist)){
    cat(paste0("Evidence of trend: ",x$trend.exist, "  (pval: ",round(x$trend.pval,3),")\n"))
  } else {
    cat("Presence of trend not tested.\n")
  }
  cat(paste0("Evidence of seasonality: ",x$season.exist, "  (pval: ",round(x$season.pval,3),")\n"))
}
