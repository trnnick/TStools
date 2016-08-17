nemenyi <- function(data, conf.int=0.95, sort=c(TRUE,FALSE), 
                    plottype=c("vline","none","mcb","vmcb","line"), 
                    console=c(FALSE,TRUE), silent=c(FALSE,TRUE), 
                    pcol=NULL, title=NULL, select=NULL, labels=NULL, ...)
{
# Friedman and Nemenyi tests for nonparametric comparisons between methods.
#
# Inputs:
#  data       A matrix that includes measures of accuracy for several methods (in columns) 
#             for several time series (rows), of size n x k.
#  conf.int   The width of the confidence interval. Default is 0.95.
#  sort       If TRUE function sorts the final values of mean ranks. If plots are request 
#             this is forced TRUE.
#  plottype   Type of plot to produce:
#               - none: No plot
#               - mcb: MCB style plot
#               - vmcb: Vertical MCB style plot
#               - line: Line visualisation (ISF style). Numbers next to method names are the mean rank
#                 see: http://www.forecasters.org/proceedings12/HibonMicheleISF2012.pdf
#               - vline: Vertical line visualisation [Default]
#  pcol       Override default colours with single colour. Default is to use blue for MCB plots
#             and rainbow colours for line plots. If no evidence of significant differences is found, 
#             based on Friedman's test, then all colours become grey
#  console    If FALSE, no console output is given.
#  silent     If TRUE, console is overriden to FALSE and plottype to "none". Use to have only numerical
#             output.
#  title      Override default title
#  select     Highlight selected model (column). Number 1 to k. Select NULL for no highlighting.
#  labels     Optional labels for models. If NULL column names of 'data' will be used.
#             If number of labels < k, then it is assumed as NULL. 
#  
# Outputs:
#   means - mean rank of each method
#   intervals - nemenyi intervals for each method
#   fpval - friedman test p-value
#   cd - critical distance for nemenyi
#
# Example:
#  N <- 50
#  M <- 4
#  data <- matrix( rnorm(N*M,mean=0,sd=1), N, M) 
#  data[,2] <- data[,2]+1
#  data[,3] <- data[,3]+0.7
#  data[,4] <- data[,4]+0.5
#  colnames(data) <- c("Method A","Method B","Method C - long name","Method D")
#  out <- nemenyi(data,conf.int=0.95,plottype="vline")
#  
# Ivan Svetunkov, Nikolaos Kourentzes, 2014. v1.4

  if(length(dim(data))!=2){
    stop("Data is of the wrong dimension! It should be a table with methods in columns and observations in rows.",call.=FALSE);
  }
  
# Defaults
  sort <- sort[1]
  plottype <- plottype[1]
  console <- console[1]
  silent <- silent[1]
  
  data <- na.exclude(data)
  rows.number <- nrow(data)
  cols.number <- ncol(data)

# Check select < cols.number
  if (!is.null(select)){
    if (select > cols.number){
      select <- NULL
    }
  }

# Check if silent is requested
  if (silent == TRUE){
    plottype = "none"
    console = FALSE
  }
  
# If plot is asked, always sort the results
  if (plottype != "none"){
    sort <- TRUE
  }

# Checks for labels
  if (!is.null(labels)){
    if (length(labels)<cols.number){
      labels <- NULL
    } else {
      labels <- labels[1:cols.number]
    }
  }

# First run Friedman test. If insignificant then ignore Nemenyi (Weaker)
  fried.pval <- friedman.test(data)$p.value
  if (fried.pval <= 1-conf.int){
    fried.H <- "Different" # At least one method is different
  } else {
    fried.H <- "Same"      # No evidence of differences between methods
  }
  
# Matrix of ranks
  ranks.matrix <- matrix(NA, nrow=rows.number, ncol=cols.number)
  colnames(ranks.matrix) <- colnames(data)

# Vector of mean values of ranks
  ranks.means <- c(1:cols.number)

# Nemenyi critical distance and bounds of intervals
  r.stat <- qtukey(conf.int,cols.number,Inf) * sqrt((cols.number*(cols.number+1))/(12*rows.number))
  r.stat <- c(-r.stat,r.stat)

# The matrix of intervals for mean ranks
  ranks.intervals <- matrix(NA, nrow=2, ncol=cols.number)
  colnames(ranks.intervals) <- colnames(data)
  
# Rank methods for each time series
  for (i in 1:rows.number){
    ranks.matrix[i, ] = rank(data[i,],na.last="keep",ties.method="average")
  }

# Calculate mean values
  ranks.means <- colMeans(ranks.matrix)

# Calculate intervals for each of the methods
  for (i in 1:cols.number){
    ranks.intervals[, i] = ranks.means[i]+r.stat
  }

# Sort interval matrix and means
  if(sort==TRUE){
    order.idx <- order(ranks.means)
    ranks.means <- ranks.means[order.idx]
    ranks.intervals <- ranks.intervals[,order(ranks.intervals[2,])]
    if (!is.null(labels)){
      labels <- labels[order.idx]
    }
    if (!is.null(select)){
      select <- which(order.idx == select)
    }
  }

# Produce Console Output
  if(console==TRUE){
    writeLines("Friedman and Nemenyi Tests")
    writeLines(paste("The significance level is ", (1-conf.int)*100, "%", sep=""))
    writeLines(paste("Number of observations is ", rows.number, " and number of methods is ", cols.number, sep=""))
    writeLines(paste("Friedman test p-value: ", format(round(fried.pval,4),nsmall=4) , " - ", fried.H, sep=""))
    writeLines(paste("Nemenyi critical distance: ", format(round(r.stat[2],4),nsmall=4), sep=""))
  }
    
# Create title for plots
if(is.null(title)){
  title <- paste("Friedman: ", format(round(fried.pval,3),nsmall=3), " (", fried.H, ") \n Nemenyi CD: ", format(round(r.stat[2],3),nsmall=3), sep = "")
}

# Labels for plots
if (is.null(labels)){
  labels <- colnames(ranks.intervals)
}

# Produce plots
  # MCB style plot
  if(plottype == "mcb"){
    # Choose colour depending on Friedman test result
    if (fried.pval > 1-conf.int){pcol <- "gray"} else {if (is.null(pcol)){pcol <- "blue"}}
    # Find min max
    ymax <- max(ranks.intervals)
    ymin <- min(ranks.intervals)
    ymax <- ymax + 0.1*(ymax-ymin)
    ymin <- ymin - 0.1*(ymax-ymin)
    # Plot
    plot(1:cols.number,ranks.means,
        xlab="", ylab="Mean ranks", main=title,
        type="p", axes=FALSE, pch=20, lwd=4, col="black", ylim=c(ymin,ymax), ...)
    if (!is.null(select)){
      lines(c(0,cols.number+1),rep(ranks.intervals[1,select],times=2), col="gray", lty=2)
      lines(c(0,cols.number+1),rep(ranks.intervals[2,select],times=2), col="gray", lty=2)    
    } else {
      lines(c(0,cols.number+1),rep(ranks.intervals[1,1],times=2), col="gray", lty=2)
      lines(c(0,cols.number+1),rep(ranks.intervals[2,1],times=2), col="gray", lty=2)    
    }
    for (i in 1:cols.number){
      lines(rep(i,times=2),ranks.intervals[,i], type="b", lwd=2, col=pcol);
    }
    axis(1,at=c(1:cols.number),labels=labels);
    axis(2);
    box(which="plot", col="black");
  }

  # MCB style plot - vertical
  if(plottype == "vmcb"){
    # Choose colour depending on Friedman test result
    if (fried.pval > 1-conf.int){pcol <- "gray"} else {if (is.null(pcol)){pcol <- "blue"}}
    # Find min max
    xmax <- max(ranks.intervals)
    xmin <- min(ranks.intervals)
    xmax <- xmax + 0.1*(xmax-xmin)
    xmin <- xmin - 0.1*(xmax-xmin)
    # Find max label size
    lbl.size <- nchar(labels)
    lbl.size <- max(lbl.size)
    # Produce plot
    par(mar=c(2, lbl.size/1.5, 4, 2) + 0.1)
    plot(ranks.means,1:cols.number,
         ylab="", xlab="Mean ranks", main=title,
         type="p", axes=FALSE, pch=20, lwd=4, col="black", xlim=c(xmin,xmax), ...)
    if (!is.null(select)){
      lines(rep(ranks.intervals[1,select],times=2), c(0,cols.number+1), col="gray", lty=2)
      lines(rep(ranks.intervals[2,select],times=2), c(0,cols.number+1), col="gray", lty=2)    
    } else {
      lines(rep(ranks.intervals[1,1],times=2), c(0,cols.number+1), col="gray", lty=2)
      lines(rep(ranks.intervals[2,1],times=2), c(0,cols.number+1), col="gray", lty=2)    
    }
    for (i in 1:cols.number){
      lines(ranks.intervals[,i], rep(i,times=2), type="b", lwd=2, col=pcol);
    }
    axis(2,at=c(1:cols.number),labels=labels,las=2)
    axis(1);
    box(which="plot", col="black");
    par(mar=c(5, 4, 4, 2) + 0.1)
  }

  # Line style plot (as in ISF)
  if(plottype == "line"){
    # Find groups
    rline <- matrix(NA, nrow=cols.number, ncol=2)
    for (i in 1:cols.number){
      tloc <- which((abs(ranks.means-ranks.means[i])<r.stat[2]) == TRUE)
      rline[i,] <- c(min(tloc),max(tloc))
    }
    # Get rid of duplicates and single member groups
    rline <- unique(rline)
    rline <- rline[apply(rline,1,min) != apply(rline,1,max),]
    # Re-convert to matrix if necessary and find number of remaining groups
    if (length(rline)==2){
      rline <- as.matrix(rline)
      rline <- t(rline)
    }
    k <- nrow(rline)
    # Choose colour depending on Friedman test result
    if (fried.pval > 1-conf.int){pcol <- rep("gray",times=k)} else {
      if (is.null(pcol)){pcol <- rainbow(k)} else {pcol <- rep(pcol,times=k)}
    }
    # Prepare method labels and add mean rank to them
    lbl <- labels
    lblm <- matrix(NA,nrow=1,ncol=cols.number)
    lbl.size <- matrix(NA,nrow=1,ncol=cols.number)
    for (i in 1:cols.number){
      if (is.null(lbl)){
        lblm[i] <- i
      } else {
        lblm[i] <- lbl[i]
      }
      lblm[i] <- paste(lblm[i]," - ",format(round(ranks.means[i],2),nsmall=2),sep="")
      lbl.size[i] <- nchar(lblm[i])
    }
    lbl.size <- max(lbl.size)
    # Produce plot
    par(mar=c(lbl.size/2, 4, 4, 2) + 0.1)
    plot(1:cols.number,rep(0,times=cols.number),
         xlab="", ylab="", main=title,
         type="p", axes=FALSE, pch=20, lwd=4, col="black", ylim=c(0,k+1), ...)
    if (k>0){
      for (i in 1:k){
        lines(rline[i,],c(i,i), col=pcol[i], lwd = 4)
        lines(rep(rline[i,1],times=2),c(0,i), col="gray", lty = 2)
        lines(rep(rline[i,2],times=2),c(0,i), col="gray", lty = 2)
      }
    }
    axis(1,at=c(1:cols.number),labels=lblm,las=2)
    if (!is.null(select)){
      points(select,0,pch=20,col="red",cex=2)
    }
    # box(which="plot", col="black")
    par(mar=c(5, 4, 4, 2) + 0.1)
  }

  # Line style plot (as in ISF) - vertical
  if(plottype == "vline"){
    # Find groups
    rline <- matrix(NA, nrow=cols.number, ncol=2)
    for (i in 1:cols.number){
      tloc <- which((abs(ranks.means-ranks.means[i])<r.stat[2]) == TRUE)
      rline[i,] <- c(min(tloc),max(tloc))
    }
    # Get rid of duplicates and single member groups
    rline <- unique(rline)
    rline <- rline[apply(rline,1,min) != apply(rline,1,max),]
    # Re-convert to matrix if necessary and find number of remaining groups
    if (length(rline)==2){
      rline <- as.matrix(rline)
      rline <- t(rline)
    }
    k <- nrow(rline)
    # Choose colour depending on Friedman test result
    if (fried.pval > 1-conf.int){pcol <- rep("gray",times=k)} else {
      if (is.null(pcol)){pcol <- rainbow(k)} else {pcol <- rep(pcol,times=k)}
    }
    # Prepare method labels and add mean rank to them
    lbl <- labels
    lblm <- matrix(NA,nrow=1,ncol=cols.number)
    lbl.size <- matrix(NA,nrow=1,ncol=cols.number)
    for (i in 1:cols.number){
      if (is.null(lbl)){
        lblm[i] <- i
      } else {
        lblm[i] <- lbl[i]
      }
      lblm[i] <- paste(lblm[i]," - ",format(round(ranks.means[i],2),nsmall=2),sep="")
      lbl.size[i] <- nchar(lblm[i])
    }
    lbl.size <- max(lbl.size)
    # Produce plot
    par(mar=c(2, lbl.size/2, 4, 2) + 0.1)
    plot(rep(0,times=cols.number),1:cols.number,
         ylab="", xlab="", main=title, 
         type="p", axes=FALSE, pch=20, lwd=4, col="black", xlim=c(0,k+1), ylim=rev(c(1,cols.number)), ...)
    if (k>0){
      for (i in 1:k){
        lines(c(i,i), rline[i,], col=pcol[i], lwd = 4)
        lines(c(0,i), rep(rline[i,1],times=2), col="gray", lty = 2)
        lines(c(0,i), rep(rline[i,2],times=2), col="gray", lty = 2)
      }
    }
    axis(2,at=c(1:cols.number),labels=lblm,las=2)
    if (!is.null(select)){
      points(0,select,pch=20,col="red",cex=2)
    }
    # box(which="plot", col="black")
    par(mar=c(5, 4, 4, 2) + 0.1)
  }


  return(list(means=ranks.means,intervals=ranks.intervals,fpval=fried.pval,cd=r.stat[2]));
}