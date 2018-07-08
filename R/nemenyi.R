nemenyi <- function(data, conf.int=0.95, sort=c(TRUE,FALSE), 
                    plottype=c("vline","none","mcb","vmcb","line"), 
                    select=NULL, labels=NULL, ...)
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
    #  select     Highlight selected model (column). Number 1 to k. Select NULL for no highlighting.
    #  labels     Optional labels for models. If NULL column names of 'data' will be used.
    #             If number of labels < k, then it is assumed as NULL. 
    #   ...         Additional arguments can be passed to the plot.
    #  
    # Outputs:
    #   means - mean rank of each method
    #   intervals - Nemenyi intervals for each method
    #   fpval - Friedman test p-value
    #   fH - Friedman hypothesis result
    #   cd - critical distance for Nemenyi
    #   conf.int - the width of the confidence interval
    #   k - number of methods
    #   n - number of measurements
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
    
    # Save the current par() values
    if(plottype!="none"){
        parDefault <- par(no.readonly=TRUE);
        par.mar <- parDefault$mar;
    }
    
    data <- na.exclude(data)
    rows.number <- nrow(data)
    cols.number <- ncol(data)
    
    # Check select < cols.number
    if (!is.null(select)){
        if (select > cols.number){
            select <- NULL
        }
    }
    
    # If plot is asked, always sort the results
    if (plottype != "none"){
        sort <- TRUE
    }
    
    # Checks for labels
    if (is.null(labels)){
        labels <- colnames(data)
        if (is.null(labels)){
            labels <- 1:cols.number
        }
    } else {
        labels <- labels[1:cols.number]
    }
    
    # First run Friedman test. If insignificant then ignore Nemenyi (Weaker)
    fried.pval <- friedman.test(data)$p.value
    if (fried.pval <= 1-conf.int){
        fried.H <- "Ha: Different" # At least one method is different
    } else {
        fried.H <- "H0: Identical" # No evidence of differences between methods
    }
    
    # Nemenyi critical distance and bounds of intervals
    r.stat <- qtukey(conf.int,cols.number,Inf) * sqrt((cols.number*(cols.number+1))/(12*rows.number))
    # 0.5 is needed in order to place the mean ranks in the centre of the interval
    r.stat <- 0.5*c(-r.stat,r.stat)
    
    # Rank methods for each time series
    ranks.matrix <- matrix(NA, nrow=rows.number, ncol=cols.number)
    colnames(ranks.matrix) <- labels
    for (i in 1:rows.number){
        ranks.matrix[i, ] <- rank(data[i,],na.last="keep",ties.method="average")
    }
    
    # Calculate mean rank values
    ranks.means <- colMeans(ranks.matrix)
    
    # Calculate intervals for each of the methods
    ranks.intervals <- rbind(ranks.means + r.stat[1],ranks.means + r.stat[2])
    colnames(ranks.intervals) <- labels
    
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
    
    # Create title for plots
    if(is.null(title)){
        title <- paste("Friedman: ", format(round(fried.pval,3),nsmall=3), " (", fried.H, ") \n Nemenyi CD: ", format(round(r.stat[2],3),nsmall=3), sep = "")
    }
    
    # Labels for plots
    if (is.null(labels)){
        labels <- colnames(ranks.intervals)
    }
    
    # Produce plots
    
    # # For all plots
    if (plottype != "none"){
        args <- list(...)
        args.nms <- names(args)
        # par.nms <- c("ask", "fig", "fin", "lheight", "mai", "mar", "mex", "mfcol", "mfrow", "mfg", "new", "oma", "omd", "omi", "pin", "plt", "ps", "pty", "usr", "xlog", "ylog", "ylbias")
        # loc <- which(args.nms %in% par.nms)
        # par.mar <- args[loc]
        # args <- args[setdiff(1:length(args),loc)]
    }
    
    # MCB style plot
    if(plottype == "mcb"){
        
        cmp <- brewer.pal(3,"Set1")[1:2]
        # Choose colour depending on Friedman test result
        if (fried.pval > 1-conf.int){pcol <- "gray"} else {pcol <- cmp[2]}
        # Find min max
        ymax <- max(ranks.intervals)
        ymin <- min(ranks.intervals)
        ymax <- ymax + 0.1*(ymax-ymin)
        ymin <- ymin - 0.1*(ymax-ymin)
        
        if (!("main" %in% names(args))){
            args$main <- paste0("Friedman: ", format(round(fried.pval,3),nsmall=3), " (", fried.H, ") \n MCB interval: ", format(round(r.stat[2],3),nsmall=3))
        }
        if (!("xlab" %in% names(args))){
            args$xlab <- ""
        }
        if (!("ylab" %in% names(args))){
            args$ylab <- "Mean ranks"
        }
        if (!("xaxs" %in% names(args))){
            args$xaxs <- "i"
        }
        if (!("yaxs" %in% names(args))){
            args$yaxs <- "i"
        }
        # Remaining defaults
        args$x <- args$y <- NA
        if(is.null(args$xlim)){
            args$xlim <- c(0,cols.number+1)
        }
        if(is.null(args$ylim)){
            args$ylim <- c(ymin,ymax)
        }
        args$axes <- FALSE
        
        if(all(par.mar==(c(5,4,4,2)+0.1))){
            par.mar <- c(2, 2, 0, 0) + 0.1
        }
        if(!is.null(args$main)){
            if(args$main!=""){
                par.mar <- par.mar + c(0,0,4,0)
            }
        }
        par(mar=par.mar)

        # Use do.call to use manipulated ellipsis (...)
        do.call(plot,args)
        # Plot rest
        points(1:cols.number,ranks.means,pch=20,lwd=4)
        axis(1,at=c(1:cols.number),labels=labels)
        axis(2)
        box(which="plot", col="black")
        # Intervals for best method
        if (is.null(select)){
            select <- 1
        }
        lines(c(0,cols.number+1),rep(ranks.intervals[1,select],times=2), col="gray", lty=2)
        lines(c(0,cols.number+1),rep(ranks.intervals[2,select],times=2), col="gray", lty=2)    
        # Intervals for all methods
        for (i in 1:cols.number){
            lines(rep(i,times=2),ranks.intervals[,i], type="b", lwd=2, col=pcol);
        }
        # Highlight identical
        idx <- !((ranks.intervals[2,select] < ranks.intervals[1,]) | (ranks.intervals[1,select] > ranks.intervals[2,]))
        points((1:cols.number)[idx],ranks.means[idx],pch=20,lwd=4,col=cmp[1])
    }
    
    # MCB style plot - vertical
    if(plottype == "vmcb"){
        
        cmp <- brewer.pal(3,"Set1")[1:2]
        # Find max label size
        lbl.size <- max(nchar(labels))
        # Choose colour depending on Friedman test result
        if (fried.pval > 1-conf.int){pcol <- "gray"} else {pcol <- cmp[2]}
        # Find min max
        xmax <- max(ranks.intervals)
        xmin <- min(ranks.intervals)
        xmax <- xmax + 0.1*(xmax-xmin)
        xmin <- xmin - 0.1*(xmax-xmin)
        
        if (!("main" %in% names(args))){
            args$main <- paste0("Friedman: ", format(round(fried.pval,3),nsmall=3), " (", fried.H, ") \n MCB interval: ", format(round(r.stat[2],3),nsmall=3))
        }
        if (!("ylab" %in% names(args))){
            args$ylab <- ""
        }
        if (!("xlab" %in% names(args))){
            args$xlab <- "Mean ranks"
        }
        if (!("xaxs" %in% names(args))){
            args$xaxs <- "i"
        }
        if (!("yaxs" %in% names(args))){
            args$yaxs <- "i"
        }
        # Remaining defaults
        args$x <- args$y <- NA
        if(is.null(args$xlim)){
            args$xlim <- c(xmin,xmax)
        }
        if(is.null(args$ylim)){
            args$ylim <- c(0,cols.number+1)
        }
        args$axes <- FALSE
        # Use do.call to use manipulated ellipsis (...)
        # if (!("mar" %in% names(par.mar))){
            # par.mar <- c(2, lbl.size/1.5, 4, 2) + 0.1
        # }
        # If the default mar is used, do things.
        if(all(par.mar==(c(5,4,4,2)+0.1))){
            par.mar <- c(2, lbl.size/2, 0, 0) + 0.1
        }
        else{
            par.mar <- par.mar + c(2, lbl.size/2, 0, 0)
        }
        if(!is.null(args$main)){
            if(args$main!=""){
                par.mar <- par.mar + c(0,0,4,0)
            }
        }
        par(mar=par.mar)
        
        do.call(plot,args)
        # Plot rest
        points(ranks.means,1:cols.number,pch=20,lwd=4)
        axis(2,at=c(1:cols.number),labels=labels,las=2)
        axis(1)
        box(which="plot", col="black")
        # Intervals for best method
        if (is.null(select)){
            select <- 1
        }
        lines(rep(ranks.intervals[1,select],times=2), c(0,cols.number+1), col="gray", lty=2)
        lines(rep(ranks.intervals[2,select],times=2), c(0,cols.number+1), col="gray", lty=2)    
        # Intervals for all methods
        for (i in 1:cols.number){
            lines(ranks.intervals[,i], rep(i,times=2), type="b", lwd=2, col=pcol);
        }
        # Highlight identical
        idx <- !((ranks.intervals[2,select] < ranks.intervals[1,]) | (ranks.intervals[1,select] > ranks.intervals[2,]))
        points(ranks.means[idx],(1:cols.number)[idx],pch=20,lwd=4,col=cmp[1])
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
        cmp <- colorRampPalette(brewer.pal(12,"Paired"))(k)
        if (fried.pval > 1-conf.int){pcol <- rep("gray",times=k)} else {pcol <- cmp}
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
        if (!("main" %in% names(args))){
            args$main <- paste0("Friedman: ", format(round(fried.pval,3),nsmall=3), " (", fried.H, ") \n CD: ", format(round(r.stat[2],3),nsmall=3))
        }
        if (!("ylab" %in% names(args))){
            args$ylab <- ""
        }
        if (!("xlab" %in% names(args))){
            args$xlab <- ""
        }
        if (!("xaxs" %in% names(args))){
            args$xaxs <- "i"
        }
        if (!("yaxs" %in% names(args))){
            args$yaxs <- "i"
        }
        # Remaining defaults
        args$x <- args$y <- NA
        if(is.null(args$xlim)){
            args$xlim <- c(1,cols.number)
        }
        if(is.null(args$ylim)){
            args$ylim <- c(0,k+1)
        }
        args$axes <- FALSE
        # Use do.call to use manipulated ellipsis (...)
        # temp.mar <- par()$mar
        # if (!("mar" %in% names(par.mar))){
        #   par.mar$mar <- c(lbl.size/2, 4, 4, 2) + 0.1
        # }
        
        # If the default mar is used, do things.
        if(all(par.mar==(c(5,4,4,2)+0.1))){
            par.mar <- c(lbl.size/2, 2, 0, 2) + 0.1
        }
        else{
            par.mar <- par.mar + c(lbl.size/2, 0, 0, 0)
        }
        if(!is.null(args$main)){
            if(args$main!=""){
                par.mar <- par.mar + c(0,0,4,0)
            }
        }
        par(mar=par.mar)
        
        do.call(plot,args)
        points(1:cols.number,rep(0,cols.number),pch=20,lwd=4)
        if (k>0){
            for (i in 1:k){
                lines(rline[i,],c(i,i), col=pcol[i], lwd = 4)
                lines(rep(rline[i,1],times=2),c(0,i), col="gray", lty = 2)
                lines(rep(rline[i,2],times=2),c(0,i), col="gray", lty = 2)
            }
        }
        axis(1,at=c(1:cols.number),labels=lblm,las=2)
        if (!is.null(select)){
            points(select,0,pch=20,col=brewer.pal(3,"Set1")[1],cex=2)
        }
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
        cmp <- colorRampPalette(brewer.pal(12,"Paired"))(k)
        if (fried.pval > 1-conf.int){pcol <- rep("gray",times=k)} else {pcol <- cmp}
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
        if (!("main" %in% names(args))){
            args$main <- paste0("Friedman: ", format(round(fried.pval,3),nsmall=3), " (", fried.H, ") \n CD: ", format(round(r.stat[2],3),nsmall=3))
        }
        if (!("ylab" %in% names(args))){
            args$ylab <- ""
        }
        if (!("xlab" %in% names(args))){
            args$xlab <- ""
        }
        if (!("xaxs" %in% names(args))){
            args$xaxs <- "i"
        }
        if (!("yaxs" %in% names(args))){
            args$yaxs <- "i"
        }
        # Remaining defaults
        args$x <- args$y <- NA
        if(is.null(args$xlim)){
            args$xlim <- c(0,k+1)
        }
        if(is.null(args$ylim)){
            args$ylim <- c(1,cols.number)
        }
        args$axes <- FALSE
        # Use do.call to use manipulated ellipsis (...)
        # temp.mar <- par()$mar
        # if (!("mar" %in% names(par.mar))){
            # par.mar$mar <- c(2, lbl.size/2, 4, 2) + 0.1
        # }
        # If the default mar is used, do things.
        if(all(par.mar==(c(5,4,4,2)+0.1))){
            par.mar <- c(2, lbl.size/2, 2, 0) + 0.1
        }
        else{
            par.mar <- par.mar + c(0, lbl.size/2, 0, 0)
        }
        if(!is.null(args$main)){
            if(args$main!=""){
                par.mar <- par.mar + c(0,0,4,0)
            }
        }
        par(mar=par.mar)
        
        do.call(plot,args)
        points(rep(0,cols.number),1:cols.number,pch=20,lwd=4)
        if (k>0){
            for (i in 1:k){
                lines(c(i,i), rline[i,], col=pcol[i], lwd = 4)
                lines(c(0,i), rep(rline[i,1],times=2), col="gray", lty = 2)
                lines(c(0,i), rep(rline[i,2],times=2), col="gray", lty = 2)
            }
        }
        axis(2,at=c(1:cols.number),labels=lblm,las=2)
        if (!is.null(select)){
            points(0,select,pch=20,col=brewer.pal(3,"Set1")[1],cex=2)
        }
    }
    
    # Revert to the original par() values
    if(plottype!="none"){
        par(parDefault)
    }
    
    return(structure(list("means"=ranks.means,"intervals"=ranks.intervals,"fpval"=fried.pval,"fH"=fried.H,"cd"=r.stat[2],"conf.int"=conf.int,"k"=cols.number,"n"=rows.number),class="nemenyi"))
    
}

summary.nemenyi <- function(object,...){
    print(object)
}

print.nemenyi <- function(x,...){
    
    writeLines("Friedman and Nemenyi Tests")
    writeLines(paste0("The significance level is ", (1-x$conf.int)*100, "%"))
    writeLines(paste0("Number of observations is ", x$n, " and number of methods is ", x$k))
    writeLines(paste0("Friedman test p-value: ", format(round(x$fpval,4),nsmall=4) , " - ", x$fH))
    writeLines(paste0("Nemenyi critical distance: ", format(round(x$cd,4),nsmall=4)))
    
}