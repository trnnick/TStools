abc <- function(x,prc=c(0.2,0.3,0.5)){
    # ABC analysis
    #
    # Inputs
    #   x           This can either be an array, where each column is a series, or a vector of values.
    #               If x is an array of time series, then the importance of each series is calculated
    #               as the mean of the observations (sales volume). Otherwise, value can be whatever 
    #               quantity is provided. 
    #   prc         A list of percentages to separate the items in. By default this is c(0.2,0.3,0.5),
    #               but any set of percentage values can be used as long as 0<=prc[i]<=1 and sum(prc)==1.
    #   
    # Outputs
    #   value       A vector containing the importance value of each series.
    #   class       A vector containing the class membership of each series.
    #   rank        A vector containing the rank of each series. Highest value is 1.
    #   conc        The importance concentration of each class, as percentage of total value.
    #
    # Example
    #   x <- abs(matrix(cumsum(rnorm(5400,0,1)),36,150))
    #   abc(x)
    #
    # Nikolaos Kourentzes, 2014 <nikolaos@kourentzes.com>
    # Update 2016.10: Change colour scheme, accept ellipsis, S3method
    # Update 2017.02: Removed plot options and moved to a separate function. 
    
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
    
    names(prc) <- nam.abc 
    
    return(structure(list("type"="ABC","prc"=prc,"value"=x.mean,"class"=x.class,"rank"=x.rank,"conc"=x.imp),class="abc"))
    
}

summary.abc <- function(object,...){
    print(object)
}

print.abc <- function(x,...){
    writeLines(paste(x$type,"analysis"))
    writeLines(paste("        ",colnames(x$conc),"%"))
    row.nams <- paste0(rownames(x$conc),paste0(" (",round(x$prc,2)*100,"%): "))
    row.nams <- paste0(row.nams, as.vector(round(x$conc,3)))
    for (i in 1:length(row.nams)){
        writeLines(row.nams[i])
    }
    
}

plot.abc <- function(x,cex.prc=0.8,...){
    # Plot function for ABC and XYZ
    # Inputs
    #   x           An object of class "abc"
    #   cex.prc     Font size of percentage reported in plot.
    #   ...         Additional arguments can be passed to the plot.
    #
    # Example
    #   x <- abs(matrix(cumsum(rnorm(5400,0,1)),36,150))
    #   z <- abc(x)
    #   plot(z)
    #
    # Nikolaos Kourentzes, 2017 <nikolaos@kourentzes.com>
    
    type <- x$type
    k <- length(x$prc)  
    n <- length(x$value)
    p <- round(n*x$prc)
    nam.abc <- rownames(x$conc)
    x.imp <- x$conc
    x.sort <- x$value[x$rank]
    x.sort <- (x.sort/sum(x.sort))*100
    
    # Get colours
    if (type == "ABC"){
        cmp <- colorRampPalette(brewer.pal(9,"Greens")[4:7])(k)    
    } else if (type == "XYZ"){
        cmp <- colorRampPalette(brewer.pal(9,"Oranges")[4:7])(k)
    }
    
    # Allow user to override plot defaults
    args <- list(...)
    if (!("main" %in% names(args))){
        args$main <- paste0(type," analysis")
    }
    if (!("xlab" %in% names(args))){
        args$xlab <- "Cumulative number of items"
    }
    if (!("ylab" %in% names(args))){
        if (type == "ABC"){
            args$ylab <- "Cumulative importance"
        } else if (type == "XYZ"){
            args$ylab <- "Cumulative error"
        }
    }
    if (!("xaxs" %in% names(args))){
        args$xaxs <- "i"
    }
    if (!("yaxs" %in% names(args))){
        args$yaxs <- "i"
    }
    # Remaining defaults
    args$x <- args$y <- NA
    args$xlim <- args$ylim <- c(0,100)
    # Use do.call to use manipulated ellipsis (...)
    do.call(plot,args)
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