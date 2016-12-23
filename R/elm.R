elm <- function(y,hd=50,type=c("lasso","step","lm"),reps=20,comb=c("median","mean","mode"),
                lags=NULL,difforder=-1,outplot=c(FALSE,TRUE)){
    
    # Defaults
    type <- type[1]
    comb <- comb[1]
    outplot <- outplot[1]
    
    # Check if y input is a time series
    if (class(y) != "ts"){
        stop("Input y must be of class ts.")
    }
    
    # Default lagvector
    if (is.null(lags)){
        lags <- 1:frequency(y)
    }
    
    # Find differencing order
    if (difforder == -1){
        # Identify difforder automatically
        st <- seasplot(y)
        difforder <- NULL
        if (st$trend.exist == TRUE){
            difforder <- 1
        }
        if (st$season.exist == TRUE){
            difforder <- c(difforder,frequency(y))
        }
    }
    
    # Apply differencing
    d <- length(difforder)
    y.d <- y.ud <- vector("list",d+1)
    y.d[[1]] <- y
    names(y.d)[1] <- "d0"
    if (d>0){
        for (i in 1:d){
            y.d[[i+1]] <- diff(y.d[[i]],difforder[i])
            names(y.d)[i+1] <- paste0(names(y.d)[i],"d",difforder[i])
        }
    }
    names(y.ud) <- names(y.d)
    
    # Scale target
    sc <- linscale(tail(y.d,1)[[1]],minmax=list("mn"=-.8,"mx"=0.8))
    y.sc <- sc$x
    n <- length(y.sc)
    
    # Prepare inputs & target
    y.sc.lag <- lagmatrix(y.sc,unique(c(0,lags)))
    Y <- y.sc.lag[(max(lags)+1):n,1,drop=FALSE]
    X <- y.sc.lag[(max(lags)+1):n,2:(length(lags)+1),drop=FALSE]
    colnames(X) <- paste0("X",lags)
    
    # Create network
    frm <- paste0(colnames(X),"+",collapse="")
    frm <- substr(frm,1,nchar(frm)-1)
    frm <- as.formula(paste0("Y~",frm))
    net <- neuralnet(frm,cbind(Y,X),hidden=hd,threshold=10^10,rep=reps,err.fct="sse",linear.output=FALSE)
    
    # Get hidden nodes output and weights for each repetition
    H <- W <- vector("list",reps)
    Yhat <- array(NA,c((length(y)-sum(difforder)-max(lags)),reps))
    
    for (r in 1:reps){
        H[[r]] <- as.matrix(tail(compute(net,X,i)$neurons,1)[[1]][,2:(tail(hd,1)+1)])
    
        # Calculate regression
        switch(type,
               "lasso" = {
                   fit <- cv.glmnet(H[[r]],cbind(Y))
                   cf <- as.vector(coef(fit))
               },
               {
                   reg.data <- as.data.frame(cbind(Y,H[[r]]))
                   colnames(reg.data) <- c("Y",paste0("X",1:tail(hd,1)))
                   fit <- lm(Y~.,reg.data)
                   cf <- coef(fit)
                   if (type == "step"){
                       fit <- stepAIC(fit,trace=0)
                       cf.temp <- coef(fit)
                       loc <- which(colnames(reg.data) %in% names(cf.temp))
                       cf <- rep(0,tail(hd,1)+1)
                       cf[1] <- cf.temp[1]
                       cf[loc] <- cf.temp[2:length(cf.temp)]
                   }
               })
        W[[r]] <- cbind(cf)
    
        # Produce fit
        yhat.sc <- cbind(1,H[[r]]) %*% W[[r]]
        yhat <- linscale(yhat.sc,sc$minmax,rev=TRUE)$x
    
        # # Check unscaled, but differenced fit
        # plot(1:length(tail(y.d,1)[[1]]),tail(y.d,1)[[1]], type="l")
        # lines((max(lags)+1):length(tail(y.d,1)[[1]]),yhat,col="red")
    
        # Undifference - this is 1-step ahead undifferencing
        y.ud[[d+1]] <- yhat
        if (d>0){
            for (i in 1:d){
                n.ud <- length(y.ud[[d+2-i]])
                n.d <- length(y.d[[d+1-i]])
                y.ud[[d+1-i]] <- y.d[[d+1-i]][(n.d-n.ud-difforder[d+1-i]+1):(n.d-difforder[d+1-i])] + y.ud[[d+2-i]]
            }
        }
        yout <- head(y.ud,1)[[1]]
        
        # yout <- ts(yout,end=end(y),frequency=frequency(y))
        # # Check undifferences and unscaled
        # plot(y)
        # lines(yout,col="red")
        # # plot(1:length(y),y,type="l")
        # # lines((max(lags)+1+sum(difforder)):length(y),y.ud[[1]],col="red")
    
        Yhat[,r] <- yout
        
    } # Close reps
        
    # If reps>1 combine forecasts
    if (reps>1){
        switch(comb,
               "median" = {yout <- apply(Yhat,1,median)},
               "mean" = {yout <- apply(Yhat,1,mean)},
               "mode" = {yout <- sapply(apply(Yhat,1,kdemode),function(x){x[[1]][1]})}
               )
    } else {
        yout <- Yhat[,1]
    }
    # Convert to time series
    yout <- ts(yout,end=end(y),frequency=frequency(y))
    
    MSE <- mean((y[(max(lags)+1+sum(difforder)):length(y)] - yout)^2)
    
    # Construct plot
    if (outplot==TRUE){
        
        plot(y)
        if (reps>1){
            for (i in 1:reps){
                temp <- Yhat[,i]
                temp <- ts(temp,frequency=frequency(y),end=end(y))
                lines(temp,col="grey")
            }
        }
        lines(yout,col="red")
        
    }
    
    return(structure(list("net"=net,"hd"=hd,"W"=W,"lags"=lags,"difforder"=difforder,"y"=y,"minmax"=sc$minmax,"comb"=comb,"type"=type,"fitted"=yout,"MSE"=MSE),class="elm"))
    
}

forecast.elm <- function(fit,h=NULL,outplot=c(FALSE,TRUE),y=NULL,...){
# Produce forecast with ELM
    
    outplot <- outplot[1]
    
    if (is.null(y)){
        y <- fit$y
    }
    
    if (is.null(h)){
        h <- frequency(y)
    }
    
    # Get stuff from fit list
    net <- fit$net
    hd <- fit$hd
    W <- fit$W
    lags <- fit$lags
    difforder <- fit$difforder
    minmax <- fit$minmax
    comb <- fit$comb
    fitted <- fit$fitted
    reps <- length(net$weights)
    
    # Apply differencing
    d <- length(difforder)
    y.d <- y.ud <- vector("list",d+1)
    y.d[[1]] <- y
    names(y.d)[1] <- "d0"
    if (d>0){
        for (i in 1:d){
            y.d[[i+1]] <- diff(y.d[[i]],difforder[i])
            names(y.d)[i+1] <- paste0(names(y.d)[i],"d",difforder[i])
        }
    }
    Y <- as.vector(linscale(y.d[[d+1]],minmax=minmax)$x)

    Yfrc <- array(NA,c(h,reps))
    
    # For each repetition
    for (r in 1:reps){   
        
        frc.sc <- vector("numeric",h)
        for (i in 1:h){
            if (i == 1){
                temp <- NULL
            } else {
                temp <- frc.sc[1:(i-1)]
            }
            xi <- rev(tail(c(Y,temp),max(lags)))
            xi <- rbind(xi[lags])
            H <- t(as.matrix(tail(compute(net,xi,r)$neurons,1)[[1]][,2:(tail(hd,1)+1)]))
            yhat.sc <- cbind(1,H) %*% W[[r]]
            frc.sc[i] <- yhat.sc
        }
        
        # Reverse scaling
        frc <- linscale(frc.sc,minmax,rev=TRUE)$x
        
        # Reverse differencing
        f.ud <- vector("list",d+1)
        names(f.ud) <- names(y.ud)
        f.ud[[d+1]] <- frc
        if (d>0){
            for (i in 1:d){
                temp <- c(tail(y.d[[d+1-i]],difforder[d+1-i]),f.ud[[d+2-i]])
                n.t <- length(temp)
                for (j in 1:(n.t-difforder[d+1-i])){
                    temp[difforder[d+1-i]+j] <- temp[j] + temp[difforder[d+1-i]+j]
                }
                f.ud[[d+1-i]] <- temp[(difforder[d+1-i]+1):n.t]
            }
        }
        fout <- head(f.ud,1)[[1]]    
        
        Yfrc[,r] <- fout
    }
    
    # If reps>1 combine forecasts
    if (reps>1){
        switch(comb,
               "median" = {fout <- apply(Yfrc,1,median)},
               "mean" = {fout <- apply(Yfrc,1,mean)},
               "mode" = {fout <- sapply(apply(Yfrc,1,kdemode),function(x){x[[1]][1]})}
        )
    } else {
        fout <- Yfrc[,1]
    }
    
    fout <- ts(fout,frequency=frequency(y),start=c(end(y)[1],end(y)[2]+1))
    
    if (outplot==TRUE){
        ts.plot(y,fitted,fout,col=c("black","red","red"))
        if (reps>1){
            for (r in 1:reps){
                temp <- Yfrc[,r]
                temp <- ts(temp,frequency=frequency(fout),end=end(fout))
                lines(temp,col="grey")
            }
            lines(fout,col="red")
        }
        
    }
    
    return(fout)
    
}

print.elm <- function(x, ...){
    
    difforder <- x$difforder
    d <- length(difforder)
    reps <- length(x$net$weights)
    hd <- x$hd
    if (length(hd)>1){
        hdt <- paste0(hd,",",collapse="")
        hdt <- paste0("(", substr(hdt,1,nchar(hdt)-1) ,")")
        hde <- "s"
    } else {
        hdt <- hd
        if (hd>1){
            hde <- "s"
        } else {
            hde <- ""
        }
    }
        
    writeLines(paste0("ELM fit with ", hdt," hidden node",hde," and ", reps, " repetition",if(reps>1){"s"},"."))
    if (d>0){
        writeLines(paste0("Series modelled in differences: ", paste0("D",difforder,collapse=""), "."))
    }
    if (reps>1){
        writeLines(paste0("Forecast combined using the ", x$comb, " operator."))
    }
    writeLines(paste0("Output weight estimation using: ", x$type, "."))
    writeLines(paste0("MSE: ",round(x$MSE,4),"."))
    
}