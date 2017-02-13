elm <- function(y,hd=NULL,type=c("lasso","step","lm"),reps=20,comb=c("median","mean","mode"),
                lags=NULL,difforder=NULL,outplot=c(FALSE,TRUE),sel.lag=c(TRUE,FALSE),direct=c(FALSE,TRUE),
                allow.det.season=c(TRUE,FALSE),det.type=c("auto","bin","trg"),
                xreg=NULL,xreg.lags=NULL,barebone=c(FALSE,TRUE)){
    
    # Defaults
    type <- match.arg(type,c("lasso","step","lm"))
    comb <- match.arg(comb,c("median","mean","mode"))
    outplot <- outplot[1]
    sel.lag <- sel.lag[1]
    direct <- direct[1]
    allow.det.season <- allow.det.season[1]
    det.type <- det.type[1]
    barebone <- barebone[1]

    # Check if y input is a time series
    if (!(any(class(y) == "ts") | any(class(y) == "msts"))){
      stop("Input y must be of class ts or msts.")
    }
    
    # Check xreg inputs
    if (!is.null(xreg)){
      x.n <- length(xreg[1,])
      if (!is.null(xreg.lags)){
        if (length(xreg.lags) != x.n){
          stop("Argument xreg.lags must be a list with as many elements as xreg variables (columns).")
        }
      }
    }
    
    # Get frequency
    ff.ls <- get.ff(y)    
    ff <- ff.ls$ff
    ff.n <- ff.ls$ff.n
    rm("ff.ls")
    
    # Default lagvector
    xreg.ls <- def.lags(lags,ff,xreg.lags,xreg)
    lags <- xreg.ls$lags
    xreg.lags <- xreg.ls$xreg.lags
    rm("xreg.ls")
    
    # Pre-process data (same for MLP and ELM)
    PP <- preprocess(y,m,lags,difforder,sel.lag,allow.det.season,det.type,ff,ff.n,xreg,xreg.lags)
    Y <- PP$Y
    X <- PP$X
    sdummy <- PP$sdummy
    difforder <- PP$difforder
    det.type <- PP$det.type
    lags <- PP$lags
    xreg.lags <- PP$xreg.lags
    sc <- PP$sc
    xreg.minmax <- PP$xreg.minmax
    d <- PP$d
    y.d <- PP$y.d
    y.ud <- PP$y.ud
    frm <- PP$frm
    ff.det <- PP$ff.det
    lag.max <- PP$lag.max
    rm("PP")
    
    if (is.null(hd)){
      hd <- min(100-60*(type=="step" | type=="lm"),max(4,length(Y)-2-as.numeric(direct)*length(X[1,])))
    }

    # Train ELM
    # If single hidden layer switch to fast, unless requested otherwise
    if (length(hd)==1 & barebone == TRUE){
        # Switch to elm.fast
        f.elm <- elm.fast(Y,X,hd=hd,reps=reps,comb=comb,type=type,direct=direct,linscale=FALSE,output="linear",core=TRUE)
        # Post-process output
        Yhat <- f.elm$fitted.all
        for (r in 1:reps){
            # Reverse scaling
            yhat <- linscale(Yhat[,r],sc$minmax,rev=TRUE)$x
            # Undifference - this is 1-step ahead undifferencing
            y.ud[[d+1]] <- yhat
            if (d>0){
                for (i in 1:d){
                    n.ud <- length(y.ud[[d+2-i]])
                    n.d <- length(y.d[[d+1-i]])
                    y.ud[[d+1-i]] <- y.d[[d+1-i]][(n.d-n.ud-difforder[d+1-i]+1):(n.d-difforder[d+1-i])] + y.ud[[d+2-i]]
                }
            }
            Yhat[,r] <- head(y.ud,1)[[1]]
        }
        
    } else {
        # Rely on neuralnet, very slow when number of inputs is large
        net <- neuralnet(frm,cbind(Y,X),hidden=hd,threshold=10^10,rep=reps,err.fct="sse",linear.output=FALSE)
        
        # Get weights for each repetition
        W <- W.dct <- vector("list",reps)
        B <- vector("numeric",reps)
        Yhat <- array(NA,c((length(y)-sum(difforder)-lag.max),reps))
        x.names <- colnames(X)
        
        for (r in 1:reps){
            H <- as.matrix(tail(compute(net,X,r)$neurons,1)[[1]][,2:(tail(hd,1)+1)])
            if (direct==TRUE){
                Z <- cbind(H,X)
            } else {
                Z <- H
            }
            w.out <- elm.train(Y,Z,type,X,direct,hd,output="linear")
            B[r] <- w.out[1]                                  # Bias (Constant)
            if (direct == TRUE){                              # Direct connections
                w.dct <- w.out[(1+hd+1):(1+hd+dim(X)[2]),,drop=FALSE]
                if (!is.null(x.names)){
                    rownames(w.dct) <- x.names
                }
                W.dct[[r]] <- w.dct
            }
            W[[r]] <- w.out[2:(1+hd),,drop=FALSE]             # Hidden layer
            
            # Produce fit
            yhat.sc <- H %*% W[[r]] + B[r] + if(direct!=TRUE){0}else{X %*% W.dct[[r]]}
            
            # Post-process
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
        
    } # Close elm.fast if
    
    # Combine forecasts
    yout <- frc.comb(Yhat,comb)
    
    # Convert to time series
    yout <- ts(yout,end=end(y),frequency=frequency(y))
    
    MSE <- mean((y[(lag.max+1+sum(difforder)):length(y)] - yout)^2)
    
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
        lines(yout,col="blue")
        
    }
    
    if (barebone == FALSE){
        W.in <- NULL
        class.type <- "elm"
    } else {
        net <- NULL
        hd <- f.elm$hd
        W <- f.elm$W
        B <- f.elm$b
        W.in <- f.elm$W.in
        W.dct <- f.elm$W.dct
        class.type <- c("elm","elm.fast")
    }
    
    return(structure(list("net"=net,"hd"=hd,"W.in"=W.in,"W"=W,"b"=B,"W.dct"=W.dct,
                          "lags"=lags,"xreg.lags"=xreg.lags,"difforder"=difforder,
                          "sdummy"=sdummy,"ff.det"=ff.det,"det.type"=det.type,"y"=y,"minmax"=sc$minmax,"xreg.minmax"=xreg.minmax,
                          "comb"=comb,"type"=type,"direct"=direct,"fitted"=yout,"MSE"=MSE),class=class.type))
    
}

forecast.elm <- function(object,h=NULL,y=NULL,xreg=NULL,...){
  forecast.net(object,h=h,y=y,xreg=xreg,...)
}

elm.thief <- function(y,h=NULL,...){
    # This is a wrapper function to use MLP with THieF
    
    # Remove level input from ellipsis
    ellipsis.args <- list(...)
    ellipsis.args$level <- NULL
    ellipsis.args$y <- y
    
    # Fit network
    fit <- do.call(elm,ellipsis.args) 
    
    # Default h
    if (is.null(h)){
        h <- frequency(y)
    }
    # Check if xreg was given and pass to forecast
    if ("xreg" %in% names(ellipsis.args)){
        xreg <- ellipsis.args$xreg
    } else {
        xreg <- NULL
    }
    
    # Forecast
    out <- forecast(fit,h,xreg)
    # Make fitted values span the complete sample
    n <- length(out$x)
    m <- length(out$fitted)
    if (m < n){
        out$fitted <- ts(c(rep(NA,n-m),out$fitted),frequency=frequency(out$fitted),end=end(out$fitted))
    }
    
    return(out)
    
}

plot.elm <- function(x, r=1, ...){
    plot.net(x,r)
}

print.elm <- function(x, ...){
    print.net(x,...)
}

elm.train <- function(Y,Z,type,X,direct,hd,output="linear"){
# Find output weights for ELM
  switch(type,
         "lasso" = {
             if (output=="logistic"){
                fit <- suppressWarnings(cv.glmnet(Z,cbind(Y),family="binomial"))
             } else {
                fit <- suppressWarnings(cv.glmnet(Z,cbind(Y)))
             }
             cf <- as.vector(coef(fit))
         },
         {
           reg.data <- as.data.frame(cbind(Y,Z))
           colnames(reg.data) <- c("Y",paste0("X",1:(tail(hd,1)+as.numeric(direct)*length(X[1,]))))
           # Take care of linear dependency
           alias.fit <- alias(as.formula(paste0("Y~",paste0("X",1:(tail(hd,1)+as.numeric(direct)*length(X[1,])),collapse="+"))),data=reg.data)
           alias.x <- rownames(alias.fit$Complete)
           frm <- as.formula(paste0("Y~",paste0(setdiff(colnames(reg.data)[2:(hd+1+as.numeric(direct)*length(X[1,]))],alias.x),collapse="+")))
           fit <- suppressWarnings(lm(frm,reg.data))
           if (type == "step"){
             fit <- suppressWarnings(stepAIC(fit,trace=0)) # ,direction="backward",k=log(length(Y)))) # BIC criterion
           }
           cf.temp <- coef(fit)
           loc <- which(colnames(reg.data) %in% names(cf.temp))
           cf <- rep(0,(tail(hd,1)+1+direct*length(X[1,])))
           cf[1] <- cf.temp[1]
           cf[loc] <- cf.temp[2:length(cf.temp)]
         })
  return(cbind(cf))
}