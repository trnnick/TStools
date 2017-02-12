mlp <- function(y,m=frequency(y),hd=NULL,reps=20,comb=c("median","mean","mode"),
                lags=NULL,difforder=NULL,outplot=c(FALSE,TRUE),sel.lag=c(TRUE,FALSE),
                allow.det.season=c(TRUE,FALSE),det.type=c("auto","bin","trg"),
                xreg=NULL, xreg.lags=NULL,hd.auto.type=c("set","valid","cv","elm"),
                hd.max=NULL, ...){
  
    # hd.max is only relevant to valid and cv
    
    # Defaults
    comb <- comb[1]
    outplot <- outplot[1]
    sel.lag <- sel.lag[1]
    allow.det.season <- allow.det.season[1]
    det.type <- det.type[1]
    hd.auto.type <- hd.auto.type[1]
    
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
    
    # Auto specify number of hidden nodes
    mseH <- NULL
    if (is.null(hd)){
      switch(hd.auto.type,
             "set" = {hd <- 5},
             "valid" = {hd <- auto.hd.cv(Y,X,frm,comb,reps,type="valid",hd.max)
                        mseH <- hd$mseH; hd <- hd$hd},
             "cv" = {hd <- auto.hd.cv(Y,X,frm,comb,reps,type="cv",hd.max)
                     mseH <- hd$mseH; hd <- hd$hd},
             "elm" = {hd <- auto.hd.elm(Y,X,frm)}
             )
    }
    
    # Create network
    net <- neuralnet(frm,cbind(Y,X),hidden=hd,rep=reps,err.fct="sse",linear.output=TRUE,...)
    # In case some networks did not train reduce the number of available repetitions
    reps <- length(net$weights)

    # Produce forecasts
    Yhat <- array(NA,c((length(y)-sum(difforder)-lag.max),reps))
    
    for (r in 1:reps){
        
        yhat.sc <- compute(net,X,r)$net.result
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
    
    return(structure(list("net"=net,"hd"=hd,"lags"=lags,"xreg.lags"=xreg.lags,"difforder"=difforder,"sdummy"=sdummy,"ff.det"=ff.det,
                          "det.type"=det.type,"y"=y,"minmax"=sc$minmax,"xreg.minmax"=xreg.minmax,"comb"=comb,"fitted"=yout,
                          "MSE"=MSE,"MSEH"=mseH),class="mlp"))
    
}

forecast.net <- function(object,h=NULL,y=NULL,xreg=NULL,...){ 
    # Produce forecast with NNs
    
    if (is.null(y)){
        y <- object$y
    }

    # Get frequency
    ff.ls <- get.ff(y)    
    ff <- ff.ls$ff
    ff.n <- ff.ls$ff.n
    rm("ff.ls")
    
    if (is.null(h)){
        h <- max(ff)
    }
    
    # Get stuff from object list
    cl.object <- class(object)
    is.elm.fast <- any(cl.object == "elm.fast")
    if (!is.elm.fast){
        net <- object$net
        reps <- length(net$weights)
        if (class(object) == "elm"){
            W <- object$W
            direct <- object$direct
        } else {
            W <- NULL
            direct <- FALSE
        }
    } else {
        # Get elm.fast definition
        direct <- object$direct
        W.in <- object$W.in
        W <- object$W
        W.dct <- object$W.dct
        b <- object$b
        reps <- length(b)
    }
    hd <- object$hd
    lags <- object$lags
    xreg.lags <- object$xreg.lags
    difforder <- object$difforder
    sdummy <- object$sdummy
    det.type <- object$det.type
    minmax <- object$minmax
    xreg.minmax <- object$xreg.minmax
    comb <- object$comb
    fitted <- object$fitted
    ff.det <- object$ff.det
    ff.n.det <- length(ff.det)

    # Temporal aggregation can mess-up start/end of ts, so lets fix it
    fstart <- c(end(y)[1],end(y)[2]+1)
    if (is.na(fstart[2])){  # If the second element of end(y) does not exist because it is fractional
        fstart <- end(y) + deltat(y)
    }
    
    # Check xreg inputs
    if (!is.null(xreg)){
        x.n <- length(xreg[1,])
        if (length(xreg.lags) != x.n){
            stop("Number of xreg inputs is not consistent with network specification (number of xreg.lags).")
        }
        if (length(xreg[,1]) < length(y)+h){
            stop("Length of xreg must be longer that y + forecast horizon.")
        }
    } else {
        x.n <- 0
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
    Y <- as.vector(linscale(y.d[[d+1]],minmax=minmax)$x)
    
    # Scale xreg
    if (x.n > 0){
        xreg.sc <- xreg
        for (i in 1:x.n){
            xreg.sc[,i] <- linscale(xreg[,i],minmax=xreg.minmax[[i]])$x
        }
        # Starting point of xreg
        xstart <- length(y)+1
    } 
    
    if (sdummy == TRUE){
        temp <- ts(1:h,start=fstart,frequency=max(ff.det))
        Xd <- vector("list",ff.n.det)
        
        for (s in 1:ff.n.det){
            Xd[[s]] <- seasdummy(h,m=ff.det[s],y=temp,type=det.type)
            colnames(Xd[[s]]) <- paste0("D",s,".",1:length(Xd[[s]][1,]))
            if (det.type=="trg"){
                Xd[[s]] <- Xd[[s]][,1:2]
            }
        }
        Xd <- do.call(cbind,Xd)
        # Xd <- seasdummy(h,y=temp,type=det.type)
    }
    
    Yfrc <- array(NA,c(h,reps),dimnames=list(paste0("t+",1:h),paste0("NN.",1:reps)))
    if (length(lags)>0){
        ylag <- max(lags)
    } else {
        ylag <- 0
    }
    
    # For each repetition
    for (r in 1:reps){   
        
        frc.sc <- vector("numeric",h)
        for (i in 1:h){
            
            # Construct inputs
            if (i == 1){
                temp <- NULL
            } else {
                temp <- frc.sc[1:(i-1)]
            }
            xi <- rev(tail(c(Y,temp),ylag)) # Reverse for lags
            xi <- xi[lags]
            # Construct xreg inputs
            if (x.n > 0){
                Xreg <- vector("list",x.n)
                for (j in 1:x.n){
                    if (length(xreg.lags[[j]])>0){
                        xreg.temp <- xreg.sc[(xstart+i-1):(xstart-max(xreg.lags[[j]])+i-1),j] # Reversing is happening in the indices
                        Xreg[[j]] <- xreg.temp[xreg.lags[[j]]+1]
                    }
                }
                Xreg.all <- unlist(Xreg)
                xi <- c(xi,Xreg.all)
            }
            xi <- rbind(xi)
            # Construct seasonal dummies inputs
            if (sdummy == TRUE){
                xi <- cbind(xi,Xd[i,,drop=FALSE])
            }
            
            # Calculate forecasts
            if (any(class(object) == "mlp")){
                yhat.sc <- compute(net,xi,r)$net.result  
            } else {
                if (is.elm.fast){
                    yhat.sc <- predict.elm.fast.internal(xi,W.in[[r]],W[[r]],b[r],W.dct[[r]],direct)
                } else {
                    H <- t(as.matrix(tail(compute(net,xi,r)$neurons,1)[[1]][,2:(tail(hd,1)+1)]))
                    if (direct == TRUE){
                        Z <- cbind(H,xi)
                    } else {
                        Z <- H
                    }
                    yhat.sc <- cbind(1,Z) %*% W[[r]]  
                }
                
            }
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
    
    # Combine forecasts
    fout <- frc.comb(Yfrc,comb)
    
    fout <- ts(fout,frequency=frequency(y),start=fstart)
    
    # Prepare output
    out <- list("method"=class(object),"mean"=fout,
                "all.mean"=ts(Yfrc,frequency=frequency(y),start=fstart),
                "x"=y,"fitted"=fitted,"residuals"=y-fitted)
    return(structure(out,class=c("forecast.net","forecast")))
    
}

plot.forecast.net <- function(x,...){
    # Plot function for NNs
    method <- x$method
    if (any(method=="elm.fast")){
        method <- "ELM"
    }
    reps <- dim(x$all.mean)[2]
    ts.plot(x$x,x$all.mean,x$mean,
            col=c("black",rep("grey",reps),"blue"),lwd=c(1,rep(1,reps),2),
            main=paste("Forecasts from",toupper(method)))
    # If h==1 then use markers
    if (length(x$mean)==1){
        points(rep(time(x$mean),reps+1),c(x$all.mean,x$mean),
               pch=c(rep(1,reps),20),col=c(rep("grey",reps),"blue"))
    }
    
    
}
    
mlp.thief <- function(y,h=NULL,...){
# This is a wrapper function to use MLP with THieF
    
    # Remove level input from ellipsis
    ellipsis.args <- list(...)
    ellipsis.args$level <- NULL
    ellipsis.args$y <- y
    
    # Fit network
    fit <- do.call(mlp,ellipsis.args) 
    
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

forecast.mlp <- function(object,h=NULL,y=NULL,xreg=NULL,...){
  forecast.net(object,h=h,y=y,xreg=xreg,...)
}

plot.mlp <- function(x, r=1, ...){
    plot.net(x,r)
}

print.mlp <- function(x, ...){
    print.net(x,...)
}

print.net <- function(x, ...){
    
    is.elm.fast <- any(class(x)=="elm.fast")
    difforder <- x$difforder
    sdummy <- x$sdummy
    d <- length(difforder)
    if (is.elm.fast){
        reps <- length(x$b)
    } else {
        reps <- length(x$net$weights)
    }
    hd <- x$hd
    xreg.lags <- x$xreg.lags
    
    if (length(hd)>1 & !is.elm.fast){
        hdt <- paste0(hd,",",collapse="")
        hdt <- paste0("(", substr(hdt,1,nchar(hdt)-1) ,")")
        hde <- "s"
    } else {
        hdt <- hd
        if (is.elm.fast){
            hdt <- paste0(min(hdt)," up to ",max(hdt))
        }
        if (any(hd>1)){
            hde <- "s"
        } else {
            hde <- ""
        }
    }
    
    dtx <- ""
    if (any(class(x)=="elm")){
        if (x$direct == TRUE){
            dtx <- ", direct output connections"
        } 
    } 
    
    method <- class(x)
    if (any(method == "elm")){
        method <- "elm"
    }
    if (is.elm.fast){
        fst <- " (fast)"
    } else {
        fst <- ""
    }
    
    writeLines(paste0(toupper(method),fst," fit with ", hdt," hidden node", hde, dtx," and ", reps, " repetition",if(reps>1){"s"},"."))
    if (d>0){
        writeLines(paste0("Series modelled in differences: ", paste0("D",difforder,collapse=""), "."))
    }
    if (all(x$lags != 0)){
      writeLines(paste0("Univariate lags: (",paste0(x$lags,collapse=","),")"))
    }
    if (!is.null(xreg.lags)){
      null.xreg <- lapply(xreg.lags,length)==0
      p <- length(xreg.lags) - sum(null.xreg)
      if (p > 0){
          if (p == 1){
              rge <- ""
          } else {
              rge <- "s"
          }
          writeLines(paste0(p," regressor",rge," included."))
          pi <- 1
          for (i in which(!null.xreg)){
              writeLines(paste0("- Regressor ",pi," lags: (",paste0(xreg.lags[[i]],collapse=","),")"))
              pi <- pi + 1
          }
      }
    }
    if (sdummy == TRUE){
        writeLines(paste0("Deterministic seasonal dummies included."))
    }
    if (reps>1){
        writeLines(paste0("Forecast combined using the ", x$comb, " operator."))
    }
    if (any(class(x)=="elm")){
        writeLines(paste0("Output weight estimation using: ", x$type, "."))
    }
    writeLines(paste0("MSE: ",round(x$MSE,4),"."))
    
}

get.ff <- function(y){
  # Get time series frequency
  if (any(class(y) == "msts")){
    ff <- attributes(y)$msts
    ff.n <- length(ff)
  } else {
    ff <- frequency(y)
    ff.n <- 1
  }
  return(list("ff"=ff,"ff.n"=ff.n))
}

def.lags <- function(lags,ff,xreg.lags,xreg){
  # Default lagvector
  if (is.null(lags)){
    if (max(ff)>3){
      lags <- 1:max(ff)
    } else {
      lags <- 1:4
    }
  }
  if (!is.null(xreg) && is.null(xreg.lags)){
    x.n <- length(xreg[1,])
    xreg.lags <- rep(list(lags),x.n)
  }
  return(list("lags"=lags,"xreg.lags"=xreg.lags))
}

preprocess <- function(y,m,lags,difforder,sel.lag,allow.det.season,det.type,ff,ff.n,xreg,xreg.lags){
# Pre-process data for MLP and ELM
  
  # Check seasonality & trend
  st <- seasplot(y,m=max(ff),outplot=0)
  if (is.null(st$season.exist)){
    st$season.exist <- FALSE
  }
  
  # Specify differencing order
  difforder <- ndiffs.net(difforder,y,ff,st)
  
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
  
  # Scale xregs and trim initial values for differencing of y
  if (!is.null(xreg)){
    x.n <- dim(xreg)[2]
    xreg.sc <- array(NA,c(dim(xreg)[1]-sum(difforder),x.n))
    xreg.minmax <- vector("list",x.n)
    dstart <- sum(difforder)
    xreg <- xreg[(sum(difforder)+1):dim(xreg)[1],,drop=FALSE]
    for (i in 1:x.n){
      xreg.sc.temp <- linscale(xreg[,i],minmax=list("mn"=-.8,"mx"=0.8))
      xreg.sc[,i] <- xreg.sc.temp$x
      xreg.minmax[[i]] <- xreg.sc.temp$minmax
    }
  } else {
    x.n <- 0
    xreg.minmax <- NULL
    xreg.sc <- xreg
  }
  
  net.inputs <- create.inputs(y.sc, xreg.sc, lags, xreg.lags, n)
  Y <- net.inputs$Y
  X <- net.inputs$X
  Xreg <- net.inputs$Xreg
  lag.max <- net.inputs$lag.max
  rm("net.inputs")
    
  # Create seasonal dummies
  seas.dum <- seas.dum.net(st,difforder,det.type,ff,ff.n,Y,y,allow.det.season)
  Xd <- seas.dum$Xd
  det.type <- seas.dum$det.type
  sdummy <- seas.dum$sdummy
  rm("seas.dum")

  # Select lags
  if (sel.lag == TRUE){
    if (x.n>0){
      Xreg.all <- do.call(cbind,Xreg)
      Xreg.all <- Xreg.all[1:length(Y),,drop=FALSE]
    } else {
      Xreg.all <- NULL
    }
    reg.isel <- as.data.frame(cbind(Y,X,Xreg.all))
    # colnames(reg.isel) <- c("Y",paste0("X",lags),paste0("Xreg",))
    if (sdummy == FALSE){
      # Check if there are no inputs at all
      if (all(colnames(reg.isel) == "Y")){
        stop("Cannot build a network with no univariate or exogenous lags and no deterministic seasonality. Increase the maximum lags.")  
      } else {
        fit <- lm(formula=Y~.,data=reg.isel)
        if (sdummy == FALSE){
          ff.det <- NULL  
        }  
      }
    } else {
      lm.frm <- as.formula(paste0("Y~.+",paste(paste0("Xd[[",1:ff.n,"]]"),collapse="+")))
      fit <- lm(formula=lm.frm,data=reg.isel)
    }
    # Make selection of lags robust to sample size issues if all
    # other checks fail by using lasso
    cf.temp <- tryCatch({
      fit <- stepAIC(fit,trace=0,direction="backward")
      # Get useful lags
      cf.temp <- coef(fit)
    }, error = function(e) {
      lasso.Y <- reg.isel[,1]
      lasso.X <- data.matrix(reg.isel[,2:dim(reg.isel)[2],drop=FALSE])
      if (sdummy == TRUE){
        for (i in 1:ff.n){
          tempX <- Xd[[i]]
          colnames(tempX) <- paste0(paste0("Xd[[",i,"]]"),colnames(tempX))
          lasso.X <- cbind(lasso.X,tempX)
        }
      }
      fit.lasso <- suppressWarnings(cv.glmnet(x=lasso.X,y=lasso.Y))
      cf.temp <- as.vector(coef(fit.lasso))
      names(cf.temp) <- rownames(coef(fit.lasso))
      cf.temp <- cf.temp[cf.temp!=0] 
    })
    X.loc <- lags[which(colnames(X) %in% names(cf.temp))]
    if (x.n>0){
      Xreg.loc <- xreg.lags
      for (i in 1:x.n){
        Xreg.loc[[i]] <- xreg.lags[[1]][which(colnames(Xreg[[i]]) %in% names(cf.temp))]
      }
      xreg.lags <- Xreg.loc
    }

    # Check if deterministic seasonality has remained in the model
    if (sdummy == TRUE){
      still.det <- rep(TRUE,ff.n)
      # Trigonometric dummies will not be retained by linear regression
      # so do not allow rejection by stepwise!
      if (det.type == "bin"){
        for (i in 1:ff.n){
          still.det[i] <- any(grepl(paste0("Xd[[",i,"]]"),names(cf.temp),fixed=TRUE))
        }
      } 
    }
    
    # Although there is an error above to avoid having no inputs, it
    # may still happen if regession rejects all lags. Give a warning!
    # Check if there are any lags
    if (x.n>0){
      # If not univariate and exogenous
      if (sum(c(length(X.loc),unlist(lapply(Xreg.loc,length))))==0){
        # If no deterministic seasonal
        if (sdummy == FALSE){
          warning("No inputs left in the network after pre-selection, forcing AR(1).")
          X.loc <- 1
        }
      }
    } else {
      # If no univariate lags
      if (length(X.loc)==0){
        # If no deterministic seasonal
        if (sdummy == FALSE){
          warning("No inputs left in the network after pre-selection, forcing AR(1).")
          X.loc <- 1
        }
      }
    }
    if (length(X.loc)>0){
      lags <- X.loc
    }
    
    # Recreate inputs
    net.inputs <- create.inputs(y.sc, xreg.sc, lags, xreg.lags, n)
    Y <- net.inputs$Y
    X <- net.inputs$X
    Xreg <- net.inputs$Xreg
    lag.max <- net.inputs$lag.max
    rm("net.inputs")
    
    # Recreate seasonal dummies
    if (sdummy == TRUE){
        # Re-create seasonal dummies
        if (sum(still.det)==0){
            sdummy <- FALSE
            ff.det <- NULL
        } else {
            ff.det <- ff[still.det]
            ff.n.det <- length(ff.det)
            Xd <- vector("list",ff.n.det)
            for (s in 1:ff.n.det){
                Xd[[s]] <- seasdummy(length(Y),y=ts(Y,end=end(y),frequency=ff[s]),type=det.type)
                colnames(Xd[[s]]) <- paste0("D",s,".",1:length(Xd[[s]][1,]))
                if (det.type=="trg"){
                    Xd[[s]] <- Xd[[s]][,1:2]
                }
            }
        }
    }
    
  } else {
    # If no selection is done, match frequencies of dummies with frequencies of time series
    ff.det <- ff
  }
  
  # Merge lags and deterministic seasonality to create network inputs
  if (x.n>0){
    Xreg.all <- do.call(cbind,Xreg)
  } else {
    Xreg.all <- NULL
  }
  X.all <- cbind(X,Xreg.all[1:length(Y),,drop=FALSE])
  if (sdummy == TRUE){
    Xd <- do.call(cbind,Xd)
    X.all <- cbind(X.all,Xd)
  } 

  # Network formula
  frm <- paste0(colnames(X.all),collapse="+")
  frm <- as.formula(paste0("Y~",frm))
  
  return(list("Y"=Y,"X"=X.all,"sdummy"=sdummy,"difforder"=difforder,"det.type"=det.type,"lags"=lags,"xreg.lags"=xreg.lags,"lag.max"=lag.max,"sc"=sc,"xreg.minmax"=xreg.minmax,"d"=d,"y.d"=y.d,"y.ud"=y.ud,"frm"=frm,"ff.det"=ff.det))
  
}

ndiffs.net <- function(difforder,y,ff,st){
  # Find differencing for neural nets
  # NULL is automatic
  # 0 is no differencing
    
  # Find differencing order
  if (all(difforder != 0)){
    if (is.null(difforder)){
      # Identify difforder automatically
      difforder <- 0
      if (st$trend.exist == TRUE){
        difforder <- 1
      }
      if (frequency(y)>1){
        if (st$season.exist == TRUE){
          # difforder <- c(difforder,frequency(y))
          
          # Remove trend appropriately
          cma <- cmav(y,ma=max(ff))
          if (length(y)/max(ff) < 3){
            m.seas <- TRUE
          } else {
            # Test can only run if there are at least three seasons
            m.seas <- mseastest(y,m=max(ff),cma=cma)$is.multiplicative
          }
          if (m.seas == TRUE){
            y.dt <- y/cma
          } else {
            y.dt <- y-cma
          }
          # Check if unit-root stochastic
          d.order <- nsdiffs(ts(y.dt,frequency=max(ff)),test="ch")
          if (d.order > 0){
            difforder <- c(difforder,max(ff))
          } 
        }
      }
    }
  }
    
    # To remove differencing from the remaining it should be set to NULL
    if (any(difforder == 0)){
        difforder <- NULL
    }
  
  return(difforder)
  
}

seas.dum.net <- function(st,difforder,det.type,ff,ff.n,Y,y,allow.det.season){
# Create seasonal dummies for networks
  
  if (if(ff.n > 1){TRUE}else{!any(difforder == max(ff))}
      & frequency(y)>1 & st$season.exist==TRUE & allow.det.season==TRUE){
    sdummy <- TRUE
    # Set type of seasonal dummies
    if (det.type == "auto"){
      if (ff.n == 1 && ff[1] <= 12){
        det.type <- "bin"
      } else {
        det.type <- "trg"
      }
    }
    Xd <- vector("list",ff.n)
    for (s in 1:ff.n){
      Xd[[s]] <- seasdummy(length(Y),y=ts(Y,end=end(y),frequency=ff[s]),type=det.type)
      colnames(Xd[[s]]) <- paste0("D",s,".",1:length(Xd[[s]][1,]))
      if (det.type=="trg"){
        Xd[[s]] <- Xd[[s]][,1:2]
      }
    }
    # Xd <- do.call(cbind,Xd)
    # X <- cbind(X,Xd) 
  } else {
    sdummy <- FALSE
    Xd <- NULL
  }
  
  return(list("Xd"=Xd,"det.type"=det.type,"sdummy"=sdummy))
  
}

create.inputs <- function(y.sc,xreg.sc,lags,xreg.lags,n){
  # Prepare inputs & target
  if (length(lags)>0){
    ylags <- max(lags)
  } else {
    ylags <- 0
  }
  if (!is.null(xreg.sc)){
    xlags <- unlist(lapply(xreg.lags,function(x){if(length(x)){max(x)}else{0}}))
    lag.max <- max(c(ylags,xlags))
  } else {
    lag.max <- ylags
  }
  # Univariate
  if (all(ylags != 0)){
    y.sc.lag <- lagmatrix(y.sc,unique(c(0,lags)))
    Y <- y.sc.lag[(lag.max+1):n,1,drop=FALSE]
    colnames(Y) <- "Y"
    X <- y.sc.lag[(lag.max+1):n,2:(length(lags)+1),drop=FALSE]
    colnames(X) <- paste0("X",lags)
  } else {
    Y <- matrix(y.sc[(lag.max+1):n],ncol=1)
    colnames(Y) <- "Y"
    X <- NULL
  }
  # Exogenous
  if (!is.null(xreg.sc)){
    x.p <- length(xreg.sc[1,])
    x.n <- length(xreg.sc[,1])
    Xreg <- vector("list",x.p)
    for (i in 1:x.p){
      if (length(xreg.lags[[i]]>0)){
        Xreg[[i]] <- lagmatrix(xreg.sc[,i],xreg.lags[[i]])[(lag.max+1):x.n,,drop=FALSE]
        colnames(Xreg[[i]]) <- paste0("Xreg.",i,".",xreg.lags[[i]])
      } else {
        Xreg[[i]] <- NULL
      }
    }
    x.n.all <- lapply(Xreg,function(x){length(x[,1])})
    x.n.all <- x.n.all[x.n.all>0]
    if (any(x.n.all < length(Y))){
      stop("Length of xreg after construction of lags smaller than training sample.")
    }
  } else {
    Xreg <- NULL
  }
  
  return(list("Y"=Y,"X"=X,"Xreg"=Xreg,"lag.max"=lag.max))
  
}

auto.hd.elm <- function(Y,X,frm){

  # Use ELM to find hidden nodes
  sz.elm <- max(10,min(length(X[1,])+2,length(Y)-2))
  reps.elm <- 20
  # sz.elm <- min(40,max(1,length(Y)-2))
  
  net <- neuralnet(frm,cbind(Y,X),hidden=sz.elm,threshold=10^10,rep=reps.elm,err.fct="sse",linear.output=FALSE)
  hd.elm <- vector("numeric",reps.elm)
  for (r in 1:reps.elm){
    Z <- as.matrix(tail(compute(net,X,r)$neurons,1)[[1]][,2:(sz.elm+1)])
    
    type <- "step"
    # Calculate regression
    switch(type,
           "lasso" = {
             fit <- suppressWarnings(cv.glmnet(Z,cbind(Y)))
             cf <- as.vector(coef(fit))
             hd.elm[r] <- sum(cf != 0)-1 # -1 for intercept
           },
           {
             reg.data <- as.data.frame(cbind(Y,Z))
             colnames(reg.data) <- c("Y",paste0("X",1:sz.elm))
             # Take care of linear dependency
             alias.fit <- alias(Y~.,data=reg.data)
             alias.x <- rownames(alias.fit$Complete)
             frm.elm <- as.formula(paste0("Y~",paste0(setdiff(colnames(reg.data)[2:(sz.elm+1)],alias.x),collapse="+")))
             fit <- suppressWarnings(lm(frm.elm,reg.data))
             if (type == "step"){
               fit <- suppressWarnings(stepAIC(fit,trace=0,direction="backward"))
             }
             hd.elm[r] <- sum(summary(fit)$coefficients[,4]<0.05,na.rm=TRUE)-(summary(fit)$coefficients[1,4]<0.05)
           })
    
  }
  hd <- round(median(hd.elm))
  if (hd<1){
    hd <- 1
  }
  
  return(hd)

}

auto.hd.cv <- function(Y,X,frm,comb,reps,type=c("cv","valid"),hd.max=NULL){
  # Find number of hidden nodes with CV
  
  # Setup
  type <- type[1]
  K <- 5                                            # Number of folds
  val.size <- 0.2                                   # Size of validation set
  reps <- min(c(20,max(c(2,reps))))                 # Number of NN reps, maximum 20
  if (is.null(hd.max)){
    hd.max <- max(2,min(length(X[1,])+2,length(Y)-2)) # Maximum number of hidden nodes
  }
  
  # Setup folds or validation set
  n <- length(Y)
  if (type == "cv"){
    # Create folds
    if (K >= n){
      stop("Too few observations to perform cross-validation for specification of hidden nodes.")
    }
    # Create fold indices
    idx.all <- sample(1:n)
    cv.cut <- seq(0,n,length.out=K+1)
    idx <- vector("list",K)
    for (i in 1:K){
      idx[[i]] <- which((idx.all <= cv.cut[i+1]) & (idx.all > cv.cut[i]))
    }  
  } else {
    # Validation set
    K <- 1
    idx <- list(sample(1:n,round(val.size*n)))
  }
  
  # Now run CV (1-step ahead)
  err.h <- array(NA,c(hd.max,1),dimnames=list(paste0("H.",1:hd.max),"MSE"))
  for (h in 1:hd.max){
    # For each fold
    err.cv <- vector("numeric",K)
    for (i in 1:K){
      Y.trn <- Y[setdiff(1:n,idx[[i]]),,drop=FALSE]
      X.trn <- X[setdiff(1:n,idx[[i]]),,drop=FALSE]
      Y.tst <- Y[idx[[i]],,drop=FALSE]
      X.tst <- X[idx[[i]],,drop=FALSE]
      net <- neuralnet(frm,cbind(Y.trn,X.trn),hidden=h,rep=reps,err.fct="sse",linear.output=TRUE)
      reps <- length(net$weights) # In case some network is untrained
      
      # For each training repetition
      frc <- array(NA,c(length(Y.tst),reps))
      for (r in 1:reps){
        frc[,r] <- compute(net,X.tst,r)$net.result  
      }
      frc <- frc.comb(frc,comb)
      err.cv[i] <- mean((Y.tst - frc)^2)
    }
    err.h[h] <- mean(err.cv)
  }
  hd <- which(err.h == min(err.h))[1]

  return(list("hd"=hd,"mseH"=err.h))
  
}

frc.comb <- function(Yhat,comb){
  # Combine forecasts
  r <- length(Yhat[1,])
  if (r>1){
    switch(comb,
           "median" = {yout <- apply(Yhat,1,median)},
           "mean" = {yout <- apply(Yhat,1,mean)},
           "mode" = {yout <- sapply(apply(Yhat,1,kdemode),function(x){x[[1]][1]})}
    )
  } else {
    yout <- Yhat[,1]
  }
  
  return (yout)
}