sma <- function(x, order=NULL, h=10, optimize=FALSE,
                min.order=1, max.order=12, holdout=FALSE,
                silent=FALSE, legend=TRUE) {
# This function computes a one-step-ahead simple moving average for a time series and returns the forecasts.
# If optimize = TRUE then the best moving average order is determined by minimizing the Mean Squared Error on the training set.

# Check if the data is vector
  if(!is.numeric(x) & !is.ts(x)){
    stop("The provided data is not a vector or ts object! Can't do anything!", call.=FALSE)
  }
  
# Check if horizon exists
  if(!is.numeric(h)){
    stop ("The forecasting horizon should be a number", call.=FALSE)
  }
  else if(h %% 1 != 0 | h <= 0){
    stop ("The forecasting horizon should be a non-negative integer", call.=FALSE)
  }

  if(is.null(order) & optimize==FALSE){
    optimize <- TRUE;
  }
  
# Define n.all, the overal number of observations (in-sample + holdout)
  n.all <- length(x) + (1 - holdout)*h

# Define n, the number of observations of in-sample
  n <- length(x) - holdout*h
  
  if(optimize == TRUE){
# Check if Minimum order is correct
    if(is.numeric(min.order) == FALSE || min.order %% 1 != 0 || min.order <= 0){ 
      stop("The minimum order should be a non-negative integer", call.=FALSE)
    }
# Check if Maximum order is correct
    if(is.numeric(max.order) == FALSE || max.order %% 1 != 0 || max.order <= 0){
      stop("The maximum order should be a non-negative integer", call.=FALSE)  
    }
# Check that the maximum order is greater than the minimum order
    if(min.order >= max.order){
      stop("The minimum order should be less than the maximum order", call.=FALSE)
    }
# Check if minimum and maximum order do not exceed half the order of the data. This is just a precaution
    if(min.order > n/2 | max.order > n/2){
      stop("The minimum and maximum can not exceed half the order of the data", call.=FALSE)
    }

    l_span <- seq(min.order, max.order, by = 1); #Create a vector to store all the sma orders
    MSE <- NULL;
    min_MSE <- NULL;
    best_order <- NULL;
    trn.sma.best <- NULL;

# Start optimisation cycle
    for(i in l_span){
# Create a vector to store insample forecasts
      trn.sma <- rep(NA,n);

      for(j in 1:(n - i)){
          trn.sma[j + i] <- mean(x[(max.order+j):(i + j - 1)]);
      }

# Calculate residuals
      errs <- na.omit(x[1:n] - trn.sma);
      MSE <- mean(errs^2);

      if(i == min.order){
# Initialize
        best_order <- i;
        min_MSE <- MSE;
        trn.sma.best <- trn.sma;
      }

      if(MSE < min_MSE){
        best_order <- i;
        min_MSE <- MSE;
        trn.sma.best <- trn.sma;
      }
    }
      
    order <- best_order;
    trn.sma <- trn.sma.best;
  }
    
  if(optimize == FALSE){
# Check if order is correct
    if(order > n/2){
      stop("The order can not exceed half the order of the data.")
    }
    trn.sma <- rep(NA, n);
# Create a vector to store insample forecasts
    for(j in 1:(n - order)){
      trn.sma[j + order] <- mean(x[(j):(order + j - 1)]);
    }
  }

  trn.sma <- ts(trn.sma,start=start(x),frequency=frequency(x))
  tst.frcst <- ts(rep(mean(x[(n - order + 1):(n)]),h),start=time(x)[n]+deltat(x),frequency=frequency(x))

  if(holdout==TRUE){
    x.holdout <- ts(x[(n+1):n.all],start=start(tst.frcst),frequency=frequency(x));
    errormeasures <- c(MAPE(as.vector(x.holdout),as.vector(tst.frcst),round=5),
                       MASE(as.vector(x.holdout),as.vector(tst.frcst),mean(abs(diff(as.vector(x)[1:n])))),
                       MASE(as.vector(x.holdout),as.vector(tst.frcst),mean(abs(as.vector(x)[1:n]))),
                       MPE(as.vector(x.holdout),as.vector(tst.frcst),round=5),
                       SMAPE(as.vector(x.holdout),as.vector(tst.frcst),round=5));
    names(errormeasures) <- c("MAPE","MASE","MASALE","MPE","SMAPE");
  }
  else{
    x.holdout <- NA;
    errormeasures <- NA;
  }

  if(silent==FALSE){
    print(paste0("Method used: SMA(",order,")"));
    if(holdout==TRUE){
        print(paste0("MPE: ",errormeasures["MPE"]*100,"%"));
        print(paste0("MAPE: ",errormeasures["MAPE"]*100,"%"));
        print(paste0("SMAPE: ",errormeasures["SMAPE"]*100,"%"));
        print(paste0("MASE: ",errormeasures["MASE"]));
    }
    
    plot.range <- range(min(x[!is.na(x)],trn.sma[!is.na(trn.sma)],tst.frcst[!is.na(tst.frcst)]),
                        max(x[!is.na(x)],trn.sma[!is.na(trn.sma)],tst.frcst[!is.na(tst.frcst)]))

    par(mfrow=c(1,1), mar=c(5,3,2,1))
    plot(x,type="l",xlim=range(time(x)[1],time(tst.frcst)[h]),
         ylim=plot.range,xlab="Time", ylab="")
    lines(trn.sma,col="purple",lwd=2,lty=2)
    
    if(h>1){
      lines(tst.frcst,col="blue",lwd=2)

      if(legend==TRUE){
# Define where to position the legend
        if(mean(c(trn.sma[!is.na(trn.sma)],tst.frcst)[1:round(n.all/3,0)])<(plot.range[2]+plot.range[1])/2){
          leg.place = "topleft"
        }
        else{
          leg.place = "bottomleft"
        }

        legend(x=leg.place,
               legend=c("Series","Fitted values","Point forecast","Forecast origin"),
               col=c("black","purple","blue","red"),
               lwd=c(1,2,2,2),
               lty=c(1,2,1,1))
      }
    }
    else{
      points(tst.frcst,col="blue",lwd=2,pch=4)
      if(legend==TRUE){
# Define where to position the legend
        if(mean(c(trn.sma[!is.na(trn.sma)],tst.frcst)[1:round(n.all/3,0)])<(plot.range[2]+plot.range[1])/2){
          leg.place = "topleft"
        }
        else{
          leg.place = "bottomleft"
        }

        legend(x=leg.place,
               legend=c("Series","Fitted values","Point forecast","Forecast origin"),
               col=c("black","purple","blue","red"),
               lwd=c(1,2,2,2),
               lty=c(1,2,NA,1),
               pch=c(NA,NA,4,NA))
      }
    }
  }

  return(list(order = order, fitted = trn.sma, forecast = tst.frcst));
}