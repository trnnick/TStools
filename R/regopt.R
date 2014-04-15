regopt <- function(y,X,constant=c(TRUE,FALSE),cost=c("MAE","MdAE","MSE","MdSE","ME","MdE"),
                   outplot=c(FALSE,TRUE))
# Identify regression beta using different cost functions
#
# Inputs:
#   y           Vector of target data, can be ts object
#   X           Matrix of explanatory data, each variable is a column
#   constant    If TRUE then a constant is added to the model
#   cost        Cost function to use:
#                 MAE  - Mean Absolute Error [Default]
#                 MdAE - Median Absolute Error
#                 MSE  - Mean Squared Error
#                 MdSE - Median Squared Error
#                 ME   - Mean Error
#                 MdE  - Median Error
#   outplot     If TRUE plot regression fit
#
# Outputs:
#   b           Regression coefficients
#               yhat <- X %*% b
#
# Example:
#   y <- wineind
#   X <- c(1:length(y))
#   X <- as.matrix(X,nrow=length(y),ncol=1)
#   
#   b <- regopt(y,X,cost="MAE",outplot=TRUE)
#   
#   # Add a constant to the X matrix and calculate yhat
#   yhat <- cbind(as.matrix(rep(1,times=length(y)),nrow=length(y),ncol=1),X)%*%b
#   
#   # Calculate some errors
#   ME <- mean(y-yhat)
#   MAE <- mean(abs(y-yhat))
#   RMSE <- sqrt(mean((y-yhat)^2))
#   
#   print(paste("ME:",round(ME,2),"MAE",round(MAE,2),"RMSE",round(RMSE,2)))
#
# Nikolaos Kourentzes, 2014
{

  # Set defaults
  cost <- cost[1]
  constant <- constant[1]
  outplot <- outplot[1]
  
  nY <- length(y)
  nX <- nrow(X)
  k <- ncol(X)
  
  # Check if sample size of y and X is the same
  if (nY != nX){
    stop("Length of y and rows of X do not match.")
  }
  
  # Add constant if requested
  if (constant == TRUE){
    C <- as.matrix(rep(1,times=nX),nrow=nX,ncol=1)
    X <- cbind(C,X)
  }
  
  p <- k + constant
  
  # Initial beta values
  b0 <- solve(t(X)%*%X)%*%t(X)%*%y
  # b0 <- as.matrix(rep(0,times=p),nrow=p,ncol=1)
  
  if (cost != "MSE"){
    # Cost functions
    if (cost == "MAE"){
      # Mean absolute error 
      costf <- function(b,y,X)mean(abs(y - X%*%b))
    } else if (cost == "MdAE"){
      # Median absolute error
      costf <- function(b,y,X)median(abs(y - X%*%b))
    } else if (cost == "MdSE"){
      # Median squared error
      costf <- function(b,y,X)median((y - X%*%b)^2)
    } else if (cost == "ME"){
      # Mean error
      costf <- function(b,y,X)abs(mean(y - X%*%b))
    } else if (cost == "MdE"){
      # Median error
      costf <- function(b,y,X)abs(median(y - X%*%b))
    }
    # Optimise beta
    opt <- optim(par=b0, fn=costf, method = "Nelder-Mead", y=y, X=X, control = list(maxit = 2000))
    b <- opt$par
    # cf <- opt$value
    # cf
    
  } else {
    # Mean square error
    b <- b0
  }
  
  b <- as.matrix(b,nrow=p,ncol=1)
  
  # Produce plot
  if (outplot==TRUE){
    yhat <- X%*%b
    yols <- X%*%b0
    if (class(y) == "ts"){
      yhat <- ts(yhat,frequency=frequency(y),start=start(y))
      yols <- ts(yols,frequency=frequency(y),start=start(y))
      st <- start(y)
    } else {
      st <- 1
    }
    plot(y)
    lines(yhat,col="blue")
    lines(yols,col="red")
    legend("topleft",c("Data","y-hat","y-ols"),lty=c(1,1),lwd=c(2.5,2.5),col=c("black","blue",'red'),bty="n")
  }
  
  return(b)
  
}