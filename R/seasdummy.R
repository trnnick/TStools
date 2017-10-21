seasdummy <- function(n,m=NULL,y=NULL,type=c("bin","trg"),full=c(FALSE,TRUE)){
# Create deterministic seasonality dummies

    # Default
    type <- match.arg(type,c("bin","trg"))
    full <- full[1]

    # Get time series information if available
    if (!is.null(y)){
      if (is.null(m)){
        m <- frequency(y)
      }
      start <- tail(start(y),1)
    } else {
      start <- 1
    }
    
    if (start >= n){
      n.sim <- n + start
    } else {
      n.sim <- n
    }
    
    if (is.null(m)){
      stop("Seasonality not specified.")
    }
    
    # Create dummies
    if (type == "bin"){
        x <- matrix(rep(diag(rep(1,m)),ceiling(n.sim/m)),ncol=m,byrow=TRUE)[1:n.sim,,drop=FALSE]
    } else { # trg
        x <- array(0,c(n.sim,m))
        t <- 1:n.sim
        for (i in 1:(m/2)){
            x[,1+(i-1)*2] <- cos((2*t*pi*i)/m)
            x[,2+(i-1)*2] <- sin((2*t*pi*i)/m)
        }
    }
    
    # Trim co-linear dummy
    if (full==FALSE){
        x <- x[,1:(m-1),drop=FALSE]
    }
    
    # Shift for start
    if (start > 1){
        x <- rbind(x[start:n, ,drop=FALSE], x[1:(start - 1), ,drop=FALSE])
    }
    
    # If n.sim is larger, just retain n observations
    x <- x[1:n,,drop=FALSE]
    
    return(x)
    
}