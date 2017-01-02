seasdummy <- function(n,m=12,y=NULL,type=c("bin","trg"),full=c(FALSE,TRUE)){
# Create deterministic seasonality dummies

    # Default
    type <- type[1]
    full <- full[1]
    
    # Get time series information if available
    if (!is.null(y)){
        m <- frequency(y)
        start <- tail(start(y),1)
    } else {
        start <- 1
    }
    
    # Create dummies
    x <- array(0,c(n,m))
    if (type == "bin"){
        x[seq(1,n,m),1] <- 1
        for (i in 2:m){
            x[,i] <- c(tail(x[,i-1],1),x[1:(n-1),i-1])
        }
    } else { # trg
        t <- 1:n
        for (i in 1:(m/2)){
            x[,1+(i-1)*2] <- cos((2*t*pi*i)/m)
            x[,2+(i-1)*2] <- sin((2*t*pi*i)/m)
        }
    }
    
    # Trim co-linear dummy
    if (full==FALSE){
        x <- x[,1:(m-1)]
    }
    
    # Shift for start
    if (start > 1){
        x <- rbind(x[start:n,],x[1:(start-1),])
    }
    
    return(x)
    
}