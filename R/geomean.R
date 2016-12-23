geomean <- function(x,na.rm=c(FALSE,TRUE)){
    
    na.rm <- na.rm[1]
    gm <- exp(mean(log(x),na.rm=na.rm))
    return(gm)
    
}