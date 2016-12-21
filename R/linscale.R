linscale <- function(x,minmax=NULL,rev=c(FALSE,TRUE)){
# Apply minmax linear scaling to a vector
#
# Inputs:
#   x       Input vector.
#   minmax  minmax must be a list with elements "mn", "mx", "mn.orig" and "mx.orig", 
#           where "mn" and "mx" refer to the target min and max, and the remaining two 
#           refer to the current vector min and max. By default mn=-1 and mx=1. mn.orig 
#           and mx.orig can be missing, unless the scaling is reversed. 
#   rev     Reverse scaling back to original.
#
# Outputs:
#   x       Scaled vector
#   minmax  List with resulting mn, mx, mn.orig and mx.orig. Can be used as input to 
#           reverse scaling.
#
# Example:
#   y <- rnorm(20)*100
#   sc <- linscale(y)
#   x <- sc$x
#   print(c(min(y),max(y)))
#   print(c(min(x),max(x)))
#   sc.rev <- linscale(x,minmax=sc$minmax,rev=TRUE)
#   print(c(min(sc.rev$x),max(sc.rev$x)))
#
# Nikolaos Kourentzes, 2016 <nikolaos@kourentzes.com>
    
    rev <- rev[1]
    
    # For reversing the scaling both scaled and original mn/mx are required to allow scaling
    # samples with mn/mx different from the one used to infer the scaling.
    
    if (is.null(minmax)){
        # minmax is not given
        if (rev == TRUE){
            stop("To reverse scaling minmax input is required.")
        }
        minmax <- list(mn=-1,mx=1,mn.orig=NULL,mx.orig=NULL)
        mn <- minmax$mn
        mx <- minmax$mx
    } else {
        # minmax is given
        names.minmax <- names(minmax)
        if (rev == TRUE){
            # reverse scaling
            if (any(names.minmax == "mn.orig") & any(names.minmax == "mx.orig") &
                any(names.minmax == "mn") & any(names.minmax == "mx")){
                mn <- minmax$mn
                mx <- minmax$mx
                mn.orig <- minmax$mn.orig
                mx.orig <- minmax$mx.orig
            } else {
                stop("Provided minmax list is not of correct type. It must contain mn, mx, mn.orig and mx.orig to reverse scaling.")
            }
        } else {
            # apply scaling
            if (any(names.minmax == "mn") & any(names.minmax == "mx")){
                mn <- minmax$mn
                mx <- minmax$mx
            } else {
                stop("Provided minmax list is not of correct type. It must contain mn and mx to apply scaling.")
            }
        }
    }
    
    if (rev == FALSE){
        # Apply scaling
        mx.orig <- max(x)
        mn.orig <- min(x)
        x.sc <- (mx-mn)*(x-mn.orig)/(mx.orig-mn.orig)+mn
        return(list("x"=x.sc,"minmax"=list("mn"=mn,"mx"=mx,"mn.orig"=mn.orig,"mx.orig"=mx.orig)))
    } else {
        # Reverse scaling
        x.orig <- (mx.orig-mn.orig)*(x-mn)/(mx-mn)+mn.orig
        return(list("x"=x.orig,"minmax"=list("mn"=mn.orig,"mx"=mx.orig,"mn.orig"=mn,"mx.orig"=mx)))
    }
    
}

