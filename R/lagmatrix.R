lagmatrix <- function(x,lag){
# Create matrix with lead/lags of an input vector
# 
# Inputs:
#   x       Input.
#   lag     Vector or leads and lags. Positive numbers are lags, negative are leads. O is the original x.
#
# Output:
#   lmat    Resulting matrix.
#
# Example:
#   x <- rnorm(10)
#   lagmatrix(x,c(0,1,-1))
#
# Nikolaos Kourentzes, 2016 <nikolaos@kourentzes.com>
    
    # Construct matrix with lead and lags
    n <- length(x)
    k <- length(lag)
    # How much to expand for leads and lags
    mlg <- max(c(0,lag[lag>0]))
    mld <- max(abs(c(0,lag[lag<0])))
    # Assign values
    lmat <- array(NA,c(n+mlg+mld,k))
    for (i in 1:k){
        lmat[(1+lag[i]+mld):(n+lag[i]+mld),i] <- x
    }
    # Trim lmat for expansion
    lmat <- lmat[(mld+1):(mld+n),,drop=FALSE]
    return(lmat)
}