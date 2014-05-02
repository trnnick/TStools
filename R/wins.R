wins <- function (x, prc=0.05, n=NULL){
# Winsorise a vector
#
# Inputs
#   x       Vector to be winsorised
#   prc     Percentage of observations to be winsorised. 
#           Values about 0.5 are overriden to 0.5
#   n       Number of observations to be winsorised. If given this overrides prc. 
#           If n > floor((length(x)-1)/2) then it is set equal to floor((length(x)-1)/2).
#
# Output
#   x.out   Winsorised vector.
#
# Example
#   x <- rnorm(100,mean=0,sd=1)
#   x.out <- wins(x)
#
# Nikolaos Kourentzes, 2014 <nikolaos@kourentzes.com>
  
  l <- length(x)
  
  # Convert percentage to number of observations  
  # n overrides prc
  if (is.null(n)){
    n <- ceiling(l*prc)
  }
  
  # At maximum floor((l-1)/2) observations are taken from each side
  if (n>floor((l-1)/2)){
    n <- floor((l-1)/2)
  }
  
  # Sort and find wmin and wmax
  x.sort <- sort(x)
  x.wmin <- x.sort[n+1]
  x.wmax <- x.sort[l-n]
  
  # Winsorise
  x.out <- x
  x.out[x<x.wmin] <- x.wmin
  x.out[x>x.wmax] <- x.wmax
  
  return(x.out)
  
}

colWins <- function (x, prc=0.05, n=NULL){
# Winsorise by columns
#
# Inputs
#   x       Array/matrix to be winsorised
#   prc     Percentage of observations to be winsorised per column.
#           Values about 0.5 are overriden to 0.5
#   n       Number of observations to be winsorised per column. If given this overrides prc. 
#           If n > floor((length(x)-1)/2) then it is set equal to floor((length(x)-1)/2).
#
# Output
#   x.out   Winsorised array/matrix.
#
# Example
#   x <- matrix(rnorm(10*20,mean=0,sd=1),10,20) 
#   x.out <- colWins(x)
#
# Nikolaos Kourentzes, 2014 <nikolaos@kourentzes.com>

  k <- dim(x)[2]
  v.wins <- Vectorize(function(i){temp <- wins(x[,i], prc, n)})
  x.out <- v.wins(1:k)
  x.out <- x.out + 0*x
  
}

rowWins <- function (x, prc=0.05, n=NULL){
# Winsorise by rows
#
# Inputs
#   x       Array/matrix to be winsorised
#   prc     Percentage of observations to be winsorised per row.
#           Values about 0.5 are overriden to 0.5
#   n       Number of observations to be winsorised per row. If given this overrides prc. 
#           If n > floor((length(x)-1)/2) then it is set equal to floor((length(x)-1)/2).
#
# Output
#   x.out   Winsorised array/matrix.
#
# Example
#   x <- matrix(rnorm(10*20,mean=0,sd=1),10,20) 
#   x.out <- rowWins(x)
#
# Nikolaos Kourentzes, 2014 <nikolaos@kourentzes.com>
  
  d <- dim(x)
  v.wins <- Vectorize(function(i){temp <- wins(x[i,], prc, n)})
  x.out <- v.wins(1:d[1])
  x.out <- array(t(x.out),d)
  x.out <- x.out + 0*x
  
}