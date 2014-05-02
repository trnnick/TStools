wins <- function (x, p=0.05){
# Winsorise a vector
#
# Inputs
#   x       Vector to be winsorised
#   p       Percentage or number of observations to be winsorised. 
#           If value is <1 then it is used as a percentages. Otheriwse it is the number
#           of observations to winsorise. If the resulting p > floor((length(x)-1)/2) 
#           then it is set equal to floor((length(x)-1)/2).
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
  n <- sum(is.na(x))
  c <- l - n
  
  # Convert percentage to number of observations  
  if (p<1){
    p <- ceiling(c*p)
  }
  
  # At maximum floor((l-1)/2) observations are taken from each side
  if (p>floor((c-1)/2)){
    p <- floor((c-1)/2)
  }
  
  # Sort and find wmin and wmax
  x.sort <- sort(x)
  x.wmin <- x.sort[p+1]
  x.wmax <- x.sort[c-p]
  
  # Winsorise (ignore NAs)
  x.out <- x
  idx <- !is.na(x)
  x.w <- x.out[idx]
  x.w[x[idx]<x.wmin] <- x.wmin
  x.w[x[idx]>x.wmax] <- x.wmax
  x.out[idx] <- x.w
  
  return(x.out)
  
}

colWins <- function (x, p=0.05){
# Winsorise by columns
#
# Inputs
#   x       Array/matrix to be winsorised
#   p       Percentage or number of observations to be winsorised. 
#           If value is <1 then it is used as a percentages. Otheriwse it is the number
#           of observations to winsorise. If the resulting p > floor((length(x)-1)/2) 
#           then it is set equal to floor((length(x)-1)/2).#
# Output
#   x.out   Winsorised array/matrix.
#
# Example
#   x <- matrix(rnorm(10*20,mean=0,sd=1),10,20) 
#   x.out <- colWins(x)
#
# Nikolaos Kourentzes, 2014 <nikolaos@kourentzes.com>

  k <- dim(x)[2]
  v.wins <- Vectorize(function(i){temp <- wins(x[,i], p)})
  x.out <- v.wins(1:k)
  x.out <- x.out + 0*x
  
  return(x.out)
  
}

rowWins <- function (x, p=0.05){
# Winsorise by rows
#
# Inputs
#   x       Array/matrix to be winsorised
#   p       Percentage or number of observations to be winsorised. 
#           If value is <1 then it is used as a percentages. Otheriwse it is the number
#           of observations to winsorise. If the resulting p > floor((length(x)-1)/2) 
#           then it is set equal to floor((length(x)-1)/2).#
# Output
#   x.out   Winsorised array/matrix.
#
# Example
#   x <- matrix(rnorm(10*20,mean=0,sd=1),10,20) 
#   x.out <- rowWins(x)
#
# Nikolaos Kourentzes, 2014 <nikolaos@kourentzes.com>
  
  d <- dim(x)
  v.wins <- Vectorize(function(i){temp <- wins(x[i,], p)})
  x.out <- v.wins(1:d[1])
  x.out <- array(t(x.out),d)
  x.out <- x.out + 0*x
  
  return(x.out) 
  
}