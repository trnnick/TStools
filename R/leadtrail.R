leadtrail <- function(x,rm=c("zeros","na"),lead=c(TRUE,FALSE),trail=c(TRUE,FALSE)){
# Function to remove leading/trailing zeros or NAs
#
# Inputs: 
#   x         Vector of values to check.
#   rm        What to remove, can be "zeros" or "na".
#   lead      Remove leading values: TRUE or FALSE.
#   trail     Remove trailing values: TRUE or FALSE.
#
# Output:
#   y         Resulting vector
#
# Example:
#   x <- c(rep(0,5),rnorm(100),rep(0,5))
#   y <- leadtrail(x)
#   print(x)
#   print(y)
#
# Nikolaos Kourentzes, 2016 <nikolaos@kourentzes.com>
  
  # Defaults
  rm <- rm[1]
  lead <- lead[1]
  trail <- trail[1]
  
  # Select what to remove
  if (rm=="zeros" | rm=="z"){
    idx <- which(x == 0)
  } else if (rm=="na" | rm=="n"){
    idx <- which(is.na(x))
  } else {
    stop("Incorrect rm choice.")
  }
  
  n <- length(x)
  l <- length(idx)
  
  # Handle leading observations
  if (lead==TRUE){
    
    if (idx[1]==1){
      d.idx <- diff(idx)
      loc <- which(d.idx > 1)[1]
      if (is.na(loc)){
        loc <- l
      }
      lead.rm <- 1:loc
    } else {
      lead.rm <- NULL
    }
    
  } else {
    lead.rm <- NULL
  }
  
  # Handle trailing observations
  if (trail==TRUE){
    
    if (tail(idx,1)==n){
      d.idx <- diff(rev(idx))
      loc <- which(d.idx != -1)[1]
      if (is.na(loc)){
        loc <- l
      }
      trail.rm <- (n-loc+1):n
    } else {
      trail.rm <- NULL
    }
    
  } else {
    trail.rm <- NULL
  }
  
  keep <- rep(TRUE,n)
  keep[lead.rm] <- FALSE
  keep[trail.rm] <- FALSE
  
  y <- x[keep]
  return(y)
  
}

