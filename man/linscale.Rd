\name{linscale}
\alias{linscale}

\title{ Apply minmax linear scaling to a vector}

\description{Apply minmax linear scaling to a vector.}

\usage{
linscale(x,minmax=NULL,rev=c(FALSE,TRUE))
}

\arguments{
  \item{x}{
    Input vector.
    }
  \item{minmax}{
    minmax must be a list with elements "mn", "mx", "mn.orig" and "mx.orig", where "mn" and "mx" refer to the target min and max, and the remaining two refer to the current vector min and max. By default mn=-1 and mx=1. mn.orig and mx.orig can be missing, unless the scaling is reversed.
    }
  \item{rev}{
    Reverse scaling back to original: TRUE or FALSE.
    }
}
\value{
\item{x}{Scaled vector.}
\item{minmax}{List with resulting mn, mx, mn.orig and mx.orig. Can be used as input to reverse scaling.}
}
\author{
Nikolaos Kourentzes
}
\examples{
y <- rnorm(20)*100
sc <- linscale(y)
x <- sc$x
print(c(min(y),max(y)))
print(c(min(x),max(x)))
sc.rev <- linscale(x,minmax=sc$minmax,rev=TRUE)
print(c(min(sc.rev$x),max(sc.rev$x)))
}