\name{leadtrail}
\alias{leadtrail}

\title{Remove leading/trailing zeros or NAs}

\description{Function to remove leading/trailing zeros or NAs.}

\usage{
leadtrail(x,rm=c("zeros","na"),lead=c(TRUE,FALSE),trail=c(TRUE,FALSE))
}

\arguments{
  \item{x}{
    Vector of values to check.
    }
  \item{rm}{
    What to remove, can be "zeros" or "na".
    }
  \item{lead}{
    Remove leading values: TRUE or FALSE.
    }
  \item{trail}{
    Remove trailing values: TRUE or FALSE.
    }
}
\value{
\item{y}{Resulting vector.}
}
\author{
Nikolaos Kourentzes
}
\examples{
x <- c(rep(0,5),rnorm(100),rep(0,5))
y <- leadtrail(x)
print(x)
print(y)
}