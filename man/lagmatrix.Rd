\name{lagmatrix}
\alias{lagmatrix}

\title{Create matrix with lead/lags of an input vector}

\description{Create matrix with lead/lags of an input vector.}

\usage{
lagmatrix(x,lag)
}

\arguments{
  \item{x}{
    Input.
    }
  \item{lag}{
    Vector or leads and lags. Positive numbers are lags, negative are leads. O is the original x.
    }
}
\value{
\item{lmat}{Resulting matrix.}
}
\author{
Nikolaos Kourentzes
}
\examples{
x <- rnorm(10)
lagmatrix(x,c(0,1,-1))
}