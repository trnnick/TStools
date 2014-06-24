\name{nzlength}
\alias{nzlength}

\title{Non-zero length}

\description{Find dataset series non-zero length and identify series over or equal to threshold.}

\usage{
nzlength(y,t=NULL,outplot=c(FALSE,TRUE))
}

\arguments{
  \item{y}{
    Dataset, each column is a series.
}
  \item{t}{
    Threshold to compare length with. Use NULL if not needed.
}
  \item{outplot}{
    If TRUE produce plot of the cumulative distribution of lengths.
}
}
\value{
\item{data.length}{ector with non-zero length of each time series.}
\item{index}{Logical vector of series that are above or equal to threshold.}
}
\author{
Nikolaos Kourentzes
}
\examples{
x <- abs(matrix(cumsum(rnorm(5400,0,1)),36,150))
x[1:20,1:50] <- 0
nzlength(x,30,outplot=TRUE)
}

