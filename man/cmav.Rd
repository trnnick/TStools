\name{cmav}
\alias{cmav}

\title{Centred moving average}

\description{Calculate centred moving average for a time series.}

\usage{
cmav(y,ma=NULL,fill=c(TRUE,FALSE),outplot=c(FALSE,TRUE),fast=c(TRUE,FALSE))
}

\arguments{
  \item{y}{
    A time series. Can be a time series object.
}
  \item{ma}{
    Length of centred moving average.
    If y is a ts object then this is taken from its frequency, unless overriden.
}
  \item{fill}{
    If TRUE then fill first and last ma/2 observations using ETS.
}
  \item{outplot}{
    If TRUE then produce plot.
}
  \item{fast}{
    If TRUE then only a limited set of models are evaluated for CMA extrapolation.
}
}
\value{
\item{cma}{Centred moving average. If y is a ts object, then cma has the same properties}
}
\author{
Nikolaos Kourentzes
}
\examples{
cma <- cmav(referrals)
}
\keyword{ ~ets }
