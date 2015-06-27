\name{stlboot}
\alias{stlboot}

\title{STL/Loess bootstrapping}

\description{STL/Loess bootstrapping by Bergmeir, Hyndman and Benitez, 2014.}

\usage{
stlboot(ts,k=1,test.season=c(TRUE,FALSE),outplot=c(FALSE,TRUE))
}

\arguments{
  \item{ts}{
    Time series to bootstrap. Must be ts object.
  }
  \item{k}{
    Number of bootstraps.
  }
  \item{test.season}{
    If TRUE then test for presence of seasonality. If FALSE then all non-yearly are decomposed using STL, assuming seasonality.
  }
  \item{outplot}{
    If TRUE provide a plot of the bootstrapped series.
}
}
\value{
\item{ts.recon}{Array of bootstrapped series. Each column is a time series.}
}
\references{
Bergmeir C.,  Hyndman R. J., Benitez J. M., Bagging Exponential Smoothing Methods using STL Decomposition and Box-Cox Transformation. 2014, Working Paper.
http://robjhyndman.com/working-papers/bagging-ets/
}
\author{
Nikolaos Kourentzes
}
\examples{
stlboot(referrals,k=20,outplot=TRUE)
}
\keyword{ ~stlboot }
\seealso{
\code{\link{stl}}, \code{\link{loess}}.
}