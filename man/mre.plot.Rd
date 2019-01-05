\name{mre.plot}
\alias{mre.plot}

\title{Mean Root Error Bias Plot}

\description{Plots the `Bias Plot' for MRE.}

\usage{
mre.plot(mre,main=NULL,plot.legend=c(TRUE,FALSE),maxmag=NULL,ts=1,group=NULL)
}

\arguments{
  \item{mre}{
    k Root or Mean Root Errors to be plotted. These must already be complex numbers.
}
  \item{main}{
    Main title for plot.
}
  \item{plot.legend}{
    If TRUE plot legend when k>1.
}
  \item{maxmag}{
    Maximum MRE magnitude.
}
  \item{ts}{
    Text size cex.
}
  \item{group}{
    Vector containing either characters or numerics to group errors.
}
}
\references{
Kourentzes N., Trapero J. R., Svetunkov I. Measuring the behaviour of experts on demand forecasting: a complex task. Working paper
http://kourentzes.com/forecasting/2014/12/17/measuring-the-behaviour-of-experts-on-demand-forecasting-a-complex-task/
}
\author{
Nikolaos Kourentzes
}
\examples{
# Create some random MRE
mre <- runif(10,0,5) + runif(10,0,5)*1i
mre.plot(mre)
}
\keyword{ ~mre }
\seealso{
\code{\link{bias.coeff}}, \code{\link{mre}}.
}