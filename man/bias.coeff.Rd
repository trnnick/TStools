\name{bias.coeff}
\alias{bias.coeff}

\title{Bias Coefficient Calculation}

\description{Calculate the bias coefficient and plot it.}

\usage{
bias.coeff(mre,outplot=c(0,1,2))
}

\arguments{
  \item{mre}{
    k Root Errors.
}
  \item{outplot}{
    Optional plot: 0 - No; 1 - Histogram; 2 - Boxplot.
}
}
\value{
\item{bias}{k Bias Coefficients.}
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
bias.coeff(mre,outplot=2)
}
\keyword{ ~mre }
\seealso{
\code{\link{mre}}, \code{\link{mre.plot}}.
}