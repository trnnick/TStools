\name{mre}
\alias{mre}

\title{Calculate Root Error}

\description{Given k errors mre calculates the Mean Root Error, the Sum Root Error or the Mod(Geometric Squared Mean Root Error) = GRMSE.}

\usage{
mre(e,op=c("mean","sum","gm"))
}

\arguments{
  \item{e}{
    k forecast error(s).
}
  \item{op}{
    Aggregation operator of the Root Error: 1) "mean" - Mean; 2) "sum" - Sum; 3) "gm" - Magnitude of Geometric Mean
}
}
\value{
\item{mre}{Root Error.}
}
\references{
Kourentzes N., Trapero J. R., Svetunkov I. Measuring the behaviour of experts on demand forecasting: a complex task. Working paper
http://kourentzes.com/forecasting/2014/12/17/measuring-the-behaviour-of-experts-on-demand-forecasting-a-complex-task/
}
\author{
Nikolaos Kourentzes
}
\examples{
# Generate some random errors
e <- runif(10,-10,10)
mre(e)
}
\keyword{ ~mre }
\seealso{
\code{\link{bias.coeff}}, \code{\link{mre.plot}}.
}