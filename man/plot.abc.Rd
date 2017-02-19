\name{plot.abc}
\alias{plot.abc}

\title{Plot ABC/XYZ analysis.}

\description{Produces a plot of the ABC or XYZ analysis.}

\usage{
\method{plot}{abc}(x, cex.prc=0.8,...)
}
\arguments{
  \item{x}{
    \code{"abc"} class object, produced using \code{\link{abc}} or \code{\link{xyz}}.
    }
  \item{cex.prc}{
    Font size of percentage reported in plot.
  }
  \item{...}{
    Additional arguments can be passed for the plot. User can override default values. For example main="" will remove the title from the plot.
  }
}
\value{None. Function produces a plot.}
\author{
Nikolaos Kourentzes, \email{nikolaos@kourentzes.com}
}
\seealso{
\code{\link{abc}}, \code{\link{xyz}}, \code{\link{abcxyz}}.
}
\examples{
x <- abs(matrix(cumsum(rnorm(5400,0,1)),36,150))
z <- abc(x)
plot(z)
}
\keyword{ ~ts }