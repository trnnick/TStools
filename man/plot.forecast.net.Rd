\name{plot.forecast.net}
\alias{plot.forecast.net}

\title{Plot neural network forecasts.}

\description{Plot MLP/ELM ensemble forecasts.}

\usage{
\method{plot}{forecast.net}(x, ...)
}

\arguments{
  \item{x}{
    An object of class \code{"forecast.net"}, produced using \code{\link{forecast.mlp}} or \code{\link{forecast.elm}}.
    }
  \item{...}{
    Unused input. 
    }
}
\value{None. Function produces a plot.}
\author{
Nikolaos Kourentzes, \email{nikolaos@kourentzes.com}
}
\seealso{
\code{\link{forecast.mlp}}, \code{\link{forecast.elm}}, \code{\link{mlp}}, \code{\link{elm}}.}
\examples{
fit <- mlp(AirPassengers)
print(fit)
plot(fit)
frc <- forecast(fit,h=36)
plot(frc)
}
\keyword{ ~mlp }
\keyword{ ~thief }
\keyword{ ~ts }