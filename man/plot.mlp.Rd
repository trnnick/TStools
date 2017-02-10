\name{plot.mlp}
\alias{plot.mlp}

\title{Plot MLP network.}

\description{Produces a plot of the MLP network architecture.}

\usage{
\method{plot}{mlp}(x, r=1, ...)
}
\arguments{
  \item{x}{
    MLP network object, produced using \code{\link{mlp}}.
    }
  \item{r}{
    Ensemple member to plot.
    }
  \item{...}{
    Unused argument. 
  }
}
\value{None. Function produces a plot.}
\note{Neurons are coloured with \code{"lightgrey"}. Seasonal dummies are coloured with \code{"lightpink"} and xreg with \code{"lightblue"}.}
\author{
Nikolaos Kourentzes, \email{nikolaos@kourentzes.com}
}
\seealso{
\code{\link{elm}}, \code{\link{mlp}}.
}
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