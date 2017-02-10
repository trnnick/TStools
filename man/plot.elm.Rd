\name{plot.elm}
\alias{plot.elm}

\title{Plot ELM network.}

\description{Produces a plot of the ELM network architecture.}

\usage{
\method{plot}{elm}(x, r=1, ...)
}
\arguments{
  \item{x}{
    ELM network object, produced using \code{\link{elm}}.
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
fit <- elm(AirPassengers)
print(fit)
plot(fit)
frc <- forecast(fit,h=36)
plot(frc)
}
\keyword{ ~elm }
\keyword{ ~thief }
\keyword{ ~ts }