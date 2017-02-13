\name{plot.elm.fast}
\alias{plot.elm.fast}

\title{Plot ELM (fast) network.}

\description{Produces a scatter plot of the ELM fit.}

\usage{
\method{plot}{elm.fast}(x, ...)
}
\arguments{
  \item{x}{
    ELM network object, produced using \code{\link{elm.fast}}.
    }
  \item{...}{
    Unused argument. 
  }
}
\value{None. Function produces a plot.}
\author{
Nikolaos Kourentzes, \email{nikolaos@kourentzes.com}
}
\seealso{
\code{\link{elm.fast}}, \code{\link{predict.elm.fast}}.
}
\examples{
p <- 2000
n <- 150
X <- matrix(rnorm(p*n),nrow=n)
b <- cbind(rnorm(p))
Y <- X \%*\% b
fit <- elm.fast(Y,X)
plot(fit)
}
\keyword{ ~elm }