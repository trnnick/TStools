\name{plot.net}
\alias{plot.net}

\title{Plot MLP or ELM neural network.}

\description{Plot MLP or ELM neural network.}

\usage{
plot.net(fit,r=1)
}

\arguments{
  \item{fit}{
    Output of mlp or elm functions.
    }
  \item{r}{
    Training repetition to plot. 
    }
}
\author{
Nikolaos Kourentzes, \email{nikolaos@kourentzes.com}
}
\examples{
fit <- elm(AirPassengers)
plot.net(fit)
}
