\name{plot.net}
\alias{plot.net}

\title{Plot MLP or ELM neural network.}

\description{Plot MLP or ELM neural network.}

\usage{
    \method{plot}{net}(x, r=1, ...)
}

\arguments{
  \item{x}{
    Output of mlp or elm functions.
    }
  \item{r}{
    Training repetition to plot. 
    }
  \item{...}{
    Parameters passed to plot.
    }
}
\author{
Nikolaos Kourentzes, \email{nikolaos@kourentzes.com}
}
\examples{
fit <- elm(AirPassengers)
plot.net(fit)
}
