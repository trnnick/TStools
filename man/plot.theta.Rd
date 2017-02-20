\name{plot.theta}
\alias{plot.theta}

\title{Plot Theta method fit.}

\description{PLot Theta method fit and Theta lines.}

\usage{
\method{plot}{theta}(x, thetalines=c(TRUE,FALSE), ...)
}
\arguments{
  \item{x}{
    Theta method object, produced using \code{\link{theta}}.
    }
  \item{thetalines}{
    If TRUE then Theta0 and Theta2 lines are plotted.
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
\code{\link{theta}}, \code{\link{forecast.theta}}, \code{\link{theta.thief}}.
}
\examples{
fit <- theta(referrals)
plot(fit)
}