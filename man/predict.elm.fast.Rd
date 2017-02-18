\name{predict.elm.fast}
\alias{predict.elm.fast}

\title{Predictions for ELM (fast) network.}

\description{Calculate predictions for ELM (fast) network.}

\usage{
\method{predict}{elm.fast}(object, newx, na.rm=c(FALSE,TRUE),...)
}
\arguments{
  \item{object}{
    ELM network object, produced using \code{\link{elm.fast}}.
    }
  \item{newx}{
    Explanatory variables. Each column is a variable. 
    }
  \item{na.rm}{
    If TRUE remove columns and object produces an ensemble forecast, then remove any members that give NA in their forecasts.
  }
  \item{...}{
    Unused argument. 
  }
}
\value{
\item{Y.hat}{Ensemble prediction.}
\item{Y.all}{Predictions of each training repetition.}
}
\author{
Nikolaos Kourentzes, \email{nikolaos@kourentzes.com}
}
\seealso{
\code{\link{elm.fast}}, \code{\link{plot.elm.fast}}.
}
\examples{
p <- 2000
n <- 150
X <- matrix(rnorm(p*n),nrow=n)
b <- cbind(rnorm(p))
Y <- X \%*\% b
fit <- elm.fast(Y,X)
predict(fit,X)
}
\keyword{ ~elm }