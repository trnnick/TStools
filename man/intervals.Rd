\name{intervals.ts}
\alias{intervals.ts}
\alias{intervals.hm}
\alias{intervals.qr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Prediction intervals
}
\description{
    These function return prediction intervals based on the estimated "smooth" model.
}
\usage{
    intervals.ts(model, level=0.95)
    intervals.hm(model, level=0.95)
    intervals.qr(model, level=0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{Object of the class 'smooth'.}
  \item{level}{Confidence level for the intervals construction.}
}
\details{
    \itemize{
        \item \code{intervals.ts()} - intervals produced using two sigma values;
        \item \code{intervals.hm()} - intervals produced using half moment;
        \item \code{intervals.qr()} - intervals produced using quantile regressions;
    }
}
\value{
  \item{lower}{Lower bound of prediction interval.}
  \item{upper}{Upper bound of prediction interval.}
  \item{sigma}{Matrix with lower and upper sigmas.}
  \item{hmValues}{Vector with values of half moments for each horizon.}
}
\references{
    Nothing here yet.
}
\author{
    Ivan Svetunkov
}

\seealso{
    Nothing
}
\examples{
    library(smooth)

    ourModel <- es(rnorm(100,100,10))
    intervals.ts(ourModel)
    intervals.hm(ourModel)
    intervals.qr(ourModel)
}
\keyword{ htest }
