\name{intervals.twoSigmas}
\alias{intervals.twoSigmas}
\alias{intervals.hm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Prediction intervals based on two standard deviations
}
\description{
    This function returns prediction intervals based on two separate sigmas for left and right sides of the distribution
}
\usage{
    intervals.twoSigmas(model, level=0.95)
    intervals.hm(model, level=0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{Object of the class 'smooth'.}
  \item{level}{Confidence level for the intervals construction.}
}
\details{
    Nothing here yet.
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
    intervals.twoSigmas(ourModel)
    intervals.hm(ourModel)
}
\keyword{ htest }
