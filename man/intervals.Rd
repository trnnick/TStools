\name{intervals}
\alias{intervals.sd}
\alias{intervals.mad}
\alias{intervals.mae}
\alias{intervals.ts}
\alias{intervals.hm}
\alias{intervals.qr}
\alias{intervals.qr2}
\alias{intervals.kde}
\alias{intervals.empir}
\alias{intervals.lowess}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Prediction intervals
}
\description{
    These function return prediction intervals based on the estimated "smooth" model.
}
\usage{
    intervals.sd(model, level=0.95, centre=TRUE)
    intervals.mae(model, level=0.95, centre=TRUE)
    intervals.mad(model, level=0.95, centre=TRUE)
    intervals.ts(model, level=0.95, centre=TRUE)
    intervals.hm(model, level=0.95, centre=TRUE)
    intervals.qr(model, level=0.95)
    intervals.qr2(model, level=0.95)
    intervals.kde(model, level=0.95, type=c("multi","single"),
                          kdetype=c("diffusion","SJ","nrd0"))
    intervals.empir(model, level=0.95, type=c("multi","single"), quantiletype=7)
    intervals.lowess(model, level=0.95, type=c("multi","single"), quantiletype=7)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{Object of the class 'smooth'.}
  \item{level}{Confidence level for the intervals construction.}
  \item{centre}{Defines whether the centering is needed for errors. If \code{FALSE}, then the errors are assumed to have zero mean.}
  \item{type}{Type of multi-step interval calculation: `multi' for using multiple steps ahead errors and `single' for using one step ahead error only.}
  \item{kdetype}{Type of kernel used for intervals.kde. See \code{\link{kdemode}} for details.}
  \item{quantiletype}{Quantile algorithm used. See \code{\link{quantile}} for details.}
}
\details{
    \itemize{
        \item \code{intervals.sd()} - intervals produced using sigma values for each horizon;
        \item \code{intervals.mae()} - intervals produced using MAE (Mean Absolute Error) for each horizon;
        \item \code{intervals.mad()} - intervals produced using MAD (Median Absolute Deviation) for each horizon;
        \item \code{intervals.ts()} - intervals produced using two sigma values;
        \item \code{intervals.hm()} - intervals produced using half moment;
        \item \code{intervals.qr()} - intervals produced using quantile regressions Taylor and Bunn (1999) style;
        \item \code{intervals.qr2()} - intervals produced using quantile regressions based on power function;
        \item \code{intervals.kde()} - intervals produced using kernel density estimation;
        \item \code{intervals.empir()} - intervals produced using empirical distribution;
        \item \code{intervals.lowess()} - intervals produced using lowess smoothing of edf;
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
    Ivan Svetunkov and Nikolaos Kourentzes
}

\seealso{
    Nothing
}
\examples{
    library(smooth)

    ourModel <- es(rnorm(100,100,10))
    intervals.sd(ourModel)
    intervals.mae(ourModel)
    intervals.mad(ourModel)
    intervals.ts(ourModel)
    intervals.hm(ourModel)
    intervals.qr(ourModel)
    intervals.qr2(ourModel)
    intervals.kde(ourModel)
    intervals.empir(ourModel)
    intervals.lowess(ourModel)
}
\keyword{ htest }
