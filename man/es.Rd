\name{es}
\alias{es}
\alias{sim.ets}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Exponential Smoothing in SSOE state-space model
}
\description{
This function is no longer available in TStools. It is a part of "smooth" package. Use \code{install.packages("smooth")} in order to install it.
}
\usage{
es(...)
sim.ets(...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{...}{ Any parameter you want to pass.
  }
}
\details{
    This function is just a place holder that gives you warning that function is no longer available.
    
    Note that sim.ets() has been renamed to sim.es().
}
\value{
    Only warnings are returned...
}
\references{
    \enumerate{
        \item{Hyndman, R.J., Koehler, A.B., Ord, J.K., and Snyder, R.D. (2008) Forecasting with exponential smoothing: the state space approach, Springer-Verlag. \url{http://www.exponentialsmoothing.net}.}
        \item{Taylor, J.W. and Bunn, D.W. (1999) A Quantile Regression Approach to Generating Prediction Intervals. Management Science, Vol 45, No 2, pp 225-237.}
        \item{Kolassa, S. (2011) Combining exponential smoothing forecasts using Akaike weights. International Journal of Forecasting, 27, pp 238 - 251.}
    }
}
\author{
Ivan Svetunkov, \email{ivan@svetunkov.ru}
}

\seealso{
\code{\link[forecast]{ets}, \link[forecast]{forecast}, \link[stats]{ts}, \link[smooth]{es}, \link[smooth]{sim.es}}
}
\examples{
library(Mcomp)

# See it for yourself!
\dontrun{es(M3$N1245$x,model="AAdN",h=8,holdout=FALSE,trace=FALSE)}
}
\keyword{ exponential smoothing }
\keyword{ ETS }
\keyword{ forecasting }
\keyword{ trace likelihood }
