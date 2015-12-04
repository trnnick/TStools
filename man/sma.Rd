\name{sma}
\alias{sma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simple Moving Average
}
\description{
Function constructs moving average for the provided data and produces forecast.
}
\usage{
sma(x, order=NULL, h=10, optimize=FALSE,
    min.order=1, max.order=12, holdout=FALSE,
    silent=FALSE, legend=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ The vector or ts object, containing the data needed to be forecasted.
  }
  \item{order}{ The order of the moving average.
  }
  \item{h}{ The length of the forecasting horizon.
  }
  \item{optimize}{ Binary that determines if the order of moving average should be defined automatically.
  }
  \item{min.order}{ The minimum order, that should be used in the optimisation.
  }
  \item{max.order}{ The maximum order, that should be used in the optimisation.
  }
  \item{holdout}{ If \code{TRUE}, the holdout sample of size \code{h} is taken from the end of the data.
  }
  \item{silent}{ If \code{TRUE}, no output or graphs are printed out.
  }
  \item{legend}{ If \code{TRUE}, the legend is drawn in the outplot.
  }
}
\value{
  The list of the following values is returned:

  \itemize{
    \item \code{order} - the order of SMA used in the construction.
    \item \code{fitted} - the fitted values.
    \item \code{forecast} - the point forecast of SMA.
  }
}
\references{
    \enumerate{
        \item{Hyndman, R.J., Athanasopoulos G. (2014) Forecasting: principles and practice. \url{https://www.otexts.org/fpp/6/2}.}
    }
}
\author{
Patrick Saoud, \email{p.saoud@lancaster.ac.uk}
}

\seealso{
\code{\link[forecast]{ma}, \link[TStools]{es}}
}
\examples{
library(Mcomp)

sma(M3$N1812$x,order=8,h=18,holdout=TRUE) -> test
sma(M3$N1812$x,h=18,optimize=TRUE,holdout=TRUE) -> test
}
\keyword{ simple moving average }
\keyword{ sma }
\keyword{ forecasting }