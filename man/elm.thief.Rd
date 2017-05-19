\name{elm.thief}
\alias{elm.thief}

\title{ELM network for THieF.}

\description{Function for ELM forecasting with Temporal Hierarchies.}

\usage{
elm.thief(y,h=NULL,...)
}
\arguments{
  \item{y}{
    Input time series. Can be ts or msts object. 
    }
  \item{h}{
    Forecast horizon. If NULL then h is set to match frequency of time series.
    }
  \item{...}{
    Additional arguments passed to \code{\link{elm}}.
    }
}
\value{An object of classes "\code{forecast.net}" and "\code{forecast}".
The function \code{plot} produces a plot of the forecasts.
An object of class \code{"forecast.net"} is a list containing the following elements:
\item{method}{The name of the forecasting method as a character string}
\item{mean}{Point forecasts as a time series}
\item{all.mean}{An array h x reps of all ensemble members forecasts, where reps are the number of ensemble members.}
\item{x}{The original time series (either \code{fit} used to create the network.}
\item{fitted}{Fitted values. Any values not fitted for the initial period of the time series are imputted with NA.}
\item{residuals}{Residuals from the fitted network.}
}
\references{
    \itemize{
        \item{For combination operators see: Kourentzes N., Barrow B.K., Crone S.F. (2014) \href{http://kourentzes.com/forecasting/2014/04/19/neural-network-ensemble-operators-for-time-series-forecasting/}{Neural network ensemble operators for time series forecasting}. \emph{Expert Systems with Applications}, \bold{41}(\bold{9}), 4235-4244.}
        \item{For ELMs see: Huang G.B., Zhou H., Ding X. (2006) Extreme learning machine: theory and applications. \emph{Neurocomputing}, \bold{70}(\bold{1}), 489-501.}
}
}
\note{This function is created to work with Temporal Hierarchied (\href{https://cran.r-project.org/package=thief}{thief} package). For conventional ELM networks use \code{\link{elm}}.}
\author{
Nikolaos Kourentzes, \email{nikolaos@kourentzes.com}
}
\seealso{
\code{\link{elm}}, \code{\link{plot.forecast.net}}.
}
\examples{
## **Not run:** 
library(thief)
frc <- thief(AirPassengers,forecastfunction=elm.thief)
plot(frc)
## End(**Not run**)
}
\keyword{ ~elm }
\keyword{ ~thief }
\keyword{ ~ts }