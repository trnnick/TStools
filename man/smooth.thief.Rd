\name{smooth.thief}
\alias{smooth.thief}

\title{Smooth models for THieF.}

\description{Function produces forecasts with Temporal Hierarchies using models from the smooth package.}

\usage{
smooth.thief(data,h,type=c("es","ssarima","ces","ges"),...)
}
\arguments{
  \item{data}{
    Input time series. Can be ts object. 
    }
  \item{h}{
    Forecast horizon. 
    }
  \item{type}{
    Can be "es" to use Exponential Smoothing, "ssarima" to use state-space ARIMA, "ces" to use Complex Exponential Smoothing, "ges" to use Generalised Exponential Smoothing. 
  }
  \item{...}{
    Additional arguments passed forecasting model.
    }
}
\value{An object of class "\code{forecast}".
The function \code{plot} produces a plot of the forecasts.
It has a list containing the following elements:
\item{model}{The estimated model.}
\item{fitted}{Fitted values of the model.}
\item{x}{Actuals provided in the call of the model.}
\item{mean}{Point forecasts of the model.}
}
\note{This function is created to work with Temporal Hierarchies (\href{https://cran.r-project.org/package=thief}{thief} package).}
\author{
Nikolaos Kourentzes, \email{nikolaos@kourentzes.com}
}
\seealso{
\code{\link{es}}, \code{\link{auto.ssarima}}, \code{\link{auto.ces}}
}
\examples{
## **Not run:** 
library(thief)
frc <- thief(referrals,forecastfunction=smooth.thief,type="es")
plot(frc)
## End(**Not run**)
}
\keyword{ ~thief }
\keyword{ ~ts }