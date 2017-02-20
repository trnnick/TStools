\name{forecast.elm}
\alias{forecast.elm}

\title{Forecast using ELM neural network.}

\description{Create forecasts using ELM neural networks.}

\usage{
\method{forecast}{elm}(object,h=NULL,y=NULL,xreg=NULL,...)
}
\arguments{
  \item{object}{
    ELM network object, produced using \code{\link{elm}}.
    }
  \item{h}{
    Forecast horizon. If NULL then h is set to match frequency of time series.
    }
  \item{y}{
    Optionally forecast using different data than what the network was trained on. Expected to create havoc and do really bad things!
    }
   \item{xreg}{
    Exogenous regressors. Each column is a different regressor and the sample size must be at least as long as the target in-sample set plus the forecast horizon, but can be longer. Set it to NULL if no xreg inputs are used.
    }
  \item{...}{
    Unused argument. 
  }
}
\value{
An object of classes "\code{forecast.net}" and "\code{forecast}".
The function \code{plot} produces a plot of the forecasts.
An object of class \code{"forecast.net"} is a list containing the following elements:
\item{method}{The name of the forecasting method as a character string}
\item{mean}{Point forecasts as a time series}
\item{all.mean}{An array h x reps of all ensemble members forecasts, where reps are the number of ensemble members.}
\item{x}{The original time series used to create the network.}
\item{fitted}{Fitted values.}
\item{residuals}{Residuals from the fitted network.}
}
\author{
Nikolaos Kourentzes, \email{nikolaos@kourentzes.com}
}
\seealso{
\code{\link{elm}}, \code{\link{plot.forecast.net}}, \code{\link{elm.thief}}, \code{\link{mlp}}.
}
\examples{
fit <- elm(AirPassengers)
print(fit)
plot(fit)
frc <- forecast(fit,h=36)
plot(frc)
}
\keyword{ ~elm }
\keyword{ ~thief }
\keyword{ ~ts }