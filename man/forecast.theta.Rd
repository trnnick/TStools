\name{forecast.theta}
\alias{forecast.theta}

\title{Forecast using Theta method.}

\description{Create forecasts using ELM neural networks.}

\usage{
\method{forecast}{theta}(object,h=NULL,...)
}
\arguments{
  \item{object}{
    Theta method network object, produced using \code{\link{theta}}.
    }
  \item{h}{
    Forecast horizon. If NULL then h is set to match frequency of time series.
    }
  \item{...}{
    Unused argument. 
  }
}
\value{
An object of class "\code{forecast}".
The function \code{plot} produces a plot of the forecasts.
It has a list containing the following elements:
\item{method}{The name of the forecasting method as a character string}
\item{mean}{Point forecasts as a time series}
\item{frc.theta0}{Point forecasts for Theta0 line.}
\item{frc.theta2}{Point forecasts for Theta2 line.}
\item{frc.season}{Point forecasts for seasonal profile.}
\item{x}{The original time series used to fit the method.}
\item{fitted}{Fitted values.}
\item{residuals}{Residuals.}
}
\author{
Nikolaos Kourentzes, \email{nikolaos@kourentzes.com}
}
\seealso{
\code{\link{theta}}, \code{\link{plot.theta}}, \code{\link{theta.thief}}.
}
\examples{
fit <- theta(referrals)
frc <- forecast(fit)
plot(frc)
}
\keyword{ ~ts }