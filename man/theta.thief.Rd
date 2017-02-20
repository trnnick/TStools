\name{theta.thief}
\alias{theta.thief}

\title{Theta method for THieF.}

\description{Function using the Theta method with Temporal Hierarchies.}

\usage{
theta.thief(y,h=NULL,...)
}
\arguments{
  \item{y}{
    Input time series. Can be ts object. 
    }
  \item{h}{
    Forecast horizon. If NULL then h is set to match frequency of time series.
    }
  \item{...}{
    Additional arguments passed to \code{\link{theta}}.
    }
}
\value{An object of class "\code{forecast}".
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
\note{This function is created to work with Temporal Hierarchied (\href{http://cran.r-project.org/package=thief}{thief} package). For conventional Theta method use \code{\link{theta}}.}
\author{
Nikolaos Kourentzes, \email{nikolaos@kourentzes.com}
}
\seealso{
\code{\link{theta}}.
}
\examples{
## **Not run:** 
library(thief)
frc <- thief(referrals,forecastfunction=theta.thief)
plot(frc)
## End(**Not run**)
}
\keyword{ ~thief }
\keyword{ ~ts }