\name{theta}
\alias{theta}

\title{Theta method}

\description{Theta method for time series forecasting.}

\usage{
theta(y,m=NULL,h=10,outplot=0,sign.level=0.05,
      cost0=c("MSE","MdSE","MAE","MdAE"),
      cost2=c("MSE","MdSE","MAE","MdAE"),
      costs=c("MSE","MdSE","MAE","MdAE"),
      multiplicative=c(TRUE,FALSE),cma=NULL)
}

\arguments{
  \item{y}{
    Time series to model. Can be either a vector or a ts object.
}
  \item{m}{
    Periods in a season of the time series. If insample is a ts object then this is taken from its frequency, unless overriden.
}
  \item{h}{
    Forecast horizon. Default is 10.
}
  \item{outplot}{
    Provide plot: 0 - No plot; 1 - Series and forecast; 2 - As above with theta lines.
}
  \item{sign.level}{
    Significance level for trend and seasoanlity statistical tests.
}
  \item{cost0}{
    Cost function of theta0 line. Can be: "MSE", "MdSE", "MAE", "MdAE".
}
  \item{cost2}{
    Cost function of theta2 line. Can be: "MSE", "MdSE", "MAE", "MdAE".
}
  \item{costs}{
    Cost function of seasonal element. Can be: "MSE", "MdSE", "MAE", "MdAE".
}
  \item{multiplicative}{
    If TRUE then multiplicative decomposition is performed. Otherwise additive is used.
}
  \item{cma}{
    Input pre-calculated centred moving average. Use NULL to calculate internally. 
}
}
\value{
\item{frc}{Forecasts.}
\item{exist}{exist[1] is the result for trend, exist[2] is for season.}
\item{theta0}{Forecasted values of theta0 line.}
\item{theta2}{Forecasted values of theta2 line.}
\item{season}{Forecasted values of seasonal element.}
\item{a}{SES parameters of theta2.}
\item{b}{Regression parameters of theta0.}
\item{g}{Pure seasonal exponential smoothing parameters of season.}
}
\note{This implementation of Theta method tests automatically for seasonality and trend. Seasonal decomposition can be done either additively or multiplicatively and the seasonality is treated as a pure seasonal model. The various components of Theta can be optimised using different cost functions. The originally proposed Theta method always assume multiplicative seasonality and presence of trend, while all theta lines are optimised using MSE. Seasonality is estimated using classical decomposition.}
\author{
Nikolaos Kourentzes
}
\examples{
theta(referrals,outplot=2)
}