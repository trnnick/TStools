\name{theta}
\alias{theta}

\title{Theta method}

\description{Theta method for time series forecasting.}

\usage{
theta(y,m=NULL,sign.level=0.05,
      cost0=c("MSE","MdSE","MAE","MdAE"),
      cost2=c("MSE","MdSE","MAE","MdAE"),
      costs=c("MSE","MdSE","MAE","MdAE"),
      multiplicative=c(TRUE,FALSE),cma=NULL,
      outliers=NULL)
}

\arguments{
  \item{y}{
    Time series to model. Can be either a vector or a ts object.
}
  \item{m}{
    Periods in a season of the time series. If y is a ts object then this is taken from its frequency, unless overriden.
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
  \item{outliers}{
    Optional. Provide vector of location of observations that are considered outliers. These will be included in theta0 estimation. To consider no outliers then use NULL.
}
}
\value{An object of class "\code{theta}".
The function \code{plot} produces a plot of the in-sample fit.
An object of class \code{"theta"} is a list containing the following elements:
\item{method}{Forecasting method.}
\item{y}{Input time series.}
\item{m}{Periods in a season of the time series.}
\item{exist}{exist[1] is the result for trend, exist[2] is for season.}
\item{multiplicative}{If TRUE seasonality is modelled multiplicatively.}
\item{theta0}{Fitted theta0 line values.}
\item{theta2}{Fitted theta2 line values.}
\item{season}{Seasonal profile.}
\item{a}{SES parameters of theta2.}
\item{b}{Regression parameters of theta0.}
\item{p}{Coefficients of outliers from theta0, theta2 and seasonality estimation.}
\item{g}{Pure seasonal exponential smoothing parameters of season.}
\item{fitted}{Fitted values.}
\item{residuals}{In-sample residuals.}
\item{MSE}{In-sample Mean Squared Error.}
}
\note{This implementation of Theta method tests automatically for seasonality and trend. Seasonal decomposition can be done either additively or multiplicatively and the seasonality is treated as a pure seasonal model. The various components of Theta can be optimised using different cost functions. The originally proposed Theta method always assumes multiplicative seasonality and presence of trend, while all theta lines are optimised using MSE. Seasonality is estimated using classical decomposition.
To use Theta with Temporal Hierarchies (\href{http://cran.r-project.org/package=thief}{thief} package) see \code{\link{theta.thief}}.
}
\author{
Nikolaos Kourentzes, \email{nikolaos@kourentzes.com}
}
\seealso{
\code{\link{forecast.theta}}, \code{\link{plot.theta}}, \code{\link{theta.thief}}.
}
\examples{
theta(referrals)
}
\keyword{ ~ts }