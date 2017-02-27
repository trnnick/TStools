\name{xyz}
\alias{xyz}

\title{XYZ analysis}

\description{Perform XYZ analysis based on coefficient of variation, naive and seasonal naive fit, exponential smoothing or custom quantities.}

\usage{
xyz(x,m=NULL,prc=c(0.2,0.3,0.5),type=c("naive","ets","cv"))
}

\arguments{
  \item{x}{
    This can either be an array, where each column is a series, or a vector of values. If x is an array of time series, then the importance of each series is calculated as the mean of the observations (sales volume). Otherwise, value can be whatever quantity is provided. 
  }
  \item{m}{
    Seasonal length for time series. Required when type is "naive" or "ets".
  }
  \item{prc}{
    A list of percentages to separate the items in. By default this is c(0.2,0.3,0.5), but any set of percentage values can be used as long as 0<=prc[i]<=1 and sum(prc)==1.
  }
  \item{type}{
    The type of forecastability calculation: 1: "naive" - fit naive and seasonal naive and calculate fit RMSE/mean level; 2: "ets" - fit ets and calculate fit RMSE/mean level; 3: "cv" - use coefficient of variation as a proxy of forecastability
  }
}
\value{An object of class "\code{abc}". 
The function \code{plot} produces a plot of the analysis.
An object of class \code{"abc"} is a list containing the following elements:
\item{type}{The type of analysis.}
\item{prc}{The requested percentage for each class.}
\item{value}{A vector containing the importance value of each series.}
\item{class}{A vector containing the class membership of each series.}
\item{rank}{A vector containing the rank of each series. Highest value is 1.}
\item{conc}{The error concentration of each class, as percentage of total value.}
\item{model}{Fitted model for each series.}
}
\author{
Nikolaos Kourentzes, \email{nikolaos@kourentzes.com}
}
\seealso{
\code{\link{abc}}, \code{\link{plot.abc}}, \code{\link{abcxyz}}.
}
\examples{
x <- abs(matrix(cumsum(rnorm(5400,0,1)),36,150))
xyz(x,m=12)
}