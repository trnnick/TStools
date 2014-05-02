\name{decomp}
\alias{decomp}

\title{Decomposition of series}

\description{Decompose series into its components.}

\usage{
decomp(y,m=NULL,s=NULL,trend=NULL,outplot=c(FALSE,TRUE),
       decomposition=c("multiplicative","additive"),
       h=0,type=c("mean","median","pure.seasonal"),w=NULL)
}

\arguments{
  \item{y}{
    Vector of target data, can be ts object.
}
  \item{m}{
    Seasonal period. If y is a ts object then the default is its frequency.
}
  \item{s}{
    Starting period in the season. If y is a ts object then default is read.
}
  \item{trend}{
   Vector of level/trend of the time series. If NULL then it is calculated using CMA.
}
  \item{outplot}{
    If TRUE provide a plot of the decomposed components.
}
  \item{decomposition}{
    Type of seasonal decomposition: "multiplicative" or "additive".
}
  \item{h}{
    Forecast horizon for seasonal component.
}
  \item{type}{
    Type of calculation for seasonal component: "mean" - The mean of each seasonal period; "median" - The median of each seasonal period; "pure.seasonal" - Model using a pure seasonal model.
}
  \item{w}{
    Percentage or number of observations to winsorise in the calculation of mean seasonal indices. If w>1 then it is the number of observations, otherwise it is a percentage. If type != "mean" then this is ignored.
}
}
\value{
\item{trend}{Trend component.}
\item{season}{Seasonal component.}
\item{irregular}{Irregular component.}
\item{f.season}{Forecasted seasonal component if h>0.}
\item{g}{Purse seasonal model parameters.}
}

\author{
Nikolaos Kourentzes
}
\examples{
decomp(referrals,outplot=TRUE)
}