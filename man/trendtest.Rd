\name{trendtest}
\alias{trendtest}

\title{Test series for the presence of trend}

\description{Test a time series for trend uing AICc of ETS models or Cox-Struart, using Multiple Temporal Aggregation or not.}

\usage{
trendtest(y,extract=c("FALSE","TRUE"),type=c("aicc","cs"),mta=c(FALSE,TRUE))
}

\arguments{
  \item{y}{
    A time series. This must be of class "ts".
}
  \item{extract}{
    IF TRUE, then use Centred Moving Average to extract the trend of the series and test that.
}
  \item{type}{
    Type of trend test. Use "aicc" to make the choice depending on the AICc of ETS models. Use "cs" to use the Cox-Stuart test at 5\% significance level. 
}
  \item{mta}{
    If TRUE then Multiple Temporal Aggregation is used and the tests are performed at multiple views of the series. The result is the majority vote. 
}
}
\value{
\item{H}{Result of the test: TRUE there is evidence of trend, FALSE otherwise.}
}
\author{
Nikolaos Kourentzes
}
\examples{
trendtest(AirPassengers,extract=TRUE)
}