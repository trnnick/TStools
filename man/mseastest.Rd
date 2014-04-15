\name{mseastest}
\alias{mseastest}

\title{Seasonal plots}

\description{Multiplicative seasonality test.}

\usage{
mseastest(y,m=NULL,type=c("pearson","spearman","kendall"),cma=NULL,
          sn=1, alpha=0.05,outplot=c(FALSE,TRUE))
}

\arguments{
  \item{y}{
    Vector of target data, can be ts object.
}
  \item{m}{
    Seasonal period. If y is a ts object then the default is its frequency.
}
  \item{type}{
    Test is based on: 1. "pearson" correlation [Default]; 2. "spearman" correlation; 3. "kendall" correlation.
}
  \item{cma}{
    Level/trend of time series. If not given a central moving average is calculated.
}
  \item{sn}{
    Seasonal periods of decreasing magnitude to consider for the test. 
    The suggested value is 1. Default = 1.
}
  \item{alpha}{
    Significance level. Default = 0.05.
}
  \item{outplot}{
    If TRUE produce a plot of seasonal elements and level.
}
}
\value{
\item{out$is.multiplicative}{If TRUE the test found evidence of multiplicative seasonality.}
\item{out$statistic}{The test statistic.}
\item{out$pvalue}{P-value of the test.}
}
\author{
Nikolaos Kourentzes
}
\examples{
mseastest(referrals)
}