\name{seasplot}
\alias{seasplot}

\title{Seasonal plots}

\description{Seasonal plots and crude trend/season tests.}

\usage{
seasplot(y,m=NULL,s=NULL,trend=NULL,colour=NULL,alpha=0.05,outplot=c(1,0,2,3,4,5),
         decomposition=c("multiplicative","additive"),cma=NULL,labels=NULL)
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
    If TRUE then a trend is assumed and is removed using CMA. If FALSE then no trend is assumed. If NULL then trend is identified and removed if found.
}
  \item{colour}{
    Single colour override for plots.
}
  \item{alpha}{
    Significance level for statistical tests (kpss and friedman).
}
  \item{outplot}{
    Provide plot output: 0 - None; 1 - Seasonal diagramme; 2 - Seasonal boxplots; 3 - Seasonal subseries; 4 - Seasonal distribution; 5 - Seasonal density.
}
  \item{decomposition}{
    Type of seasonal decomposition: "multiplicative" or "additive".
}
  \item{cma}{
    Input precalculated level/trend for the analysis. Overrides trend=NULL.
}
  \item{labels}{
    External labels for the seasonal periods. Use NULL for default. If length(labels) < m, then this input is ignored.
}
}
\value{
\item{season}{Matrix of (detrended) seasonal elements.}
\item{season.exist}{TRUE/FALSE results of friedman test.}
\item{season.pval}{Friedman test p-value.}
\item{trend}{CMA estimate or NULL if trend == FALSE.}
\item{trend.exist}{TRUE/FALSE result of kpss test.}
\item{trend.pval}{kpss (approximate) p-value.}
}
\author{
Nikolaos Kourentzes
}
\examples{
seasplot(referrals)
}