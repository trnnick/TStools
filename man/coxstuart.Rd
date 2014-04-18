\name{coxstuart}
\alias{coxstuart}

\title{Cox-Stuart test}

\description{Cox-Stuart test for location and dispersion.}

\usage{
coxstuart(y,type=c("trend","deviation","dispersion"),alpha=0.05)
}

\arguments{
  \item{y}{
    A time series. Can be a time series object.
}
  \item{type}{
    Type of test: 1. trend - test for changes in trend [default]; 2. deviation - test for changes in deviation;3. dispersion - test for changes in dispersion (range).
}
  \item{alpha}{
    Significance level.
}
}
\value{
\item{H}{Hypothesis (H0/H1).}
\item{p.value}{P-value.}
\item{Htxt}{Description of the result.}
}
\author{
Nikolaos Kourentzes
}
\examples{
H <- coxstuart(referrals)
}

