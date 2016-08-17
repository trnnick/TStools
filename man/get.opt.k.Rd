\name{get.opt.k}
\alias{get.opt.k}

\title{Optimal temporal aggregation level}

\description{Find optimal temporal aggregation level for AR(1), MA(1), ARMA(1,1).}

\usage{
get.opt.k(y,m=12,type=c("ar","ma","arma"))
}

\arguments{
  \item{y}{
    Time series (ts object).
}
  \item{m}{
    Maximum temporal aggregation level to consider.
}
  \item{type}{
    Type of DGP which can be "ar" for AR(1), "ma" for MA(1) or "arma" for ARMA(1,1).
}
}
\value{
\item{k}{Identified optimal temporal aggregation level.}
}
\references{
Rostami-Tabar, Bahman, et al. "Demand forecasting by temporal aggregation." Naval Research Logistics (NRL) 60.6 (2013): 479-498.

Rostami-Tabar, Bahman, et al. "A note on the forecast performance of temporal aggregation." Naval Research Logistics (NRL) 61.7 (2014): 489-500.
}
\author{
Nikolaos Kourentzes & Bahman Rostami-Tabar
}
\examples{
get.opt.k(AirPassengers,12)
}