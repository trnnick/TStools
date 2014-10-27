\name{init.ets}
\alias{init.ets}

\title{ETS with preset initial values.}

\description{This function overrides the ETS initial values. The model is not reoptimised and therefore should be used in conjuction with pre-defined smoothing parameters}

\usage{
init.ets(y,level=NULL,trend=NULL,season=NULL,...)
}

\arguments{
  \item{y}{
    A numeric vector or time series.
}
  \item{level}{
    Initial level.
}
  \item{trend}{
    Initial trend.
}
  \item{season}{
    Vector of initial seasonal indices.
}
  \item{...}{
    ets input arguments.
}
}
\value{
\item{fit}{An object of class "ets".}
}
\author{
Nikolaos Kourentzes
}
\examples{
init.ets(AirPassengers,trend=1.1)
}
\seealso{
\code{\link{ets}}.
}
\keyword{ ~ets }
