\name{wins}
\alias{wins}

\title{Winsorise vector}

\description{Winsorise vector to either a percentage or number of observations.}

\usage{
wins(x, prc=0.05, n=NULL)
}

\arguments{
  \item{x}{
    Vector to be winsorised.
}
  \item{prc}{
    Percentage of observations to be winsorised. Values about 0.5 are overriden to 0.5.
}
  \item{n}{
    Number of observations to be winsorised. If given this overrides prc. If n > floor((length(x)-1)/2) then it is set equal to floor((length(x)-1)/2).
}
}
\value{
\item{x.out}{Winsorised vector.}
}
\author{
Nikolaos Kourentzes
}
\examples{
x <- rnorm(100,mean=0,sd=1)
x.out <- wins(x)
}