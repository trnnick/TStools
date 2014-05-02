\name{wins}
\alias{wins}

\title{Winsorise}

\description{Winsorise vector to either a percentage or number of observations.}

\usage{
wins(x, p=0.05)
}

\arguments{
  \item{x}{
    Vector to be winsorised.
}
  \item{p}{
    Percentage or number of observations to be winsorised. If value is <1 then it is used as a percentages. Otheriwse it is the number of observations to winsorise. If the resulting p > floor((length(x)-1)/2) then it is set equal to floor((length(x)-1)/2).
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