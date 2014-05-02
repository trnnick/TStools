\name{colWins}
\alias{colWins}

\title{Winsorise per column}

\description{Winsorise per column to either a percentage or number of observations.}

\usage{
colWins(x, prc=0.05, n=NULL)
}

\arguments{
  \item{x}{
    Array/matrix to be winsorised.
}
  \item{prc}{
    Percentage of observations to be winsorised per column. Values about 0.5 are overriden to 0.5.
}
  \item{n}{
    Number of observations to be winsorised per column. If given this overrides prc. If n > floor((length(x)-1)/2) then it is set equal to floor((length(x)-1)/2).
}
}
\value{
\item{x.out}{Winsorised array/matrix.}
}
\author{
Nikolaos Kourentzes
}
\examples{
x <- matrix(rnorm(10*20,mean=0,sd=1),10,20) 
x.out <- colWins(x)
}