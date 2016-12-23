\name{geomean}
\alias{geomean}

\title{Geometric mean}

\description{Calculate geometric mean.}

\usage{
geomean(x,na.rm=c(FALSE,TRUE))
}

\arguments{
  \item{x}{
    Input data.
    }
  \item{na.rm}{
    a logical value indication whether NA values should be stripped before the computation proceeds.
    }
}
\value{
\item{gm}{Geometric mean result.}
}
\author{
Nikolaos Kourentzes
}
\examples{
geomean(c(0.7,1.1,1.3,1.2))
}
