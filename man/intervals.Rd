\name{intervals}
\alias{intervals}
\title{
    Construction of prediction intervals
}
\description{
    Function returns quantiles of distribution using different methods based on the provided vector or matrix of random variables.
}
\usage{
intervals(data,intType=c("standard","2sd","hm"),centre=NULL,
          level=0.95,df=NULL,k=NULL)
}
\arguments{
    \item{data}{ Vector or matrix of random variable. Can be, for example, a matrix of 1 to h steps ahead errors.
    }
    \item{intType}{ Type of intervals to construct. The first letter can be used instead of the whole word.
    }
    \item{centre}{ Value from which to calculate statistics. By default \code{centre=mean(data)}.
    }
    \item{level}{ Defines width of the prediction interval.
    }
    \item{df}{ Number of degrees of freedom. If not provided, then set to n-1.
    }
    \item{k}{ Number of parameters used in model. Needed if \code{df} is NULL.
    }
}
\details{
    If \code{data} is matrix, then several intervals are constructed (for each column). Otherwise only one is done.
}
\value{
  Function returns:

  \itemize{
    \item \code{lower} - lower bound of prediction interval.
    \item \code{upper} - upper bound of prediction interval.
  }
}
\author{
    Ivan Svetunkov
}
\examples{
x <- rpois(100,1.5)
intervals(x,intType="hm",centre=median(x))

x <- matrix(exp(rnorm(1000,0,1)),100,10)
intervals(x,intType="2sd",centre=apply(x,2,median))
}
\keyword{ intervals }
\keyword{ prediction }
\keyword{ prediction intervals }
