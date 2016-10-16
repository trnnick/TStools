\name{abc}
\alias{abc}

\title{ABC analysis}

\description{Perform ABC analysis based on volumes or custom quantities.}

\usage{
abc(x,prc=c(0.2,0.3,0.5),outplot=c(TRUE,FALSE),cex.prc=0.8,...)
}

\arguments{
  \item{x}{
    This can either be an array, where each column is a series, or a vector of values. If x is an array of time series, then the importance of each series is calculated as the mean of the observations (sales volume). Otherwise, value can be whatever quantity is provided. 
  }
  \item{prc}{
    A list of percentages to separate the items in. By default this is c(0.2,0.3,0.5), but any set of percentage values can be used as long as 0<=prc[i]<=1 and sum(prc)==1.
  }
  \item{outplot}{
    If TRUE provide a visualisation of the ABC analysis result.
  }
  \item{cex.prc}{
    Font size of percentage reported in plot.
  }
  \item{...}{
    Additional arguments can be passed for the plot. User can override default values. For example main="" will remove the title from the plot.
  }
}
\value{
\item{value}{A vector containing the importance value of each series.}
\item{class}{A vector containing the class membership of each series.}
\item{rank}{A vector containing the rank of each series. Highest value is 1.}
\item{conc}{The importance concentration of each class, as percentage of total value.}
}
\author{
Nikolaos Kourentzes
}
\seealso{
\code{\link{xyz}}, \code{\link{abcxyz}}.
}
\examples{
x <- abs(matrix(cumsum(rnorm(5400,0,1)),36,150))
abc(x,outplot=TRUE)
}

