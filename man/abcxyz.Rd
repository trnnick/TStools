\name{abcxyz}
\alias{abcxyz}

\title{ABC-XYZ visualisation}

\description{Visualise the results of both ABC-XYZ analyses and provide summary.}

\usage{
abcxyz(imp,frc,outplot=c(TRUE,FALSE),error=NULL,...)
}

\arguments{
  \item{imp}{
    Output of abc. 
  }
  \item{frc}{
    Output of xyz.
  }
  \item{outplot}{
    If TRUE provide a visualisation of the ABC analysis result.
  }
  \item{error}{
    Vector of errors for each item that will be averaged in each class.
  }
  \item{...}{
    Additional arguments can be passed for the plot. User can override default values. For example main="" will remove the title from the plot.
  }
}
\value{
\item{class}{Matrix containing the number of time series in each class.}
\item{error}{Matrix containing the average error for each class.}
}
\author{
Nikolaos Kourentzes
}
\seealso{
\code{\link{abc}}, \code{\link{xyz}}.
}
\examples{
x <- abs(matrix(cumsum(rnorm(5400,0,1)),36,150))
abcxyz(abc(x,outplot=FALSE),xyz(x,type="naive",outplot=FALSE))
}

