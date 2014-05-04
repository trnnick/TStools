\name{abcxyz}
\alias{abcxyz}

\title{ABC-XYZ visualisation}

\description{Visualise the results of both ABC-XYZ analysis and provide summary.}

\usage{
abcxyz(imp,frc,outplot=c(TRUE,FALSE))
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
}
\value{
\item{abcxyz}{Matrix containing the number of time series in each class.}
}
\author{
Nikolaos Kourentzes
}
\examples{
x <- abs(matrix(cumsum(rnorm(5400,0,1)),36,150))
abcxyz(abc(x,outplot=FALSE),xyz(x,type="cv",outplot=FALSE))
}

