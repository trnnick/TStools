\name{residout}
\alias{residout}

\title{Create a control chart of residuals and identify outliers}

\description{Create a control chart of residuals and identify outliers.}

\usage{
residout(resid,t=2,outplot=c(TRUE,FALSE))
}

\arguments{
  \item{resid}{
    Vector of residuals.
}
  \item{t}{
    Threshold value over which standardised residuals are regarded as outliers.
}
  \item{outplot}{
    If TRUE then a control chart of the standardised residuals is plotted.
}
}
\value{
\item{location}{Location of outliers.}
\item{outliers}{Value of outliers.}
\item{residuals}{Standardised residuals.}
}
\author{
Nikolaos Kourentzes
}
\examples{
resid <- rnorm(50)
residout(resid, outplot=TRUE)
}